__precompile__()
@doc raw"""
    AFD2D_mod

# Finite-diference modeling of acoustic wave equation

```math
\frac{1}{c^2(\mathbf{x}} \frac{\partial P}{\partial t} = \nabla P + \dot{q}(t)
\delta(\mathbf{x} - \mathbf{x}_{s})
```
"""
module AFD2d_mod

#using Yeppp
#using Plots
using VegaLite
using DataFrames
using BenchmarkTools
using DelimitedFiles
using Printf
using Documenter
using StaticArrays

include("SU_mod.jl")
include("KaiserWindow_mod.jl")
include("LanczosInterpolation_mod.jl")
include("InvSmooth2d_mod.jl")
push!(LOAD_PATH, @__DIR__)
#push!(LOAD_PATH, pwd())

#import SU_mod
#import KaiserWindow_mod
#import LanczosInterpolation_mod
#import InvSmooth2d_mod

using .SU_mod
using .KaiserWindow_mod
using .LanczosInterpolation_mod
using .InvSmooth2d_mod
#using Profile

export AFD2d, Model2d, afd_mod

const Rabc=1.e-3
const nabc=49
const CFL=0.5
const ndrv1=2
const ndrv2=7
const dtsampling=0.004

const drv1=[ 9.0/8.0, 
            -1.0/24.0]
#
const drv2=[-3.171097979234166, 
      1.883703563645937,
     -0.392013173326410,
      0.126419594523541,
     -0.043436978346425,
      0.013633723622292,
     -0.002757740501852]
#
const KAISERWIND=4
#
@doc raw"""
    Grid2d(n1::Int64,
           n2::Int64,
           o1::Float64,
           o2::Float64,
           d1::Float64
           d2::Float64)

    Type to define a regular grid for Finite-Difference

##    Constructor Signature:

* n1: number of grid point along dimension 1
* n2: number of grid point along dimension 2
* o1: x1-coordinate of grid origin
* o2: x2-coordinate of grid origin
* d1: sampling interval along x1-coordinate
* d2: sampling interval along x2-coordinate
"""
mutable struct Grid2d
n1::Int64
n2::Int64
o1::Float64
o2::Float64
d1::Float64
d2::Float64
end
#

mutable struct Model2d
grd::Grid2d
n1fd::Int64
n2fd::Int64
vel::Array{Float64,2}
@doc raw""" 
    Model2d(file::String)
   
    Type Model2d stores the 2D velocity model and
    extends its grid boundaries for FD computations

## Constructor signature:

* file: velocity model header file
"""
function Model2d(file::String)
        hdr=readdlm(file,' ')
        n1::Int64=hdr[1,1]
        n2::Int64=hdr[1,2]
        o1::Float64=hdr[1,3]
        o2::Float64=hdr[1,4]
        d1::Float64=hdr[1,5]
        d2::Float64=hdr[1,6]
        binfile::String=hdr[2,1]

        grd::Grid2d=Grid2d(n1,n2,o1,o2,d1,d2)

        n1fd::Int64=n1+2*nabc
        n2fd::Int64=n2+2*nabc
        vel=zeros(Float64,n1fd,n2fd)
        value::Float32=Float32(0.0)
        iob::IOStream=open(binfile,"r")
        @assert isopen(iob)
        for i2=1:n2
            for i1=1:n1
                    vel[nabc+i1,nabc+i2]=Float64(read(iob,Float32))
            end
        end
        close(iob)

# left boundary
        for i2=1:nabc
        for i1=nabc+1:nabc+n1
            vel[i1,i2]=vel[i1,nabc+1]
        end 
        end
# right boundary
        for i2 in nabc+n2+1:n2fd
        for i1 in nabc+1:nabc+n1
            vel[i1,i2]=vel[i1,nabc+n2]
        end 
        end

# top boundary
        for i2 in 1:n2fd
        for i1 in 1:nabc
            vel[i1,i2]=vel[nabc+1,i2]
        end 
# bottom boundary
        for i1 in nabc+n1+1:n1fd
            vel[i1,i2]=vel[nabc+n1,i2]
        end 
        end 

        new(grd,n1fd,n2fd,vel)
end

function Model2d(n1::Int64,
                 n2::Int64,
                 o1::Float64,
                 o2::Float64,
                 d1::Float64,
                 d2::Float64
                )
grd::Grid2d=Grid2d(n1,n2,o1,o2,d1,d2)
n1fd::Int64=n1+2*nabc
n2fd::Int64=n2+2*nabc
vel=zeros(Float64,n1fd,n2fd)
new(grd,n1fd,n2fd,vel)
end
end


mutable struct AFD2dGeom
grd::Grid2d
n::Int64
i1::Array{Int64,1}
i2::Array{Int64,1}
wind1::Array{Float64,2}
wind2::Array{Float64,2}

function AFD2dGeom(grd::Grid2d,
                  x1stat::Array{Float64,1},x2stat::Array{Float64,1})
@assert length(x1stat) == length(x2stat) 
# Kaiser window parameters
bKaiser::Float64=4.14
lKaiser::Float64=Float64(KAISERWIND)

nstat::Int64 = length(x1stat)
i1::Array{Int64,1}=zeros(Int64,nstat)
i2::Array{Int64,1}=zeros(Int64,nstat)
wind1::Array{Float64,2}=zeros(Float64,2*KAISERWIND+1,nstat)
wind2::Array{Float64,2}=zeros(Float64,2*KAISERWIND+1,nstat)
for is=1:nstat
    s1::Float64=(x1stat[is]-grd.o1)/grd.d1+Float64(nabc)
    s2::Float64=(x2stat[is]-grd.o2)/grd.d2+Float64(nabc)
    @assert s1 > 0.0 && s2 > 0.0
    i1[is]=floor(Int64,s1)
    i2[is]=floor(Int64,s2)

#=
Kaiser window coefficients
for each station
=#
    dx1::Float64=s1-Float64(i1[is])
    @simd for i=-KAISERWIND:KAISERWIND
          @inbounds  wind1[i+KAISERWIND+1,is]=KaiserWindow(Float64(i)+dx1,lKaiser,bKaiser)
    end
    dx2::Float64=s2-Float64(i2[is])
    @simd for i=-KAISERWIND:KAISERWIND
           @inbounds wind2[i+KAISERWIND+1,is]=KaiserWindow(Float64(i)+dx2,lKaiser,bKaiser)
    end
#= normalization of Kaiser window =#
    wsum::Float64=0.0
    for j=-KAISERWIND:KAISERWIND
    for i=-KAISERWIND:KAISERWIND
       @inbounds wsum += wind1[i+KAISERWIND+1,is]*wind2[j+KAISERWIND+1,is]
    end
    end
    wsum=1.0/sqrt(wsum)
    @inbounds wind1[:,is] *= wsum
    @inbounds wind2[:,is] *= wsum
end
new(grd,nstat,i1,i2,wind1,wind2)
end
end


# data structure for acoustic finite difference
mutable struct AFD2d
mdl::Model2d
freq::Float64
dt::Float64
h1::Float64
h2::Float64
abc1::Array{Float64,2}
abc2::Array{Float64,2}
abch1::Array{Float64,2}
abch2::Array{Float64,2}

function AFD2d(m::Model2d,freq::Float64)
vmin::Float64=m.vel[1,1]
vmax::Float64=m.vel[1,1]
@inbounds for i2 in 1:m.n2fd
@inbounds for i1 in 1:m.n1fd
    vmin=min(vmin,m.vel[i1,i2])
    vmax=max(vmax,m.vel[i1,i2])
end
end
# check dispersion
d::Float64=vmin/(6.0*freq)

if  min(m.grd.d1,m.grd.d2) > d
    println("vmin = ",vmin)
    println("d    = ",d)
    println("d1   = ",m.grd.d1)
    exit()
end
# check stability
arg::Float64=pi/sqrt(2.0)
cfl::Float64=CFL 
sum::Float64=drv2[1]
for ic in 2:ndrv2
   @inbounds sum += 2.0*drv2[ic]*cos(Float64(ic)*arg)
end
mu::Float64 = cfl*sqrt(2.0/abs(sum))
dt::Float64 = Float64(floor(Int64,10000.0*mu*min(m.grd.d1,m.grd.d2)/vmax)/10000)

println("dt = ",dt)
ndt = floor(Int64,dtsampling/dt)
if  ndt > 1 
    dt= dtsampling/Float64(ndt)
    ndt+=1
else
    dt=dtsampling
    ndt=1
end

println("d1 = ",m.grd.d1)
println("d2 = ",m.grd.d2)
println("dt = ",dt)
println("vmin = ",vmin)
println("vmax = ",vmax)
# fill up structure

mu = (dt*dt) #/(m.grd.d1*m.grd.d2)
@simd  for i2 in 1:m.n2fd
@simd for i1 in 1:m.n1fd
        @inbounds m.vel[i1,i2] *= mu*m.vel[i1,i2]
end
end


h1::Float64=1.0/m.grd.d1
h2::Float64=1.0/m.grd.d2

# CPML absorbing boundary conditions

abc1::Array{Float64,2}=zeros(Float64,2,m.n1fd)
abc2::Array{Float64,2}=zeros(Float64,2,m.n2fd)
abch1::Array{Float64,2}=zeros(Float64,2,m.n1fd)
abch2::Array{Float64,2}=zeros(Float64,2,m.n2fd)

d01::Float64 = -1.5*log(Rabc)*vmax/(Float64(nabc)*m.grd.d1)
for i in 1:nabc
    uabc::Float64=Float64(i)/Float64(nabc)
    alpha::Float64=pi*freq*(1.0-uabc)
    d1::Float64 = d01*(uabc-sin(2.0*pi*uabc)/(2.0*pi))
    aux1::Float64 = alpha+d1
    aux2::Float64 = exp(-aux1*dt)
    abc1[1,nabc+1-i] = aux2
    abc1[2,nabc+1-i] = -d1*(1.0-aux2)/aux1
    abc1[1,m.n1fd-nabc+i] = aux2
    abc1[2,m.n1fd-nabc+i] = -d1*(1.0-aux2)/aux1
end

for i in 1:nabc
    uabc::Float64=(Float64(i-1)+0.5)/Float64(nabc)
    alpha::Float64=pi*freq*(1.0-uabc)
    d1::Float64 = d01*(uabc-sin(2.0*pi*uabc)/(2.0*pi))
    aux1::Float64 = alpha+d1
    aux2::Float64 = exp(-aux1*dt)
    abch1[1,nabc+1-i] = aux2
    abch1[2,nabc+1-i] = -d1*(1.0-aux2)/aux1

    uabc =(Float64(i)+0.5)/Float64(nabc)
    alpha=pi*freq*(1.0-uabc)
    d1   = d01*(uabc-sin(2.0*pi*uabc)/(2.0*pi))
    aux1 = alpha+d1
    aux2 = exp(-aux1*dt)
    abch1[1,m.n1fd-nabc+i] = aux2
    abch1[2,m.n1fd-nabc+i] = -d1*(1.0-aux2)/aux1
end

d02::Float64 = -1.5*log(Rabc)*vmax/(Float64(nabc)*m.grd.d2)
for i in 1:nabc
    uabc::Float64=Float64(i)/Float64(nabc)
    alpha::Float64=pi*freq*(1.0-uabc)
    d2::Float64 = d02*(uabc-sin(2.0*pi*uabc)/(2.0*pi))
    aux1::Float64 = alpha+d2
    aux2::Float64 = exp(-aux1*dt)
    abc2[1,nabc+1-i] = aux2
    abc2[2,nabc+1-i] = -d2*(1.0-aux2)/aux1
    abc2[1,m.n2fd-nabc+i] = aux2
    abc2[2,m.n2fd-nabc+i] = -d2*(1.0-aux2)/aux1
end

for i in 1:nabc
    uabc::Float64=(Float64(i-1)+0.5)/Float64(nabc)
    alpha::Float64=pi*freq*(1.0-uabc)
    d2::Float64 = d02*(uabc-sin(2.0*pi*uabc)/(2.0*pi))
    aux1::Float64 = alpha+d2
    aux2::Float64 = exp(-aux1*dt)
    abch2[1,nabc+1-i] = aux2
    abch2[2,nabc+1-i] = -d2*(1.0-aux2)/aux1

    uabc=(Float64(i)+0.5)/Float64(nabc)
    alpha=pi*freq*(1.0-uabc)
    d2   = d02*(uabc-sin(2.0*pi*uabc)/(2.0*pi))
    aux1 = alpha+d2
    aux2 = exp(-aux1*dt)
    abch2[1,m.n2fd-nabc+i] = aux2
    abch2[2,m.n2fd-nabc+i] = -d2*(1.0-aux2)/aux1
end

println("vmin = ",vmin)
println("vmax = ",vmax)
new(m,freq,dt,h1,h2,abc1,abc2,abch1,abch2)
end

end

function afd_ricker(FD::AFD2d,t::Float64)
arg::Float64=(pi*FD.freq*(t-1.0/FD.freq))^2
pulse::Float64=exp(-arg)*(1.0-2.0*arg)
pulse
end
 
function afd_update(FD::AFD2d,p1::Array{Float64,2},
                              p2::Array{Float64,2},
                              phi1::Array{Float64,2},
                              phi2::Array{Float64,2},
                              psi1::Array{Float64,2},
                              psi2::Array{Float64,2},
        )
#=

 Finite difference 
 time update

=#
    Threads.@threads  for i2 in ndrv2:FD.mdl.n2fd-ndrv2
        @inbounds for i1 in ndrv2:FD.mdl.n1fd-ndrv2
        dp1::Float64 =0.0
        dp2::Float64 =0.0
        for ic in 0:ndrv1-1
                dp1 += drv1[ic+1]*(p2[i1+ic+1,i2]-p2[i1-ic,i2])
        end
        dp1 *= FD.h1
        phi1[i1,i2]=FD.abch1[1,i1]*phi1[i1,i2]+FD.abch1[2,i1]*dp1

        for ic in 0:ndrv1-1
                dp2 += drv1[ic+1]*(p2[i1,i2+ic+1]-p2[i1,i2-ic])
        end
        dp2 *= FD.h2
        phi2[i1,i2]=FD.abch2[1,i2]* phi2[i1,i2]+FD.abch2[2,i2]*dp2
    end
    end


    Threads.@threads  for i2 in ndrv2:FD.mdl.n2fd-ndrv2
        @inbounds for i1 in ndrv2:FD.mdl.n1fd-ndrv2

        dphi1::Float64 =0.0
        dphi2::Float64 =0.0
        for ic in 0:ndrv1-1
                dphi1 += drv1[ic+1]*(phi1[i1+ic,i2]-phi1[i1-ic-1,i2])
        end
        dphi1 *= FD.h1

        for ic in 0:ndrv1-1
                dphi2 += drv1[ic+1]*(phi2[i1,i2+ic]-phi2[i1,i2-ic-1])
        end
        dphi2 *= FD.h2

        laplacian1::Float64=drv2[1]*p2[i1,i2]
        laplacian2::Float64=drv2[1]*p2[i1,i2]
        for ic in 1:ndrv2-1
             @inbounds laplacian1 += drv2[ic+1]*(p2[i1-ic,i2]+p2[i1+ic,i2])
        end
        laplacian1 *= FD.h1*FD.h1
        for ic in 1:ndrv2-1
             @inbounds laplacian2 += drv2[ic+1]*(p2[i1,i2-ic]+p2[i1,i2+ic])
        end
        laplacian2 *= FD.h2*FD.h2
        psi1[i1,i2]=FD.abc1[1,i1]* psi1[i1,i2]+FD.abc1[2,i1]*(laplacian1+dphi1)
        psi2[i1,i2]=FD.abc2[1,i2]* psi2[i1,i2]+FD.abc2[2,i2]*(laplacian2+dphi2)

        laplacian::Float64=laplacian1+laplacian2+dphi1+dphi2+psi1[i1,i2]+psi2[i1,i2]
        p1[i1,i2]=2.0*p2[i1,i2]-p1[i1,i2]+FD.mdl.vel[i1,i2]*laplacian
    end
    end

    Threads.@threads for i2 in 1:FD.mdl.n2fd
    @inbounds for i1 in 1:FD.mdl.n1fd
        swap::Float64=p1[i1,i2]
        p1[i1,i2]=p2[i1,i2]
        p2[i1,i2]=swap
    end
    end

    # Dirichlet boundary to avoid stability
    Threads.@threads for i1 in 1:FD.mdl.n1fd
            p2[i1,ndrv2]=0.0
            p2[i1,FD.mdl.n2fd-ndrv2]=0.0
    end

    Threads.@threads for i2 in 1:FD.mdl.n2fd
            p2[ndrv2,i2]=0.0
            p2[FD.mdl.n1fd-ndrv2,i2]=0.0
    end

end 

function afd_src_injection(FD::AFD2d,Src::AFD2dGeom,t::Float64,p1::Array{Float64,2})
@inbounds  Threads.@threads  for is=1:Src.n
@inbounds        for iw2=-KAISERWIND:KAISERWIND
                     i2s::Int64=min(max(Src.i2[is]+iw2,1),FD.mdl.n2fd)
                     wsrc2::Float64=Src.wind2[iw2+KAISERWIND+1]
@inbounds        for iw1=-KAISERWIND:KAISERWIND
                     i1s::Int64=min(max(Src.i1[is]+iw1,1),FD.mdl.n1fd)
                     wsrc::Float64=Src.wind1[iw1+KAISERWIND+1]*wsrc2

                    p1[i1s,i2s] -= wsrc*FD.mdl.vel[i1s,i2s]*(afd_ricker(FD,t+FD.dt)+
                         afd_ricker(FD,t-FD.dt)-
                         2.0*afd_ricker(FD,t))/(FD.dt*FD.dt)
        end
        end
    end 
    nothing
end

function afd_rcv_field(FD::AFD2d,Rcv::AFD2dGeom,is::Int64,p2::Array{Float64,2})
                 field::Float64 = 0.0
@inbounds        for iw2=-KAISERWIND:KAISERWIND
                     i2s::Int64=min(max(Rcv.i2[is]+iw2,1),FD.mdl.n2fd)
                     wrcv2::Float64=Rcv.wind2[iw2+KAISERWIND+1]
@inbounds        for iw1=-KAISERWIND:KAISERWIND
                     i1s::Int64=min(max(Rcv.i1[is]+iw1,1),FD.mdl.n1fd)
                     wrcv::Float64=Rcv.wind1[iw1+KAISERWIND+1]*wrcv2
                     field += wrcv*p2[i1s,i2s]
        end
        end
        return field 
end



@doc raw"""
    afd_mod(FD::AFD2d,
            T::Float64,
            x1s::Array{Float64,1},
            x2s::Array{Float64,1};
            x1r::Array{Float64,1},
            x2r::Array{Float64,1},
            section::Array{SU_trc,1}
            )

   
    Acoustic wavefield forward modeling

## Signature:

* FD : finite-Difference type
* T  : simulation time
* x1s: source stations x1-coordinates
* x2s: source stations x2-coordinates
* x1r: receiver stations x1-coordinates
* x2r: receiver stations x2-coordinates
"""
function afd_mod!(FD::AFD2d,adj::Bool,
                 T::Float64,
                 x1s::Array{Float64,1},
                 x2s::Array{Float64,1};
                 x1r::Array{Float64,1},
                 x2r::Array{Float64,1},
                 section::Array{SU_trc,1}
                )

nt::Int64  = floor(Int64,T/FD.dt+1.0)
ndt::Int64 = floor(Int64,dtsampling/FD.dt)+1


#= set source stations =#
Src::AFD2dGeom=AFD2dGeom(FD.mdl.grd,x1s,x2s)
#= set receiver stations =#
Rcv::AFD2dGeom=AFD2dGeom(FD.mdl.grd,x1r,x2r)

ns::Int16=floor(Int16,T/dtsampling+1.0)

scale::Float64=10.0
scalco::Int16=Int16(scale)
scalel::Int16=Int16(scale)
@inbounds @simd for ir=1:length(x1r)
    section[ir].hdr.scalco=-floor(Int16,scalco)
    section[ir].hdr.scalel=-floor(Int16,scalel)
    section[ir].hdr.tracl=Int32(ir)
    section[ir].hdr.ns=ns
    section[ir].hdr.dt=Int16(1000000.0*dtsampling)
    section[ir].hdr.sx=Int32(x2s[1]*scale)
    section[ir].hdr.sdepth=Int32(x1s[1]*scale)
    section[ir].hdr.gx=Int32(x2r[ir]*scale)
    section[ir].hdr.gelev=Int32(-x1r[ir]*scale)
    section[ir].hdr.ntr=Int32(length(x1r))
end

irec::Int64 = 0
trec::Float64=0.0
dtrec::Float64=Float64(dtsampling)

# state field allocation and initialization
p1::Array{Float64,2}=zeros(Float64,FD.mdl.n1fd,FD.mdl.n2fd)
p2::Array{Float64,2}=zeros(Float64,FD.mdl.n1fd,FD.mdl.n2fd)

# memory field for ABC allocation and initialization
phi1::Array{Float64,2}=zeros(Float64,FD.mdl.n1fd,FD.mdl.n2fd)
phi2::Array{Float64,2}=zeros(Float64,FD.mdl.n1fd,FD.mdl.n2fd)
psi1::Array{Float64,2}=zeros(Float64,FD.mdl.n1fd,FD.mdl.n2fd)
psi2::Array{Float64,2}=zeros(Float64,FD.mdl.n1fd,FD.mdl.n2fd)

zm::Array{Float64,1}=zeros(Float64,FD.mdl.n1fd)
xm::Array{Float64,1}=zeros(Float64,FD.mdl.n2fd)

@inbounds for i1=1:FD.mdl.n1fd
        zm[i1]=FD.mdl.grd.o1 + Float64(i1-1)*FD.mdl.grd.d1
end
@inbounds for i2=1:FD.mdl.n2fd
        xm[i2]=FD.mdl.grd.o2 + Float64(i2-1)*FD.mdl.grd.d2
end

ndtsampling::Int64 = floor(Int64,dtsampling/FD.dt)+1

isnap::Int64=0
tsnap::Float64 = 0.0
dtsnap::Float64 = 50*dtsampling

        iob=open("fdmovie.bin","w")
        @assert isopen(iob)

if  adj == 0 
@inbounds for it in 1:nt

    t::Float64=Float64(it-1)*FD.dt

    # source injection:

    afd_src_injection(FD,Src,t,p1)

    afd_update(FD,p1,p2,phi1,phi2,psi1,psi2)

    if abs(t-tsnap) < 0.5*FD.dt
       isnap += 1
       tsnap += dtsnap
       time::String=@sprintf("%.3f",Float64(isnap-1)*dtsnap)
       #=
      img=heatmap(xm,zm,FD.mdl.vel*maximum(p2)/maximum(FD.mdl.vel),
                  reuse=false,xlabel="x (m)",ylabel="y (m)",
                  yflip=true,colorbar=false)

      heatmap!(xm,zm,p2,xlabel="x (m)",ylabel="y
                       (m)",alpha=0.25,title="Snapshot at "*time*" s",yflip=true)
              #display(plot(img))
        =#

        for i2=1:size(p2,2)
            for i1=1:size(p2,1)
                write(iob,Float32(p2[i1,i2]))
            end
        end
        #=
        data=DataFrame(p2,:auto)
        select(data,[:x1,:x2]) |> @vlplot({:rect,color=data},x=:x1,y=:x2) |> display 
        =#
    end
    
    if  abs(t-trec) < FD.dt
        irec += 1
        trec += dtrec 
        @inbounds Threads.@threads for ir=1:length(x1r)
                section[ir].dat[irec] = Float32(afd_rcv_field(FD,Rcv,ir,p2))
        end
    end

end
        close(iob)

else
@inbounds for it in 1:nt

    t::Float64=Float64(it-1)*FD.dt

    # source injection:

    afd_src_injection(FD,Src,t,p1)

    afd_update(FD,p1,p2,phi1,phi2,psi1,psi2)

    if abs(t-tsnap) < 0.5*FD.dt
       isnap += 1
       tsnap += dtsnap
       time::String=@sprintf("%.3f",Float64(isnap-1)*dtsnap)
      img=heatmap(xm,zm,FD.mdl.vel*maximum(p2)/maximum(FD.mdl.vel),
                  reuse=false,xlabel="x (m)",ylabel="y (m)",
                  yflip=true,colorbar=false)

      heatmap!(xm,zm,p2,xlabel="x (m)",ylabel="y
                       (m)",alpha=0.25,title="Snapshot at "*time*" s",yflip=true)
           display(plot(img))

              for i2=1:size(p2,2)
                for i1=1:size(p2,1)
                write(iob,Float32(p2[i1,i2]))
            end
        end
    end
    
    if  abs(t-trec) < FD.dt
        irec += 1
        trec += dtrec 
        @inbounds Threads.@threads for ir=1:length(x1r)
                section[ir].dat[irec] = Float32(afd_rcv_field(FD,Rcv,ir,p2))
        end
    end

end
        close(iob)

end
end



#function __init__()
function Main()

freq::Float64=20.0
T::Float64=3.0
adj::Bool=0

#mdl::Model2d=Model2d("marm_ii_vp_10m.hdr")
#mdl::Model2d=Model2d("velstep.hdr")
mdl::Model2d=Model2d("velbd.hdr")
println("n1, n2 ",mdl.n1fd," ",mdl.n2fd)
#=
vels::Array{Float64,2}=zeros(mdl.n1fd,mdl.n2fd)
println("n1, n2 ",size(vels,1)," ",size(vels,2))

Smooth2d!(mdl.vel,vels,0.1)
        binfile="vels.bin"
        iob=open(binfile,"w")
        @assert isopen(iob)
        for i2=1:size(vels,2)
                for i1=1:size(vels,1)
                write(iob,Float32(vels[i1,i2]))
            end
        end
        close(iob)
=#

AFD::AFD2d=AFD2d(mdl,freq)
x1s::Array{Float64,1}=zeros(Float64,1)
x2s::Array{Float64,1}=zeros(Float64,1)
x1r::Array{Float64,1}=zeros(Float64,AFD.mdl.grd.n2)
x2r::Array{Float64,1}=zeros(Float64,AFD.mdl.grd.n2)

x1s[1]=200*AFD.mdl.grd.d1 + AFD.mdl.grd.o1
x2s[1]=Float64(AFD.mdl.grd.n2-1)*AFD.mdl.grd.d2/2.0+AFD.mdl.grd.o2

@simd for ir=1:length(x1r)
     @inbounds x1r[ir] = mdl.grd.o1 + 120*mdl.grd.d1 
end

@simd for ir=1:length(x2r)
    @inbounds x2r[ir] = mdl.grd.o2 + Float64(ir-1)*mdl.grd.d2 
end

ns::Int64=floor(Int64,T/dtsampling)+1
seis::Array{SU_trc,1}=Array{SU_trc,1}(undef,mdl.grd.n2)
Threads.@threads for itrc=1:length(x1r)
    @inbounds seis[itrc]=SU_trc(ns)
end

println("x1s = ",x1s[1])
println("x2s = ",x2s[1])

in::IOStream=open("survey.su","r")
try
    isopen(in)
catch 
   println("Error opening file survey.su file")
end
sutrace::SU_trc=SU_trc(ns)
SUGetTrace(in,sutrace)
close(in)

#=
#@assert isopen(io)
println("ntr = ",length(x1r))
println("ns  = ",ns)
@inbounds for itrc=1:length(x1r)
        seis[itrc].hdr.tracl=Int32(itrc)
        seis[itrc].hdr.ns=Int16(ns)
        seis[itrc].hdr.dt=Int16(2000)
        SUPutTrace(io,seis[itrc])
        #=
        @inbounds for is=1:ns
                write(io,Float32(seis[itrc].dat[is]))
        end
        =#
end

=#


@time afd_mod!(AFD,adj,T,x1s,x2s,x1r=x1r,x2r=x2r,section=seis)

io::IOStream=open("seis.su","w")
try
    isopen(io)
catch 
   println("Error opening file seis.su file")
end
#@assert isopen(io)
println("ntr = ",length(x1r))
println("ns  = ",ns)
@inbounds for itrc=1:length(x1r)
        seis[itrc].hdr.tracl=Int32(itrc)
        seis[itrc].hdr.ns=Int16(ns)
        seis[itrc].hdr.dt=Int16(2000)
        SUPutTrace(io,seis[itrc])
        #=
        @inbounds for is=1:ns
                write(io,Float32(seis[itrc].dat[is]))
        end
        =#
end
close(io)
nothing
end
end

