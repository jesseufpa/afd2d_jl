__precompile__()
module InvSmooth2d_mod
using LinearAlgebra
export Smooth2d! 

include("lsq_cgsolver.jl")

using  .LSQ_CGSolver

mutable struct CGsmooth2d 
#
# CG Data Structure:
n1::Int64
n2::Int64
nd::Int64
dat::Array{Float64,1}
scl::Float64
# cg precond
nlag::Int64
pcd::Float64
#
# Constructor:
function CGsmooth2d(img::Array{Float64,2})
n1::Int64=size(img,1)
n2::Int64=size(img,2)
nd::Int64 = n1*n2
dat::Array{Float64,1}=zeros(Float64,nd)
#
id::Int64=0
for i2=1:n2
    for i1=1:n1
        id+=1
        dat[id]=img[i1,i2]
    end
end
#
scl::Float64=sqrt(Float64(nd))
pnlag::Int64=5
pcd::Float64=0.0
#
@inbounds for i2::Int64 in -pnlag:pnlag
@inbounds for i1::Int64 in -pnlag:pnlag
    pcd+=1.0
end
end
pcd=1.0/pcd
#
new(n1,n2,nd,dat,scl,pnlag,pcd)
#
end
end
#
#traits
has_dataweight(cg::CGsmooth2d)=true
has_operator(cg::CGsmooth2d)=true
has_preconditioning(cg::CGsmooth2d)=true
has_modelweight(cg::CGsmooth2d)=true
has_cg_operators(cg::CGsmooth2d)=true
#
# Identity operator 
function oper!(adj::Int64,
               add::Int64,
               nm::Int64,
               mdl::Array{Float64,1},
               nd::Int64,
               dat::Array{Float64,1},
               cg::CGsmooth2d)
iexit::Int64=0
nd::Int64
if adj == 0
   if add == 0
      @inbounds dat[1:nd].=0.0
   end
@simd   for id in 1:nd
       @inbounds dat[id] = mdl[id]
   end
else
   if add == 0
      @inbounds mdl[1:nd].=0.0
   end
@simd   for id  in 1:nd
       @inbounds mdl[id] = dat[id]
   end
end
return iexit
end
#
# Box preconditioning
function pcnd!(adj::Int64,
               add::Int64,
               nm::Int64,
               mdl::Array{Float64,1},
               nd::Int64,
               dat::Array{Float64,1},
               cg::CGsmooth2d)

ic::Int64=0
iexit::Int64=0
if adj == 0
   #
   if add == 0
      @inbounds dat[1:nd] .= 0.0
   end
   #
   #
   for i2=cg.nlag+1:cg.n2-cg.nlag
   for i1=cg.nlag+1:cg.n1-cg.nlag
       #
       ic=i1+(i2-1)*cg.n1
       #
       for l2::Int64 in -cg.nlag:cg.nlag
           j2::Int64=i2+l2
           for l1::Int64 in -cg.nlag:cg.nlag
               j1::Int64=i1+l1
               if 0 < j2 && j2 <= cg.n2
               if 0 < j1 && j1 <= cg.n1
               idx::Int64=j1+(j2-1)*cg.n1
               @inbounds dat[ic] += cg.pcd*mdl[idx]
               end
               end
           end
       end
   end
   end
#
else
#
   if add == 0
      @inbounds mdl[1:nm].=0.0
   end
   #
   for i2=cg.nlag+1:cg.n2-cg.nlag
   for i1=cg.nlag+1:cg.n1-cg.nlag
       #
       ic=i1+(i2-1)*cg.n1
       #
       for l2::Int64 in -cg.nlag:cg.nlag
           j2::Int64=i2+l2
           for l1::Int64 in -cg.nlag:cg.nlag
               j1::Int64=i1+l1 
               if 0 < j2 && j2 <= cg.n2
               if 0 < j1 && j1 <= cg.n1
                  idx::Int64=j1+(j2-1)*cg.n1
                  @inbounds mdl[idx] += cg.pcd*dat[ic]
               end
               end
           end
       end
       #
   end
   end
   #
end
return iexit
end

# data weight 
wdat!(adj::Int64,
           nd::Int64,
           dat::Array{Float64,1},
           cg::CGsmooth2d)=0

# model selection weight 
wmdl!(adj::Int64,
      nc::Int64,
      mdl::Array{Float64,1},
      cg::CGsmooth2d)=0



function Smooth2d!(img::Array{Float64,2},
                   imgs::Array{Float64,2},
                   lmbd::Float64)
#
n1::Int64=size(img,1)
n2::Int64=size(img,2)
#
@assert size(imgs,1)==n1 && size(imgs,2)==n2
#
cg::CGsmooth2d=CGsmooth2d(img)
mdl::Array{Float64,1}=zeros(Float64,n1*n2)
niter::Int64=n1*n2
eps::Float64=0.000000000001
#
iexit::Int64=lsq_cgd(cg,
                     InvSmooth2d_mod.wdat!,
                     InvSmooth2d_mod.oper!,
                     InvSmooth2d_mod.wmdl!,
                     InvSmooth2d_mod.pcnd!,
                     cg.nd,mdl,
                     cg.nd,cg.dat,
                     niter,
                     lmbd*cg.scl,
                     eps)
#

nb::Int64=3*cg.nlag+1
id::Int64=0
for i2=nb+1:n2-nb
    for i1=nb+1:n1-nb
        id=i1+(i2-1)*n1
        imgs[i1,i2]=mdl[id]
    end
end
#=
 Boundary conditions:
=#

for i2=1:nb
   for i1=nb+1:n1-nb
       imgs[i1,i2]=imgs[i1,nb+1]
   end
end
for i2=n2-nb+1:n2
   for i1=nb+1:n1-nb
       imgs[i1,i2]=imgs[i1,n2-nb]
   end
end

for i2=1:n2
    for i1=1:nb
        imgs[i1,i2]      =imgs[nb+1,i2]
        imgs[n1-nb+i1,i2]=imgs[n1-nb,i2]
    end
end

scl::Float64=0.0
for i2=1:n2
    for i1=1:n1
        scl+= imgs[i1,i2]^2
    end
end


if scl > 0.0
   scl=sqrt(dot(cg.dat,cg.dat))/sqrt(scl)
end

for i2=1:n2
    for i1=1:n1
        imgs[i1,i2]=scl*imgs[i1,i2]
    end
end
return nothing

end


end #module
