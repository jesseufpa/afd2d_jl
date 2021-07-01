module LSQ_CGSolver  

using LinearAlgebra

export lsq_cgd,lsq_cgm

"""
Generic implemention of least-squares 
conjugate gradient algorithm
"""

#=
function wdat! end
implements data weight police

function oper! end
implements the action of the linear operator 

function pcnd! end
implements preconditioning by model reparameterization

function wmdl! end
implements selection or model weight police for model parameters
=#

wdat!(args...) = 0
function oper! end
function pcnd! end
wmdl!(args...) = 0


#traits

Has_DataWeight(cgsolver)      = error("Not implemented")
Has_Operator(cgsolver)        = error("Not implemented")
Has_Precontioning(cgsolver)   = error("Not implemented")
Has_ModelWeight(cgsolver)     = error("Not implemented")
Has_CG_Operators(cgsolver)    = error("Not implemented")

"""

Algorithm for least-squares CG data type 

"""

function lsq_cgd(cgsolver,
                wdat!::Function,
                oper!::Function,
                wmdl!::Function,
                pcnd!::Function,
                nm::Int64,mdl::Array{Float64,1},
                nd::Int64,dat::Array{Float64,1},
                niter::Int64,
                lmbd::Float64,
                eps::Float64) 
#=

=#
p::Array{Float64,1}=zeros(Float64,nd)
q::Array{Float64,1}=zeros(Float64,nd)
r::Array{Float64,1}=zeros(Float64,nd)
x::Array{Float64,1}=zeros(Float64,nd)

sd::Array{Float64,1}=zeros(nd)
sm::Array{Float64,1}=zeros(nm)

@simd for id=1:nd
@inbounds sd[id] = dat[id]
end
iexit::Int64=wdat!(1,nd,sd,cgsolver)
@simd for id=1:nd
@inbounds sd[id] = lmbd * sd[id]
end

@simd for id=1:nd
@inbounds r[id] = sd[id]
end

@simd for id in 1:nd
@inbounds p[id] = r[id]
end

res::Float64=dot(r,r)
resmin::Float64=res
bnorm::Float64=res

if res <= eps * res
   iexit = 0
   return iexit
end 

iexit = -1
tol::Float64=eps*res

for iter::Int64 in 1:niter

iexit=wdat!(1,nd,p,cgsolver)
iexit=oper!(1,0,nm,mdl,nd,p,cgsolver)
iexit=wmdl!(1,nm,mdl,cgsolver)
iexit=pcnd!(1,0,nm,sm,nm,mdl,cgsolver)

@simd for id in 1:nd
@inbounds  sd[id] = lmbd * p[id]
end

iexit=pcnd!(0,0,nm,sm,nm,mdl,cgsolver)
iexit=wmdl!(0,nm,mdl,cgsolver)
iexit=oper!(0,0,nm,mdl,nd,q,cgsolver)
iexit=wdat!(0,nd,q,cgsolver)

@simd for id in 1:nd
@inbounds q[id] += lmbd * sd[id]
end

alpha::Float64=res/(dot(sd,sd)+dot(sm,sm))

@simd for id in 1:nd
      @inbounds r[id] = r[id] - alpha*q[id]
end

@simd for id in 1:nd
@inbounds x[id] += alpha * p[id]
end

res0::Float64 = res
res    = dot(r,r);
resmin = min(resmin,res)


dxnorm::Float64 = abs(alpha) * sqrt(dot(p,p))
xnorm::Float64  = sqrt(dot(x,x))

if dxnorm < eps * xnorm
   iexit = 2
   break
end
if abs(res) < tol 
   iexit = 0
   break
end


beta::Float64  =res/res0 
@simd for id in 1:nd
@inbounds  p[id] =r[id] + beta*p[id]
end

if mod(iter,10) == 0
   println("iter = ",iter," res = ",res," dxnorm = ",tol)
end 

end # iter
istop::Int64=iexit
   iexit=wdat!(1,nd,x,cgsolver)
   iexit=oper!(1,0,nm,mdl,nd,x,cgsolver)
   iexit=wmdl!(1,nm,mdl,cgsolver)
   iexit=pcnd!(1,0,nm,sm,nm,mdl,cgsolver)

   iexit=pcnd!(0,0,nm,sm,nm,mdl,cgsolver)
   iexit=wmdl!(0,nm,mdl,cgsolver)

if lmbd > 0.0
@simd for id=1:nm
@inbounds    mdl[id] =mdl[id]/lmbd
end
end

println("ISTOP = ",istop)

return iexit

end #cgls_d

"""

Algorithm for least-squares CG model type 

"""

function lsq_cgm(cgsolver,
                wdat!::Function,
                oper!::Function,
                wmdl!::Function,
                pcnd!::Function,
                nm::Int64,mdl::Array{Float64,1},
                nd::Int64,dat::Array{Float64,1},
                niter::Int64,
                lmbd::Float64,
                eps::Float64)
#=

=#
p::Array{Float64,1}=zeros(Float64,nm)
r::Array{Float64,1}=zeros(Float64,nm)
s::Array{Float64,1}=zeros(Float64,nm)
x::Array{Float64,1}=zeros(Float64,nm)

qd::Array{Float64,1}=zeros(nd)
qm::Array{Float64,1}=zeros(nm)

@simd for id in 1:nd
@inbounds qd[id] = dat[id]
end


iexit::Int64=wdat!(1,nd,qd,cgsolver)
iexit=oper!(1,0,nm,mdl,nd,qd,cgsolver)
iexit=wmdl!(1,nm,mdl,cgsolver)
iexit=pcnd!(1,0,nm,r,nm,mdl,cgsolver)

@simd for id in 1:nm
@inbounds p[id] = r[id]
end

iexit = -1
res::Float64=dot(r,r)
resmin::Float64=res

if res < eps * res
   iexit = 0
   return iexit
end 

iexit = -1
tol::Float64=eps*res

for iter in 1:niter


iexit=pcnd!(0,0,nm,p,nm,mdl,cgsolver)
iexit=wmdl!(0,nm,mdl,cgsolver)
iexit=oper!(0,0,nm,mdl,nd,qd,cgsolver)
iexit=wdat!(0,nd,qd,cgsolver)

@simd for id in 1:nm
@inbounds qm[id] = lmbd * p[id]
end

iexit=wdat!(1,nd,qd,cgsolver)
iexit=oper!(1,0,nm,mdl,nd,qd,cgsolver)
iexit=wmdl!(1,nm,mdl,cgsolver)
iexit=pcnd!(1,0,nm,s,nm,mdl,cgsolver)

@simd for id in 1:nm
    @inbounds s[id] = s[id] + lmbd * qm[id]
end

alpha::Float64=res/(dot(qd,qd)+dot(qm,qm))

@simd for id in 1:nm
        @inbounds r[id] = r[id]-alpha * s[id]
end

@simd for id in 1:nm
        @inbounds x[id] = x[id]+alpha * p[id]
end

res0::Float64   = res
res    = dot(r,r);
resmin = min(resmin,res)

dxnorm::Float64 = abs(alpha) * sqrt(dot(p,p))
xnorm::Float64  = sqrt(dot(x,x))
if dxnorm < eps * xnorm
   iexit = 2
   break
end

if res < tol 
   iexit = 0
   break
end

beta::Float64  =res/res0 
@simd for id in 1:nm
@inbounds    p[id] =r[id] + beta*p[id]
end


if mod(iter,10) == 0
   println("iter = ",iter," res = ",res," dxnorm = ",tol)
end 

end # iter
istop::Int64=iexit

iexit=pcnd!(0,0,nm,x,nm,mdl,cgsolver)
iexit=wmdl!(0,nm,mdl,cgsolver)
println("ISTOP = ",istop)
return istop

end #cgls_m

end #module

