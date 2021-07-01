__precompile__()
module KaiserWindow_mod
#using Yeppp
export KaiserWindow

const XCUT=3.75
const NP=7
const NQ=9
const p=[1.0000000,
         3.5156229,
         3.0899424,
         1.2067492,
         0.2659732,
         0.03607680,
         0.004581300]

const q=[0.39894228,
         0.13285920e-1,
         0.22531900e-2,
        -0.15756500e-2,
         0.91628100e-2,
        -0.20577060e-1,
         0.26355370e-1,
        -0.16476330e-1,
         0.39237700e-2]

function bessel_i0(x::Float64)::Float64
ax::Float64=abs(x)
if ax < XCUT
   y::Float64=ax/XCUT
   y *= y
   poly::Float64=p[NP]
   for n=NP-1:-1:1
       poly=p[n]+y*poly
    end
else
        y = XCUT / ax;
        poly=q[NQ]
        for n=NQ-1:-1:1
            poly = q[n]+y*poly
        end
        poly *= exp(ax)/sqrt(ax)
end
return poly
end

function KaiserWindow(x::Float64,lx::Float64,b::Float64)::Float64
    a::Float64=x/lx
    if a <= 1.0
       return bessel_i0(b * sqrt(1.0 - a * a)) / bessel_i0(b);
    else
       return Float64(0.0)
    end
end

end
