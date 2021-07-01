__precompile__()
module LanczosInterpolation_mod

#using Yeppp 

export get_interp_lag,lanczos_interp_weight!, lanczos_interp_trace!

const NLAG=5

function sinc(x::Float64)::Float64
if abs(x) > 0
   return sin(pi*x)/(pi*x)
else
   return 1.0
end
end # sinc

function dsinc(x::Float64)::Float64
if abs(x) > 0
   return (cos(pi*x)-sinc(pi*x))/x
else
   return 0.0
end
end # dsinc

function lanczos(x::Float64)::Float64
   a::Float64=Float64(NLAG)
   if x < a
      return sinc(x)*sinc(x/a)
   else
      return 0.0
   end
end # lanczos

function get_interp_lag()
lag::Int64 = NLAG
return lag
end

function lanczos_interp_weight!( t::Float64,
                                 n::Int64,
                                 t0::Float64,
                                 dt::Float64,
                                 it::Int64,
                                 wlanczos::Array{Float64,1} )
        if t0 <= t &&  t <= t0+Float64(n-1)*dt
           x::Float64 = (t-t0)/dt
           it = floor(Int64,x)
           dx::Float64 = x-Float64(it)
           it += 1 
           for lag=-NLAG:NLAG
                   wlanczos[lag+NLAG+1]=lanczos(dx+Float64(lag))
           end
        else
           it = 0
           wlanczos[1:2*NLAG+1] .= Float64(0.0)
        end
nothing
end # lanczos_interp_weight!

function lanczos_interp_trace!( t::Float64,
                                n::Int64,
                                t0::Float64,
                                dt::Float64,
                                it::Int64,
                                trace::Array{Float64,1} )
        interp::Float64=0.0
        if t0 <= t && t <= t0+Float64(n-1)*dt
           x::Float64 = (t-t0)/dt
           it = floor(Int64,x)
           dx::Float64 = x-Float64(it)
           it += 1 
           for lag=-NLAG:NLAG
               is::Int64=Int64(min(max(it+lag,1),n))
               interp += trace[is]*lanczos(dx+Float64(lag))
           end
           return interp
        end
        return Float64(0.0)
end # lanczos_interp_trace!

end

