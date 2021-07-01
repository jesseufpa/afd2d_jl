      function percentile(perc::Float64,p::Array{Float64,2})
      
      sig::Float64=0.0
      sig0::Float64=0.0
      sig1::Float64=0.0
      dsig::Float64=0.0
      pmin::Float64=0.0
      delta::Float64=0.0
      nvec::Int64= size(p,1)*size(p,2)
      nsig::Int64= Int64(floor(nvec/2)) 
      itnwt::Int64=0

      phi::Float64=0.0
      phi0::Float64=0.0
      phi1::Float64=0.0
      dphi::Float64=0.0
      TOL::Float64 = 0.000001
      #
      pmin=p[1,1]
      delta=pmin
      for i2=1:size(p,2)
          for i1=1:size(p,1)
              delta=max(delta,p[i1,i2])
              pmin=min(pmin,p[i1,i2])
          end
      end
      delta=delta-pmin
      if delta <= 0.0 
          println("Null interval ")
          return nothing
      end

      sig0=0.0
      phi0=0.0
      for i2=1:size(p,2)
          for i1=1:size(p,1)
              if (p[i1,i2]-pmin) >= sig0 
                 phi0 += 1.0
              end       
          end
      end
      phi0 = perc - phi0/Float64(nvec)
      sig1=delta
      phi1=0.0
      for i2=1:size(p,2)
          for i1=1:size(p,1)
              if (p[i1,i2]-pmin) >= sig1 
                 phi1+=1.0
              end       
          end
      end
      
      phi1=perc - phi1/Float64(nvec)
      dsig=-phi0*(sig1-sig0)/(phi1-phi0)
      sig = sig0+dsig

      for itnwt=1:99
      phi=0.0
      for i2=1:size(p,2)
          for i1=1:size(p,1)
              if (p[i1,i2]-pmin) >= sig 
                 phi=phi+1.0
              end       
          end
      end
      phi = perc - phi/Float64(nvec)

      if phi*phi0 < 0.0
         sig1 = sig
         phi1 = phi
      else 
         sig0 = sig
         phi0 = phi
      end
      dsig=-phi0*(sig1-sig0)/(phi1-phi0)
      sig = sig0+dsig

      if  abs(phi) < TOL || abs(dsig) < 1.e-16*delta
          println("iter = ",itnwt)
          break
      end 
      end

      sig=sig+pmin
      println("sig = ",sig, " phi = ",phi, " dsig = ",dsig)

      return nothing 
      end # function
       

