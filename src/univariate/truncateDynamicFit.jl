function truncateDynamicFit(μ::Float64, σ::Float64, range; kernel = Normal(),ϵ = 10e-5)
  a,b = range[1],range[2]
  if a == -Inf
        a = μ + σ*quantile(kernel,10e-5)
    end
  if b == Inf
        b = μ + σ*quantile(kernel,1-10e-5)
    end
  m = [0, 1, μ, σ^2+μ^2]
  F(x,θ1,θ2) = cdf(kernel,(x-θ1)/θ2)
  f(x,θ1,θ2) = (1/θ2)*pdf(kernel,(x-θ1)/θ2)
  ll(z) = a-(1-z)/z
  hh(z) = b+(1-z)/z

  function rhs(du,u,parms,z)
      n0 = F(hh(z),u[1],u[2]) - F(ll(z),u[1],u[2])
      l(i) = ll(z)^i*f(ll(z),u[1],u[2])
      h(i) = hh(z)^i*f(hh(z),u[1],u[2])
      p(i) = u[2]^2*(h(i)-l(i) - i*m[i+1]*n0)
      c(i) = h(i) + l(i) - m[i+2]*(h(0)+l(0))
      p0,p1,p2,p3 = p(0),p(1),p(2),p(3)
      c1,c2 = c(1),c(2)
      m1, m2 = m[3], m[4]
      ddd = (p2*(p1*m1+p0*m2)-p1^2*m2-p0*p3*m1+p3*p1-p2^2)/u[2]^5
      CC = (1/z^2)*(1/(u[2]^3*ddd))
      L1 = c2*(p0*m1*u[1]-p1*(m1+u[1])+p2)-c1*((p0*m2-p2)*u[1]-p1*m2+p3)
      L2 = u[2]*(c2*(p0*m1-p1)-c1*(p0*m2-p2))
      du[1] = CC*L1
      du[2] = CC*L2
   end

   u0 = [μ,σ/std(kernel)]
   tspan = (ϵ,1)

   ode = ODEProblem(rhs,u0,tspan)
   sol = solve(ode,Tsit5(),reltol=1e-8,abstol=1e-8)

   return (last(sol.u)[1],last(sol.u)[2])
end
