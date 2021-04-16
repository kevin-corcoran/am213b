module DiffyQ

using LinearAlgebra

export Newtons, CompTrapezoid, CompSimpson, Euler

function CompTrapezoid(N, a, b, f)

    #=
        Integrates f(t) from a to b using composite trapezoid rule

        Input variables
        N: number of points
        a: initial point
        b: final point
        f: function to be integrated

        local variables
        h: step size
        s: solution
    =#

    h = (b-a)/N
    s = 0
    for i = 1:N
        s += h/2 * (f(a + h*(i-1)) + f(a + h*i))
    end
    return s
end

function CompSimpson(N, a, b, f)

    #=
        Integrates f(t) from a to b using composite simpsons rule

        Input variables
        N: number of points
        a: initial point
        b: final point
        f: function to be integrated

        local variables
        h: step size
        s: solution
        ts: partition points
    =#

    h = (b-a)/N
    s = f(a) + f(b)

    ts = a:h:b
    for i = 1:(Int(N/2) - 1)
       s += 2*f(ts[2*i])
    end
    for i = 1:Int(N/2)
        s += 4*f(ts[2*i-1])
    end

    s = h/3 * s
    return s
end

function Newtons(f, df, x0, s, tol = 10.0^(-8.0), maxiter=10000)
    #=
        Uses Newtons method to find the roots of f(x) = 0
        Input:
            f(x,s): function of x with optional parameter s
            df(x,s): derivative of f with respect to x
            x0: initial guess
            s: parameter
            tol: desired accuracy
            maxiter: maximum number of iterations
    =#

    xn = x0
    for i = 1:maxiter
        xn_1 = xn - f(xn,s)/df(xn,s)
        # Stop once desired accuracy is achieved 
        if (abs(xn_1 - xn) < tol)
            # println("success")
            break
        end
        xn = xn_1
    end
    return xn
end

function Euler(func, N, T, t0, u0, λ=0)
    u = zeros(N+1)
    u[1] = u0
    h = T/N
    t = t0
    for i = 1:N
        u[i+1] = u[i] + h*func(u[i],t,λ)
        t += h
    end
    return u
end

# Backwards Euler is too dependant on the specific differential equation
# function BackwardEuler(N,T,t0,u0,λ)

# Trapezoid too dependant on specific function

function Midpoint2Step(func, N, T, t0, u0, u1, λ=0)
    u = zeros(N+1)
    u[1] = u0
    u[2] = u1
    h = T/N
    t = t0
    for i = 2:N
        t += h
        u[i+1] = u[i-1] + 2*h*func(u[i],t, λ)
    end
    return u
end

function rk4(func, N, T, u0, μ=0)
#=
    4th Order Runga Kutta Method
    Input:
      func(u,t,μ): function of u and t with optional parameter μ.
      N: number of time steps
      T: size of time interval
      u0: initial value
      μ: optional parameter
=#

    h = T/N
    u = zeros(length(u0),N+1)
    t = collect(0:N)*h
    u[:,1] = u0
    for i = 1 : N
        uc = vec(u[:,i])
        tc = t[i]
        k1 = func(uc,tc,μ)
        k2 = func(uc+0.5*h*k1,tc+0.5*h,μ)
        k3 = func(uc+0.5*h*k2,tc+0.5*h,μ)
        k4 = func(uc+h*k3,tc+h,μ)
        u[:,i+1] = uc + h/6.0*(k1+2.0*k2+2.0*k3+k4)
    end
    return u
  
end

function fehlberg(func, N, T, u0, μ=0)
#=
    4th and 5th Fehlberg Method 
    Input:
      func(u,t,μ): function of u and t with optional parameter μ.
      N: number of time steps
      T: size of time interval
      u0: initial value
      μ: optional parameter
=#

    # A = [1/4 0 0 0 0; #  3/32]
    h = T/N
    # 4th order
    u = zeros(length(u0),N+1)
    u[:,1] = u0
    # 5th order
    # û = zeros(length(u0),N+1)
    # û[:,1] = u0
    Eₙ = zeros(N+1)

    t = collect(0:N)*h
    for i = 1 : N
        uc = vec(u[:,i])
        tc = t[i]
        k1 = func(uc,tc,μ)
        k2 = func(uc+(1/4)*h*k1, tc+(1/4)*h, μ)
        k3 = func(uc+h*((3/32)*k1 + (9/32)*k2), tc+(3/8)*h, μ)
        k4 = func(uc+h*((1932/2197)*k1 + (-7200/2197)*k2 + (7296/2197)*k3), tc+(12/13)*h, μ)
        k5 = func(uc+h*((439/216)*k1 + (-8)*k2 + (3680/513)*k3 + (-845/4104)*k4), tc + h, μ)
        k6 = func(uc+h*((-8/27)*k1 + (2)*k2 + (-3544/2565)*k3 + (1859/4104)*k4 + (-11/40)*k5), tc + (1/2)*h, μ)
        # u[:,i+1] = uc + h/6.0*(k1+2.0*k2+2.0*k3+k4)
        û = uc + h*((16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 + (-9/50)*k5 + (2/55)*k6)
        # û[:,i+1] = uc + h*((25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 + (-1/5)*k5)
        u[:,i+1] = uc + h*((25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 + (-1/5)*k5)
        
        Eₙ[i] = norm(u[:,i+1] - û)/h
    end
    return u, Eₙ
  
end


end
