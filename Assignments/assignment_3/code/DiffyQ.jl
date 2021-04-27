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

function s2_DIRK(func, N, T, u0, α, μ=0.0)
    h = T/N
    u = zeros(N+1)
    u[1] = u0

    A = [α 0
        1-α α]
    b = [1-α, α]
    c = [α
         1]

    k = zeros(2)
    t = collect(0:N)*h

    for i = 1:N
        # solve for k1: Newtons
        # 0 = func(u[i] + α*h*k1,t[i] + α*h,μ) - k1 := g(k1)
        tol = 1.0e-10; x0 = 2.0;
        g1(x,s) = func(u[i] + A[1,1]*h*x, t[i] + c[1]*h, μ) - x
        # hard code derivative
        dg1(x,s) = -α*h*(0.5+exp(20*cos(1.3*(t[i]+α*h))))*cosh(u[i] + α*h*x - cos(t[i]+α*h)) - 1
        k[1] = Newtons(g1, dg1, x0, 0, tol)

        # solve for k2: Newtons
        # k2 = func(u[i] + h*((1-α)*k1 + α*k2), t[i] + h, μ)
        g2(x,s) = func(u[i] + h*(A[2,1]*k[1] + A[2,2]*x), t[i] + c[2]*h, μ) - x
        dg2(x,s) = -α*h*(0.5+exp(20*cos(1.3*(t[i]+h))))*cosh(u[i] + h*((1-α)*k[1] + α*x) - cos(t[i]+h)) - 1
        k[2] = Newtons(g2, dg2, x0, 0, tol)

        # u[:,i+1] = uc + h*((1-α)*k1 + α*k2)
        u[i+1] = u[i] + h*(b'*k)
    end
    return u
end

function s2_DIRK_n(func, N, T, u0, α, μ=0.0)
    h = T/N
    u = zeros(N+1)
    u[1] = u0

    A = [α 0
        1-α α]
    b = [1-α, α]
    c = [α
         1]

    k = zeros(2)
    t = collect(0:N)*h
    for i = 1 : N
        uc = u[i]
        tc = t[i]

        # solve for k1: Newtons
        # 0 = func(uc + α*h*k1,tc + α*h,μ) - k1 := g(k1)
        tol = 1.0e-12; x0 = 2.0;
        g(x) = h*func(uc + A[1,1]*h*x, tc + c[1]*h, μ) - x
        k[1] = Newtons_n(g, x0, 0, tol)

        # solve for k2: Newtons
        # k2 = func(uc+h*((1-α)*k1 + α*k2), tc+h, μ)
        g(x) = h*func(uc + h*(A[2,1]*k[1] + A[2,2]*x), tc + c[2]*h, μ) - x # !!! Newtons_n couldn't solve this?
        k[2] = Newtons_n(g, x0, 0, tol)

        # u[:,i+1] = uc + h*((1-α)*k1 + α*k2)
        u[i+1] = uc + b'*k
    end

    return u
    
end

function BackwardEuler_n(func, N, T, u0, μ=0.0)
    h = T/N
    u = zeros(N+1)
    u[1] = u0

    t = collect(0:N+1)*h
    for i = 1 : N
        tol = 1.0e-12; x0 = 2.0;
        # g(x) = u[i] + h*(-(0.5+exp(20*cos(1.3*t[i+1]))) * sinh(x-cos(t[i+1]))) - x
        g(x) = u[i] + h*func(x, t[i+1], μ) - x
        u[i+1] = Newtons_n(g, x0, 0, tol)
    end

    return u
    
end

function Newtons_n(f, x0, s, tol = 10.0^(-8.0), maxiter=10000)
    # Use Newtons method with numerical differentiation
    # Solve f(x,s) = 0 for x where s is a parameter

    n = 0; # number of iterations
    h = 1.0e-5 # step size for numerical differentiation
    err = 1.0
    
    while ((err > tol) && isfinite(err))
        n += 1
        fp = (f(x0+h) - f(x0-h))/(2*h) # numerical differentiation
        dx = -f(x0)/fp
        err = abs(dx)
        x0 = x0 + dx
        if (n > maxiter)
            break
        end
    end

    return x0 #, n, err <= tol # (ie flag whether error is within tol)
end


end
