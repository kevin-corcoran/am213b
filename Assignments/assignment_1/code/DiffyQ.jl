module DiffyQ

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

end
