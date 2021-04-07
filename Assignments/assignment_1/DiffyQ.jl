module DiffyQ

export Newtons

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
            println("success")
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
# Backwards Euler is too dependend on the specific differential equation
# function BackwardEuler(N,T,t0,u0,λ)

end
