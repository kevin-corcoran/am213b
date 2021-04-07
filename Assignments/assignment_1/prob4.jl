# Problem 4
using Plots

# Load DiffyQ Module
using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ # Newtons Method is defined here

func(u, t, λ) = -λ*sinh(u-cos(t-1))
function Euler(N,T,t0, u0, λ)
    # u = spzeros(N)
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

function BackwardEuler(N,T,t0,u0,λ)
    u = zeros(N+1)
    u[1] = u0
    h = T/N
    t = t0
    x0 = 2.0
    for i = 1:N
        t += h
        f(x, s) = u[i] + h*(-λ*sinh(x-cos(t-1)))
        df(x,s) = -λ*h*cosh(x-cos(t-1)) - 1

        u[i+1] = Newtons(f, df, x0, 1)
    end

    return u
end

# Foward Euler
T = 2.0^(-10.0); t0 = 0; u0 = 0; λ = 10.0^(6.0);

h0 = 2.0^(-18);
N = Int(T/(h0*2.0^(-2))) 
h = T/N

u = Euler(N,T,t0,u0,λ)
tList = collect(0:N)*(T/N)
# tList = t0:h:T+t0
plot(tList, u)
# xlabel(L"t")
# ylabel(L"\theta(t)")


# Backward Euler
T = 10.0; t0 = 0.0; u0 = 0; λ = 10.0^(6.0); h = 0.1;
N = Int(T/h)

# u = DiffyQ.BackwardEuler(N,T,t0,u0,λ)
u = BackwardEuler(N,T,t0,u0,λ)
tList = collect(0:N)*(T/N)
# tList = t0:h:T+t0
plot(tList, u)