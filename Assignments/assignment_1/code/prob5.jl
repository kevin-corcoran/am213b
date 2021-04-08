using Plots
# Problem 5

# Load DiffyQ Module
using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ # Newtons Method is defined here

func(u, t, λ) = -λ*sinh(u-cos(t-1))

function Trapezoid(N,T,t0, u0, λ)
    # u = spzeros(N)
    u = zeros(N+1)
    u[1] = u0
    h = T/N
    t = t0
    for i = 1:N
        t += h
        # Use Newtons Method to solve nonlinear eqn for u[n+1]
        f(x,s) = u[i] + h/2 * (func(u[i],t-h,λ) - λ*sinh(x-cos(t-1))) - x
        df(x,s) = -λ*h/2*cosh(x-cos(t-1)) - 1
        u[i+1] = Newtons(f,df,2.0,1.0)
    end
    return u
end

T = 10; t0 = 0; u0 = 0; λ = 10.0^(6.0);
# h = 0.1;
h = 2.0^(-11)
N = Int(T/h)

u = Trapezoid(N,T,t0,u0,λ)
tList = collect(0:N)*(T/N)
plot(tList, u)

for i = 7:2:13
    local h = 2.0^(-i)
    N = Int(T/h)
    local u = Trapezoid(N,T,t0,u0,λ)
    local tList = collect(0:N)*(T/N)
    display(plot(tList, u))
end