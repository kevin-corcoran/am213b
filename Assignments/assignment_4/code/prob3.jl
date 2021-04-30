using Plots
theme(:mute)

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl")
using .DiffyQ: rk4, Newtons_n

# Differential equation
f(u, x, μ) = [u[2]; sin(x) + (1+0.5*u[2]^2)*u[1]]

# BVP conditions
a = 0.0; b = 2.0; α = 1; β = 0.5;

# Variables
h = 0.002; T = b-a; N = Int(T/h);

# Find ν for IVP
# ν = (β - α)/(b-a) # guess slope avg rate of change over interval
ν = -1.0

function G(ν)
    u0 = [α; ν]
    u = rk4(f, N, T, u0)
    return u[1,end] - β
end

ν = Newtons_n(G, ν, 0.0)

# Solve ode using ν found above
u0 = [α; ν]
u = rk4(f, N, T, u0)

t = collect(0:N)*h
plot(t, u[1,:])