using Plots
using SparseArrays
using LinearAlgebra
theme(:mute)

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl")
using .DiffyQ: rk4

a = 0.0; b = 2.0
T = b-a; N = 1000; h = T / N
xs = collect(0:N-1)*h

# Discretization
q = -(1 .+ exp.(-sin.(xs)))
# L =  1/h^2 * spdiagm(-1=>ones(N),0=>-2.0*ones(N+1),1=>ones(N))
L =  1/h^2 * spdiagm(-1=>ones(N-1),0=>-2.0*ones(N),1=>ones(N-1))
A = L + spdiagm(0=>q)

# Boundary conditions and forcing term 
α = 2.5; β = 0.5;
g = -5 .- sin.(xs).^2;
A[1,1] += (1/h^2)*(2-h)/(2+h)

# Incorporate boundary conditions in FDM
g[1] = g[1] - (1/h^2)*2*α*h/(2+h)
g[end] = g[end] - (1/h^2)*β

# Solve system
u = A\g
plot(xs, u)

# verify solution
f(u,x,μ) = [u[2]; -5 - sin(x)^2 + (1 + exp(-sin(x)))*u[1]]
u0 = [(2*α*h + (2-h)*u[1])/(2+h); 2*(u[1] - α)/(2+h)]
u = rk4(f, N, T, u0)
plot!(xs[1:N], u[1,1:N])