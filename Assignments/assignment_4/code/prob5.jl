using Plots
using SparseArrays
using LinearAlgebra
theme(:mute)

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl")
using .DiffyQ: rk4

a = 0.0; b = 2.0
T = b-a; N = 1000; 
# h = T / N;
# xs = collect(0:N-1)*h
h = T / (N-0.5);
xs = zeros(N+1)
for i = 1:N+1
    xs[i] = a + ((i-1)-0.5)*h
end

# Discretization
q = -(1 .+ exp.(-sin.(xs)))
L =  1/h^2 * spdiagm(-1=>ones(N),0=>-2.0*ones(N+1),1=>ones(N))
A = L + spdiagm(0=>q)
A[1,1] = -1/h^2

# Boundary conditions and forcing term 
α = 2.5; β = 0.5;
g = -5 .- sin.(xs).^2;

# Incorporate boundary conditions in FDM
g[1] = g[1] + (1/h)*α
g[end] = g[end] - (1/h^2)*β

u = A\g

plot(xs, u)

# verify solution
f(u,x,μ) = [u[2]; -5 - sin(x)^2 + (1 + exp(-sin(x)))*u[1]]
u0 = [u[1] - α*h/2; α]
u = rk4(f, N, T, u0)
plot!(xs, u[1,1:N+1])