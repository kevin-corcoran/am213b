using Plots
using SparseArrays
using LinearAlgebra
theme(:mute)


T = 2.0; N = 1000; h = T / N
xs = collect(0:N-1)*h

# Laplacian
L = 1/h^2 * spdiagm(-1=>ones(N-1),0=>-2.0*ones(N),1=>ones(N-1))

# Differential equation
α = 2.5; β = 0.5;
g = -5 .- sin.(x).^2;

# Boundary conditions
g[1] = g[1] - (1/h^2)*α
g[N] = g[N] - (1/h^2)*β

u = (L .- (1 .+ exp.(-sin.(x)))*I)\g

u_exact(x) = x + (1+exp(-50))/(1-exp(-100))*(exp.(-25*x) - exp.(25*(x.-2)))
plot(x, u)
# plot!(x,u_exact(x))