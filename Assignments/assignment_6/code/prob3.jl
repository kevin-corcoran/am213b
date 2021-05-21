using Plots
using LaTeXStrings
using SparseArrays
theme(:mute)

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl")
using .DiffyQ: ForwardEuler_sys

# Space
x0 = 0.0; x = 2.0;
Δx = 0.01; L = x-x0; N = Int(L/Δx); 

A = 1/Δx^2 * spdiagm(-1=>ones(N-1),0=>-2.0*ones(N),1=>ones(N-1))

# Time
t0 = 0.0; t = 0.2; 
T = t-t0;
# Δt = (Δx)^2/2 * 1/0.99; M = Int(T/Δt); # unstable
Δt = (Δx)^2/2 * 1/1.01; M = Int(T/Δt); # stable

# initial condition
f(x) = convert(Array{Float64}, (x .< 1) .& (x .> 0))

# boundary condition
b = zeros(N); b[1] = 1.0; b[end] = 0.0; b = 1/(Δx^2) * b;

xs = collect(0:N-1)*Δx
u0 = f(xs)

u = ForwardEuler_sys(M, T, u0, A, b)

xs = collect(0:N)*Δx

# plot(xs[1:N], u0[1:N])
# plot(xs[1:N], u[1:N, end])

ts = (0:3) .^2 * 1e-2; append!(ts, 0.2)
ms = Int.(floor.((ts .- t0)/Δt)); ms[1] = ms[1]+1;  # indices for ts
# p = plot(xs[1:N], u0[1:N], label = latexstring("t = ", 0 ))
j = 1
for i ∈ ms
    # M = Int(T/ts[i])
    j = j + 1
    plot!(xs[1:N], u[1:N, i], label = latexstring("t = ",ts[j]))
end

plot!(xs[1:N], u[1:N, end], label = latexstring("t = ", 0.2))

# anim = @animate for i ∈ 1:20:size(u)[2]
#     plot(xs[1:N], u[1:N, i])
# end
# gif(anim, "anim.gif", fps = 10) 