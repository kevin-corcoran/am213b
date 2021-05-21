using Plots
using SparseArrays
theme(:mute)

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl")
using .DiffyQ: Trapezoid_sys, BTCS

# Space
x0 = 0.0; x = 2.0;
Δx = 0.01; L = x-x0; N = Int(L/Δx); 

A = 1/Δx^2 * spdiagm(-1=>ones(N-1),0=>-2.0*ones(N),1=>ones(N-1))

# Time
t0 = 0.0; t = 0.2; 
T = t-t0;
Δt = 0.01; M = Int(T/Δt);

# initial condition
f(x) = convert(Array{Float64}, (x .< 1) .& (x .> 0))

# boundary condition
b = zeros(N); b[1] = 1.0; b[end] = 0.0; b = 1/(Δx^2) * b;

xs = collect(0:N-1)*Δx
u0 = f(xs)






# Backward Euler
u = BTCS(M, T, u0, A, b)


plot(xs[1:N], u0[1:N])
plot!(xs[1:N], u[1:N, end])

anim = @animate for i ∈ 1:size(u)[2]
    plot(xs[1:N], u[1:N, i])
end
gif(anim, "BTCS.gif", fps = 3)



# Crank-Nicolson
u = Trapezoid_sys(M, T, u0, A, b)

plot(xs[1:N], u0[1:N])
plot!(xs[1:N], u[1:N, end])

anim = @animate for i ∈ 1:size(u)[2]
    plot(xs[1:N], u[1:N, i])
end
gif(anim, "crank.gif", fps = 3)