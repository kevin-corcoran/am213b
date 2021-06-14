using Plots
using SparseArrays
theme(:mute)

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl")
using .DiffyQ: ForwardEuler_tsys

# Space x
x0 = 0.0; x = 8.0;
Δx = 0.08; L = x-x0; Nx = Int(L/Δx); 

Ax = 1/Δx^2 * spdiagm(-1=>ones(Nx-1),0=>-2.0*ones(Nx),1=>ones(Nx-1))

# Space y
y0 = 0.0; y = 8.0;
Δy = 0.08; L = y-y0; Ny = Int(L/Δy); 

Ay = 1/Δy^2 * spdiagm(-1=>ones(Ny-1),0=>-2.0*ones(Ny),1=>ones(Ny-1))

# Time
t0 = 0.0; t = 2.0; 
T = t-t0;
Δt = 1.25*10.0^(-3); M = Int(T/Δt);

# initial condition
f(x) = zeros(size(x))

# boundary conditions
function g_x(t)
    b = zeros(N);
    b[1] = ;
    b[end] = sin(2*t);
    return 1/(Δx^2) * b;
end

xs = collect(0:Nx-1)*Δx
u0 = 0.5*f(xs)

# 2s DIRK
α = 1 - 1/sqrt(2)
u = s2_DIRK_sys(M, T, u0, α, A, b_)


plot(xs[1:N], u0[1:N])
plot!(xs[1:N], u[1:N, end])

anim = @animate for i ∈ 1:size(u)[2]
    plot(xs[1:N], u[1:N, i])
end
gif(anim, "2s_DIRK.gif", fps = 3)


