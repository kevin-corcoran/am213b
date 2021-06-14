using Plots
using SparseArrays
theme(:mute)

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl")
using .DiffyQ: ForwardEuler_tsys

# solve ibvp over different spatial grids
function ibvp(N)
    # Space
    x0 = 0.0; x = 2.0;
    # Δx = 0.01; L = x-x0; N = Int(L/Δx); 
    L = x-x0; Δx = L/N;


    # Time
    t0 = 0.0; t = 3.0; 
    T = t-t0;
    Δt = 10.0^(-5); M = Int(T/Δt);

    # initial condition
    p(x) = (1 .-0.5*x).^2

    # boundary condition
    α = 0.4
    q(t) = 2*sin(t)^2
    A = spdiagm(-1=>ones(N-1),0=>-2.0*ones(N),1=>ones(N-1))
    A[1,1] = A[1,1] + (2-α*Δx)/(2+α*Δx)
    A = 1/(Δx^2)*A
    b = zeros(N); b[end] = q(N*Δt); b = 1/(Δx^2) * b;

    function b_(t)
        b = zeros(N); b[end] = 2*sin(t)^2;
        return 1/(Δx^2) * b;
    end

    xs = collect(0:N-1)*Δx
    u0 = p(xs)

    # FTCS
    u = ForwardEuler_tsys(M,T,u0,A,b_)
    return u, xs
end

u1, xs1 = ibvp(100)
u2, xs2 = ibvp(200)
u4, xs4 = ibvp(400)

plot(xs1, u1[:,end])
plot!(xs2, u2[:, end])
plot!(xs4, u4[:, end])

# anim = @animate for i ∈ 1:100:size(u)[2]
#     plot(xs[1:N], u[1:N, i], ylims = [0.0, 2.0])
# end
# gif(anim, "radheatloss.gif", fps = 100)