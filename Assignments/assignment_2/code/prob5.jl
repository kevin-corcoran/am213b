using Plots
using LaTeXStrings

# Load DiffyQ Module and required functions
using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ: rk4

# differential equation
function func(u, t, μ)
    return [u[2]; μ*(2-exp(u[2]^2)*u[2] - u[1])]
end

# defining variables
u0 = [3.0; 0.5] # y(0) and y'(0)
T = 30; μ = 4.0;
N = Int(T/h);
t = collect(0:N)*h

# Part 1
hs = 1 ./(2.0 .^(3:10))
E_max = zeros(size(hs))
for i = 1 : length(hs)
    N = Int(T/hs[i])
    u_nh = rk4(func, N, T, u0, μ)
    u_2nh2 = rk4(func, 2*N, T, u0, μ)

    # subtract every two from u_2nh2
    E_max[i] = maximum(abs.(u_nh[1,1:N] - u_2nh2[1,1:2:2*N])./(1-0.5^4))
end

plot(hs, E_max, label = L"E_{\max}", xaxis = :log, yaxis = :log, marker = (:square,5), add_marker = true)
xlabel!(L"h")
ylabel!(L"\max||E_n(h)||")

# Part 2
# E_max[E_max .< 5.0*(10.0^(-8))]
# Find index where E_max(h) < 5 x 10^-8
index = findfirst(x->x < 5.0*(10.0^(-8)), E_max)
hs[index]