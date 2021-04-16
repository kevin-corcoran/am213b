using Plots
using LaTeXStrings

# Load DiffyQ Module and required functions
using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ: fehlberg

# differential equation
function func(u, t, μ)
    return [u[2]; μ*(2-exp(u[2]^2)*u[2] - u[1])]
end

# defining variables
u0 = [3.0; 0.5] # y(0) and y'(0)
T = 30; μ = 4.0; h = 0.025;
N = Int(T/h);
t = collect(0:N)*h

# Part 1
(u, Eₙ) = fehlberg(func, N, T, u0, μ)

# E = Eₙ[Eₙ .> 10.0^(-17)]
E = Eₙ[Eₙ .!= 0.0]
plot(t[1:length(E)], E, label = L"E^{Fehlberg}_{n}", yaxis = :log)
xlabel!(L"t")

# Part 2
(u2, Eₙ) = fehlberg(func, 2*N, T, u0, μ)

E = abs.(u[1,1:N] - u2[1,1:2:2*N])/(1-0.5^5)
E = E[E .!= 0.0]
plot(t[1:length(E)], E, label = L"E_n", yaxis =:log)
xlabel!(L"t")