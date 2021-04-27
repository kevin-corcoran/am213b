# Load DiffyQ Module and required functions
using Plots
using LaTeXStrings

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ: s2_DIRK


f(u,t,μ) = -(0.5*exp(20*cos(1.3*t)) * sinh(u-cos(t)))

α = 1 - 1/sqrt(2);
T = 30.0; h = 2.0^(-5); N = Int(T/h);
u0 = 0.0;

u = s2_DIRK(f, N, T, u0, α)

# Part 1
tList = collect(0:N)*(T/N)
plot(tList, u, label = L"u(t)", thickness_scaling =1.25)
xlabel!(L"t")
ylabel!(L"u(t)")
title!(latexstring("2s-DIRK,h=2^{-5},T=",T))
plot!(tList, cos.(tList), label = L"cos(t)")

# Part 2
g(t) = 0.5*exp(20*cos(1.3*t))
plot(abs.(u - cos.(tList)), g.(tList), xaxis=:log, yaxis=:log, legend = false) 
# plot(g.(tList), abs.(u - cos.(tList)), xaxis=:log, yaxis=:log, legend = false) 
title!("loglog plot")