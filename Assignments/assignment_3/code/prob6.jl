using Plots
using LaTeXStrings

using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ: s2_DIRK, BackwardEuler_n

f(u,t,μ) = -(0.5*exp(20*cos(1.3*t)) * sinh(u-cos(t)))

α = 1 - 1/sqrt(2); T = 30.0; u0 = 0.0;

############################
####       Part 1       ####
############################
hs = 1 ./(2 .^(3:8))
for i = 1 : length(hs)
    N = Int(T/hs[i])
    tList = collect(0:N)*(T/N)

    ## Backward Euler
    u_euler = BackwardEuler_n(f, N, T, u0)
    u_eexact = BackwardEuler_n(f,2*N,T,u0)
    euler_error = abs.(u_euler[1:N] - u_eexact[1:2:2*N])./(1-0.5^1) # first order method
    p1 = plot(tList[2:N], euler_error[2:N], label = "Backward Euler", xaxis=:log, yaxis=:log, marker = (:square,5))
    xaxis!(L"t")
    yaxis!(L"u(t)")
    title!(latexstring("Error Estimate,h=",hs[i]))

    ## s2-DIRK
    u_s2_DIRK = s2_DIRK(f, N, T, u0, α)
    u_Dexact = s2_DIRK(f, 2*N, T, u0, α)
    DIRK_error = abs.(u_s2_DIRK[1:N] - u_Dexact[1:2:2*N])./(1-0.5^2) # second order method
    p2 = plot!(tList[2:N], DIRK_error[2:N], label = "2s-DIRK error", xaxis=:log, yaxis=:log, marker = (:square,5))
    display(p2)
end

############################
####       Part 2       ####
############################
h = 2.0^(-7)
N = Int(T/h)
tList = collect(0:N)*(T/N)

## Backward Euler
u_euler = BackwardEuler_n(f, N, T, u0)
u_eexact = BackwardEuler_n(f,2*N,T,u0)
euler_error = abs.(u_euler[1:N] - u_eexact[1:2:2*N])./(1-0.5^1) # first order method
p1 = plot(tList[2:N], euler_error[2:N], label = "Backward Euler", xaxis=:log, yaxis=:log, marker = (:square,5))
xaxis!(L"t")
yaxis!(L"u(t)")
title!(latexstring("Error Estimate,h=",h))

## s2-DIRK
u_s2_DIRK = s2_DIRK(f, N, T, u0, α)
u_Dexact = s2_DIRK(f, 2*N, T, u0, α)
DIRK_error = abs.(u_s2_DIRK[1:N] - u_Dexact[1:2:2*N])./(1-0.5^2) # second order method
p2 = plot!(tList[2:N], DIRK_error[2:N], label = "2s-DIRK error", xaxis=:log, yaxis=:log, marker = (:square,5))
display(p2)
