using Plots

# Load DiffyQ Module
using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ: Euler, Midpoint2Step # Newtons Method is defined here

func(u,t, λ) = -u

T = 2; h = 0.2; t0 = 0;
N = Int(T/h)
uexact(t) = exp(-t)
u0 = 1
u1 = uexact(h)

umid = Midpoint2Step(func, N, T, t0, u0, u1)
ueul = Euler(func, N, T, t0, u0)

tList = collect(0:N)*(T/N)
plot(tList,umid)
plot!(tList,ueul)