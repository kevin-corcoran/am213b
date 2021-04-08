using Plots

# Load DiffyQ Module
using Pkg
Pkg.activate("DiffyQ")
include("DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ # Newtons Method is defined here

func(u,t, λ) = -u
function Midpoint2step(N, T, t0, u0, u1, λ=0)
    u = zeros(N+1)
    u[1] = u0
    u[2] = u1
    h = T/N
    t = t0
    for i = 2:N
        t += h
        u[i+1] = u[i-1] + 2*h*func(u[i],t, λ)
    end
    return u
end

T = 2; h = 0.2; t0 = 0;
N = Int(T/h)
uexact(t) = exp(-t)
u0 = 1
u1 = uexact(h)

umid = Midpoint2step(N, T, t0, u0, u1)
ueul = DiffyQ.Euler(func, N, T, t0, u0)

tList = collect(0:N)*(T/N)
plot(tList,umid)
plot!(tList,ueul)