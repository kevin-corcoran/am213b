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
h = 0.025; T = 30;
N = Int(T/h);
t = collect(0:N)*h

# run rk4 for various values of μ and plot the results
μs = [0.5, 2, 4]
for μ in μs
    u = rk4(func, N, T, u0, μ)
    p1 = plot(t,u[1,:], label = L"y(t)", thickness_scaling = 1.5)
    p2 = plot!(t,u[2,:], label = L"y'(t)")
    xlabel!("t")
    title!(latexstring("rk4,h=0.025,\\mu=",μ))
    display(p2)

    p3 = plot(u[1,:], u[2,:], thickness_scaling = 1.5, legend = false)
    xlabel!(L"y(t)")
    ylabel!(L"y'(t)")
    title!(latexstring("rk4,h=0.025,\\mu=",μ))
    display(p3)
end

# μ = 0.5; 

# u = rk4(func, N, T, u0, μ)

# plot(t,u[1,:], label = L"y(t)")
# plot!(t,u[2,:], label = L"y'(t)")
# xlabel!("t")
# title!(latexstring("rk4,h=0.0.025,\\mu=",μ))

# plot(u[1,:], u[2,:])