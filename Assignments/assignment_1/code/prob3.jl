# Problem 3
using Plots
using LaTeXStrings

using Pkg
Pkg.activate("DiffyQ") # allows import
include("code/DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ: Newtons


α = 0.9
β = 50000.0
f(x,s) = x - α + β*sinh(x-cos(s-1))
df(x,s) = 1 + β*cosh(x-cos(s-1))
x0 = 2.0

ss = 0.0:0.1:20.0

xs = []
for s in ss
    push!(xs,Newtons(f,df,x0,s))
end


# Subplots? label axis
plot(ss,xs, label = "Newtons",marker = (:dot,2), add_marker=true)
plot!(ss, cos.(ss.-1), label = L"\cos(s-1)")
