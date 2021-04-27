# Problem 4: Plot Region of Absolute Stability
using Plots
theme(:mute)

using Pkg
Pkg.activate("RAS")
include("RAS.jl") # Makes sure the module is run before using it
using .RAS: RAS_stabf, RASrk

# "Yellow" is the stability region

# stability function for Heun's
Φ(z) = 1.0 + z + z^2
xs, Z = RAS_stabf(Φ)

contourf(xs, xs, Z, levels = 1)

# RK4
A = [0 0 0 0
    1/2 0 0 0
    0 1/2 0 0
    0 0 1 0]
b = [1/6, 1/3, 1/3, 1/6]

xs, Z = RASrk(A,b)
# plotly()
contourf(xs,xs,Z)

# 2s-DIRK
α = 1-1/sqrt(2)
A = [α 0
    1-α α]
b = [1-α, α]

xs, Z = RASrk(A,b)
contourf(xs,xs,Z)

# 2s-DIRK
α = 0.5
A = [α 0
    1-α α]
b = [1-α, α]

xs, Z = RASrk(A,b)
contourf(xs,xs,Z)


# # "Meshgrid"
# x1d = -5:5
# y1d = -6:6
# X = x1d' .* ones(13)
# Y = ones(11)' .* y1d