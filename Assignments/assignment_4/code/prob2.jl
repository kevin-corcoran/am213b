using Plots
theme(:mute)

using Pkg
Pkg.activate("StabilityRegion")
include("StabilityRegion.jl")
using .StabilityRegion: RAS, RASlmm, RASrk

# "Light" is the stability region

# 3-step Adams-Moulton
α = [0.0, 0.0, -1.0, 1.0]; β = [1/24, -5/14, 19/24, 9/24]
xs,Z = RASlmm(α, β)

contourf(xs, xs, Z, levels = 1)

# 2-step 4th order LMM
α = [-1.0, 0.0, 1.0]; β = [1/3, 4/3, 1/3];
xs,Z = RASlmm(α, β)

contourf(xs, xs, Z, levels = 1)

# 2-step BDF (BDF2)
α = [1/2, -2.0, 3/2]; β = [0.0, 0.0, 1.0];
xs,Z = RASlmm(α, β)

contourf(xs, xs, Z, levels = 1)

# 3-step BDF (BDF3)
α = [-1/3, 3/2, -3.0, 11/6]; β = [0.0, 0.0, 0.0, 1.0];
xs,Z = RASlmm(α, β)

contourf(xs, xs, Z, levels = 1)








# Forward Euler (Test)
α = [-1.0,1.0]; β = [1.0,0.0];
xs,Z = RASlmm(α, β)

contourf(xs, xs, Z, levels = 1)

# stability function for Heun's
Φ(z) = 1.0 + z + z^2
xs, Z = RAS(Φ)

contourf(xs, xs, Z, levels = 1)