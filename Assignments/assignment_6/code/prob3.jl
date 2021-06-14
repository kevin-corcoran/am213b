using Plots
theme(:mute)

hs = 10 .^(1:0.2:14)
x = 1.45
q_1(x) = (sin.(x .+ hs) .- sin.(x))./hs
q_2(x) = (sin.(x.+hs) .- sin.(x.-hs))./2hs

E_1(x) = abs.(q_1(x) .- cos(x))
E_2(x) = abs.(q_2(x) .- cos(x))

plot(hs, E_1(x), xaxis = :log, yaxis = :log)
plot!(hs, E_2(x), xaxis = :log, yaxis = :log)