using Plots
using LaTeXStrings

# Load DiffyQ Module and required functions
using Pkg
Pkg.activate("DiffyQ")
include("code/DiffyQ.jl") # Makes sure the module is run before using it
using .DiffyQ: CompTrapezoid, CompSimpson

# function to be integrated from a to b
f(t) = sqrt(2 + cos(t)^3)*exp(sin(t))

a = 1; b = 3; #= N = 2^(10); =# T = b-a;

# Plotting error routine
NList = 2 .^(2:10)
errTrapeList = zeros(size(NList))
errSimpList = zeros(size(NList))
for i = 1 : length(NList)
    N = NList[i]

    utrape = CompTrapezoid(N,a,b,f)
    utexact = CompTrapezoid(2*N,a,b,f)

    usimps = CompSimpson(N,a,b,f)
    usexact = CompSimpson(2*N,a,b,f)

    # estimate error
    errTrapeList[i] = abs(utrape-utexact)./(1-(1/2^2))
    errSimpList[i] = abs(usimps-usexact)./(1-(1/2^4))
end

plot(T./NList, errTrapeList,label="Trapezoid",xaxis=:log,yaxis=:log, marker = (:dot,5), add_marker = true)
plot!(T./NList, errSimpList,label="Simpson",xaxis=:log,yaxis=:log, marker = (:square,5), add_marker = true)
xlabel!(L"h")
ylabel!("Approximate Error")
