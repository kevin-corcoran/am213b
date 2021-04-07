using Plots
using LaTeXStrings

# Problem 3

function Newtons(f, df, x0, s, tol = 10.0^(-8.0), maxiter=10000)
    #=
        Uses Newtons method to find the roots of f(x) = 0
        Input:
            f(x,s): function of x with optional parameter s
            df(x,s): derivative of f with respect to x
            x0: initial guess
            s: parameter
            tol: desired accuracy
            maxiter: maximum number of iterations
    =#

    xn = x0
    for i = 1:maxiter
        xn_1 = xn - f(xn,s)/df(xn,s)
        # Stop once desired accuracy is achieved 
        if (abs(xn_1 - xn) < tol)
            println("success")
            break
        end
        xn = xn_1
    end
    return xn
end

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
