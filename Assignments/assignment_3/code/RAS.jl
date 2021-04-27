module RAS

using LinearAlgebra

function RAS_stabf(Φ,N=200)
    # Return region for plotting given stability function ϕ(z)

    # "Meshgrid"
    x1d = range(-5.0,5.0,length=N)
    X = x1d' .* ones(N); # generally, use ones(size(y1d)) and
    Y = ones(N)' .* x1d; # ones(size(x1d))'
    Z = zeros(size(X))

    for j = 1:length(X)
        z = X[j] + im*Y[j]
        T = Φ(z)
        if (abs(T[1]) <= 1.0)
            Z[j] = 1.0
        else
            Z[j] = 0.0
        end
    end

    return x1d, Z
end

function RASrk(A, b, N = 200)
    # Plot stability region for any Runga Kutta method

    # "Meshgrid"
    x1d = range(-5.0,5.0,length=N)
    X = x1d' .* ones(N); # generally, use ones(size(y1d)) and
    Y = ones(N)' .* x1d; # ones(size(x1d))'
    Z = zeros(size(X))

    e = ones(length(b))
    R(z) = 1.0 + z*b'*(I - z*A)^(-1)*e
    for j = 1:length(X)
        z = X[j] + im*Y[j]
        T = R(z)
        if (abs(T[1]) <= 1.0)
            Z[j] = 1.0
        else
            Z[j] = 0.0
        end
    end

    return x1d, Z
end
    
end