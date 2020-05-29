using LinearAlgebra
using QuadGK

"""
# logm_iss_mp.jl

This code computes matrix logarithm with arbitrary precision arithmetic.
This code is written based on `logm_mp_full.m` from [mplogm](https://github.com/mfasi/mplogm) by M. Fasi.
This code includes only the diagonal Padé approximation (no Taylor version).

## Usage
```julia
X = logm_mp_full(A)
```

or

```
params = (50, 100)
# The first parameter is the maxium number of computing square roots.
# The second parameter is the maxium degree for the Padé approximation.
X = logm_mp_full(A, params)
```
"""


struct Params
    max_sqrt::Int
    max_deg::Int
end


function α(A::Array{T,2}, k::Int, m::Int) where {T <: Number}
    A_f64 = convert(Array{Float64,2}, A)
    p = Int(floor((1 + sqrt(4*(m+k)+5)) / 2))
    A_to_p = A_f64^p
    x = max(norm(A_to_p)^(1/p), norm(A_to_p*A_f64)^(1/(p+1)))
    return convert(typeof(A[1,1]), x)
end


function scalar_error_diagonal(x::Number, m::Int)
    y =  m / (4m-2) * x
    for j = m-1:-1:1
        y = j*x / ((4j+2)*(1+y))
        y = j*x / ((4j-2)*(1+y))
    end
    y = x / (1+y)
    err = abs(y - log(1+x))
    return err
end


function find_min_deg(A::Array{T,2}, ϵ::Number, m0::Int) where {T <: Number}
    m_left = 1
    m_right = m0
    is_found = false
    m = 1

    A_minus_I = A - I
    while !is_found
        m_middle = Int(ceil((m_right + m_left) / 2))
        x = - α(A_minus_I, m_middle, m_middle)

        if abs(m_left - m_right) < 2
            m = m_right
            is_found = true
        elseif x < -1 || scalar_error_diagonal(x, m_middle) >= ϵ
            m_left = m_middle
        else
            m_right = m_middle
        end
    end
    return m
end


function logm_pade(A::Array{T,2}, A_minus_I::Array{T,2}, m::Int) where {T <: Number}
    x, w = gauss(typeof(A[1,1]), m)
    A_plus_I = A + I
    G = zero(A)
    for i = 1:m
        G += w[i] * inv(A_plus_I + x[i]*A_minus_I)
    end
    return A_minus_I * G
end


function sqrtm_db!(A::Array{T,2}, Z::Array{T,2}, P::Array{T,2}, s::Int) where {T <: Number}
    ϵ = eps(typeof(A[1,1]))
    n = size(A, 1)
    tol = sqrt(n) * ϵ/2
    max_iter = 50
    do_scale = true

    X = copy(A)
    X_old = copy(A)
    M = copy(A)

    k = 0
    norm_M_minus_I = norm(convert(Array{Float64,2}, M) - I)
    while norm(M - I) > tol && k <= max_iter
        k += 1
        if do_scale
            μ = abs(det(M))^(-0.5/n)
            X .= μ * X
            M .= μ^2 * M
        end

        M_inv = inv(M)
        X .= 0.5 * X * (I + M_inv)
        M .= 0.5 * (I + 0.5*(M + M_inv))
        do_scale = norm(X - X_old) / norm(X) >= 1e-2 ? true : false
        X_old .= X
        norm_M_minus_I = norm(convert(Array{Float64,2}, M) - I)
    end

    s += 1
    P .= P * (X+I)
    return X, P, s
end


function logm_mp_full(A::Array{T,2}, params::Params) where {T <: Number}
    A_old = copy(A)
    max_sqrt = params.max_sqrt
    max_deg = params.max_deg
    ϵ = eps(typeof(A[1,1]))
    n = size(A, 1)
    ψ(A::Array{T,2}) where {T <: Number} = norm(A-I, 1)

    s = 0
    Z = A - I
    P = diagm(0=>ones(typeof(A[1,1]), n))

    a = - α(A-I, max_deg, max_deg)
    while norm(A-I,1)>= 1 || scalar_error_diagonal(a, max_deg) >= ϵ*ψ(A) && s < max_sqrt
        A, P, s = sqrtm_db!(A, Z, P, s)
        a = - α(A-I, max_deg, max_deg)
    end
    m = find_min_deg(A, ϵ*ψ(A), max_deg)
    is_required_sqrtm = true
    while is_required_sqrtm && s < max_sqrt
        m_new = max(m - 7, 1)
        a = - α(A-I, m_new, m_new) / 2
        bound_abs_2 = scalar_error_diagonal(a, m_new)

        if bound_abs_2 > ϵ*ψ(A)
            is_required_sqrtm = false
        else
            A, P, s = sqrtm_db!(A, Z, P, s)
            m = find_min_deg(A, ϵ*ψ(A), m)
        end
    end
    println("The number of the sqrtm: $(s)")
    println("The degree of Padé approximation: $(m)")

    A_minus_I = s >= 1 ? Z / P : A - I

    X = 2^s * logm_pade(A, A_minus_I, m)
    return X
end


function logm_mp_full(A::Array{T,2}) where {T <: Number}
    params = Params(100, 100)
    X = logm_mp_full(A, params)
    return X
end
