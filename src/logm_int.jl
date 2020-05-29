using Arpack
using LinearAlgebra
using QuadGK


"""
    logm_gl(A, m)

Compute log(A) by using the m-point Gauss-Legendre quadrature.
"""
function logm_gl(A::Array{Ta,2}, m::Int) where {Ta <: Number}
    t, w = gauss(m)
    n = size(A, 1)
    G = zeros(n, n)
    for i = 1:m
        G .+= w[i] * inv((1+t[i])*A + (1-t[i])*I)
    end
    return (A-I) * G
end



"""
    l, r = logm_de_getlr(A, ϵ, use_arpack=true)

Compute a finite interval [l,r] s.t. ‖log(A) - ∫ₗʳFde(x)dx‖₂/‖log(A)‖₂ ≲ ϵ.
If `use_arpack` is true, this program compute ‖A⁻¹‖₂, ‖A-I‖₂, ρ(A) approximately.
If `use_arpack` is false, we compute the parameters by using `eigvals` and `svdvals`
"""
function logm_de_getlr(A::TA, ϵ::Tϵ; arpack=false) where {TA<:AbstractMatrix, Tϵ<:Real}
    # Step 2-3 in [Alg. 1, Tatsuoka et al., 2020]
    norm_A_minus_I = 1.0
    norm_A_inv = 1.0
    ρ = 1.0
    if arpack == true
        ρ = abs(eigs(A, nev=1, tol=0.001, which=:LM, ncv=3)[1][1])
        norm_A_minus_I = svds(A-I, nsv=1, tol=0.001, ncv=3)[1].S[1]
        norm_A_minus_I = 1 / sqrt(eigs(A'*A, nev=1, tol=0.001, which=:SM, ncv=3)[1][1])
    else
        λ, σ, σ_shift = eigvals(A), svdvals(A), svdvals(A-I)
        ρ = maximum(abs.(λ))
        norm_A_minus_I = maximum(σ_shift)
        norm_A_inv = 1 / minimum(σ)
    end

    θ = 1.0
    if issymmetric(A)
        θ = max(abs(log(ρ)), abs(log(norm_A_inv)))
    else
        θ = abs(log(ρ))
    end

    # Step 4-7
    ϵmax = 3 / θ * (norm_A_minus_I * norm_A_inv) / (1 + norm_A_inv)
    if ϵ > ϵmax
        ϵ = ϵmax / 2
    end
    
    # Step 8 and 10
    a1 = θ * ϵ / 3 / norm_A_minus_I
    a2 = 1 / 2 / norm_A_minus_I
    if a1 < a2
        atanh_2a_minus_1 = (log(2a1) - log(2) - log1p(-a1)) / 2
        l = asinh(2/π * atanh_2a_minus_1)
    else
        l = asinh(2/π * atanh(2*a2-1))
    end
    
    # Step 9 and 11
    δ = θ * ϵ / 3 / norm_A_minus_I / norm_A_inv
    b1 = 1 - δ
    b2 = 2 * norm_A_inv / (2*norm_A_inv+1)
    if b1 > b2
        atanh_2b_minus_1 = (log(2) + log1p(2δ) - log(2δ)) / 2
        r = asinh(2/π * atanh_2b_minus_1)
    else
        r = asinh(atanh(2*b2-1))
    end
    
    return l, r
end



"""
    logm_de(A, m, l, r)

Compute log(A) by using the m-point DE formula on [l,r].
"""
function logm_de(A::Array{Ta,2}, m::Int, l::Ti, r::Ti) where {Ta<:Number, Ti<:Real}
    n = size(A)[1]
    
    function F(x::Float64)
        t = tanh(π/2*sinh(x))
        return π/2*cosh(x)*sech(π/2*sinh(x))^2 * inv((1+t)*A + (1-t)*I)
    end

    x = LinRange(l, r, m)
    
    w = (r-l) / (m - 1) * ones(m)
    w[1], w[m] = w[1]/2, w[m]/2

    Tl = zeros(n,n)
    Tr = zeros(n,n)
    for i = 1:m÷2
        Tl += w[i] * F(x[i])
    end
    for i = m:-1:m÷2+1
        Tr += w[i] * F(x[i])
    end
    return (A-I) * (Tl+Tr)
end