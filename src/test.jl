using LinearAlgebra
using DelimitedFiles
using MatrixMarket

cd(@__DIR__)

include("./logm_iss_mp.jl")
include("./logm_int.jl")

function generate_reference()
    matnames = ["nos4", "bcsstk04", "lund_b"]
    for matname in matnames
        A_sparse = mmread("./matrix/$(matname).mtx")
        A = collect(A_sparse)
        λ = eigvals(A)
        λmax, λmin = maximum(λ), minimum(λ)

        A_tilde_bf = convert(Array{BigFloat,2}, A/sqrt(λmax*λmin))
        Ref_tilde_bf = logm_mp_full(A_tilde_bf)
        Ref_tilde = convert(Array{Float64,2}, Ref_tilde_bf)

        writedlm("./matrix/A_$(matname).txt", A/sqrt(λmax*λmin))
        writedlm("./matrix/Ref_$(matname).txt", Ref_tilde)
    end
end



function main()
    matnames = ["lund_b", "bcsstk04", "nos4"]
    m_list = 5:5:150
    Data = Array{Any,2}(undef, length(m_list)+1, 10)
    imax, jmax = size(Data)
    for i = 1:imax
        for j = 1:jmax
            Data[i,j] = "init"
        end
    end
    Data[1,1] = "m"

    ϵ = 2.0^(-53)

    for (j, matname) in zip(2:2:6, matnames)
        println("=== $(matname) ===")
        A = readdlm("./matrix/A_$(matname).txt")
        Ref = readdlm("./matrix/Ref_$(matname).txt")
        opnorm_Ref = opnorm(Ref)
        λmax = eigmax(A)
        λmax_mat = zeros(1,1)
        λmax_mat[1,1] = λmax
        Data[1,j] = "$(matname)_mat"
        Data[1,j+1] = "$(matname)_scalar"
        Data[2:end,:1] .= m_list
        l, r = logm_de_getlr(A, ϵ)
        μ = log(λmax)^2 + 2π^2
        d0 = asin(sqrt((μ - sqrt(μ^2 - 4*π^4)) / (2*π^2)))
        speed = 2π*d0/(r-l)
        println("Convergence speed: $(speed)")
        for (i, m) in zip(2:length(m_list)+1, m_list)
            print("#")
            m = m_list[i-1]
            X = logm_de(A, m, l, r)
            Data[i,j] = opnorm(X - Ref) / opnorm_Ref
            x = logm_de(λmax_mat, m, l, r)
            Data[i,j+1] = abs(x[1][1] - log(λmax)) / abs(log(λmax))
        end
        println()
    end

    for (j, matname) in zip(8:10, matnames)
        A = readdlm("./matrix/A_$(matname).txt")
        Ref = readdlm("./matrix/Ref_$(matname).txt")
        opnorm_Ref = opnorm(Ref)
        Data[1,j] = "$(matname)_gl"
        for (i, m) in zip(2:length(m_list)+1, m_list)
            m = m_list[i-1]
            X = logm_gl(A, m)
            Data[i,j] = opnorm(X - Ref) / opnorm_Ref
        end
    end

    writedlm("Result_Test2.csv", Data, ',')
end