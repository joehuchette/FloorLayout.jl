immutable MellerB2 <: Cut{2}; groups; end
export MellerB2

function getcut(prob::Problem, vars::Variables, ::MellerB2, i, j)
    dˣ, dʸ = vars.dˣ, vars.dʸ
    zˣ, zʸ = vars.zˣ, vars.zʸ
    N = prob.N
    lbˣ, lbʸ = prob.wlb, prob.hlb
    @LinearConstraints(begin
        dˣ[i,j] ≥ (lbˣ[i]+lbˣ[j])/2*(zˣ[i,j]+zˣ[j,i])
        dʸ[i,j] ≥ (lbʸ[i]+lbʸ[j])/2*(zʸ[i,j]+zʸ[j,i])
    end)
end

immutable MellerV2 <: Cut{2}; groups; end
export MellerV2

function getcut(prob::Problem, vars::Variables, ::MellerV2, i, j)
    dˣ, dʸ = vars.dˣ, vars.dʸ
    zˣ, zʸ = vars.zˣ, vars.zʸ
    ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
    N = prob.N
    ubˣ, ubʸ = prob.wub, prob.hub
    Lˣ, Lʸ = prob.W, prob.H
    @LinearConstraints(begin
        2dˣ[i,j] ≥ ℓˣ[i] + ℓˣ[j] - min(ubˣ[i]+ubˣ[j],2Lˣ)*(1-zˣ[i,j]-zˣ[j,i])
        2dʸ[i,j] ≥ ℓʸ[i] + ℓʸ[j] - min(ubʸ[i]+ubʸ[j],2Lʸ)*(1-zʸ[i,j]-zʸ[j,i])
    end)
end
