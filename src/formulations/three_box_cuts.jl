immutable ThreeBoxSITB <: Cut{3}; groups; end
export ThreeBoxSITB

function getcut(prob::Problem, vars::Variables, ::ThreeBoxSITB, i, j, k)
    cˣ, cʸ = vars.cˣ, vars.cʸ
    ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
    zˣ, zʸ = vars.zˣ, vars.zʸ

    lbˣ, lbʸ = prob.wlb, prob.hlb
    Lˣ, Lʸ = prob.W, prob.H
    @LinearConstraints(begin
        cˣ[i] ≥ 0.5ℓˣ[i] + lbˣ[k]*zˣ[k,i] + lbˣ[j]*(zˣ[j,k] + zˣ[k,i] - 1)
        cʸ[i] ≥ 0.5ℓʸ[i] + lbʸ[k]*zʸ[k,i] + lbʸ[j]*(zʸ[j,k] + zʸ[k,i] - 1)
        cˣ[i] + 0.5ℓˣ[i] + lbˣ[k]*zˣ[i,k] + lbˣ[j]*(zˣ[i,k] + zˣ[k,j] - 1) ≤ Lˣ
        cʸ[i] + 0.5ℓʸ[i] + lbʸ[k]*zʸ[i,k] + lbʸ[j]*(zʸ[i,k] + zʸ[k,j] - 1) ≤ Lʸ
    end)
end

immutable ThreeBoxNonoverlap <: Cut{3}; groups; end
export ThreeBoxNonoverlap

function getcut(prob::Problem, vars::Variables, ::ThreeBoxNonoverlap, i, j, k)
    cˣ, cʸ = vars.cˣ, vars.cʸ
    ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
    zˣ, zʸ = vars.zˣ, vars.zʸ

    lbˣ, lbʸ = prob.wlb, prob.hlb
    Lˣ, Lʸ = prob.W, prob.H
    @LinearConstraints(begin
        cˣ[i] - cˣ[j] + 0.5*(ℓˣ[i]+ℓˣ[j]) ≤ Lˣ*(1-zˣ[i,j]) - lbˣ[k]*(zˣ[i,k]+zˣ[k,j]-1)
        cʸ[i] - cʸ[j] + 0.5*(ℓʸ[i]+ℓʸ[j]) ≤ Lʸ*(1-zʸ[i,j]) - lbʸ[k]*(zʸ[i,k]+zʸ[k,j]-1)
    end)
end
