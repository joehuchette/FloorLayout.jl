immutable TightenedSITB <: Cut{2}; groups; end
export TightenedSITB

function getcut(prob::Problem, vars::Variables, ::TightenedSITB, i, j)
	N = prob.N
	Lˣ,  Lʸ  = prob.W,   prob.H
	lbˣ, lbʸ = prob.wlb, prob.hlb
	(cˣ,cʸ,ℓˣ,ℓʸ,_,_,_,zˣ,zʸ) = unpack(vars)
	@LinearConstraints(begin
		cˣ[i] ≥      0.5ℓˣ[i] + lbˣ[j]*zˣ[j,i]
		cˣ[i] ≤ Lˣ - 0.5ℓˣ[i] - lbˣ[j]*zˣ[i,j]
		cʸ[i] ≥      0.5ℓʸ[i] + lbʸ[j]*zʸ[j,i]
		cʸ[i] ≤ Lʸ - 0.5ℓʸ[i] - lbʸ[j]*zʸ[i,j]
	end)
end
