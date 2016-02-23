immutable UpperBound98 <: Cut{2}; groups; end
export UpperBound98

function getcut(prob::Problem, vars::Variables, ::UpperBound98, i, j)
	ubˣ, ubʸ = prob.wub, prob.hub
	Lˣ, Lʸ = prob.W, prob.H

	(cˣ,cʸ,ℓˣ,ℓʸ,_,_,_,zˣ,zʸ) = unpack(vars)

	@LinearConstraints(begin
		cˣ[i] + ubˣ[j]*(1-zˣ[j,i]) ≥ 0.5ℓˣ[i] + ℓˣ[j]
		cˣ[i] + ubˣ[j]*(1-zˣ[j,i]) ≥ 0.5ℓˣ[i] + ℓˣ[j]
		cˣ[j] + ubˣ[i]*(1-zˣ[i,j]) ≥ 0.5ℓˣ[j] + ℓˣ[i]
		cˣ[j] + ubˣ[i]*(1-zˣ[i,j]) ≥ 0.5ℓˣ[j] + ℓˣ[i]
		cʸ[i] + ubʸ[j]*(1-zʸ[j,i]) ≥ 0.5ℓʸ[i] + ℓʸ[j]
		cʸ[i] + ubʸ[j]*(1-zʸ[j,i]) ≥ 0.5ℓʸ[i] + ℓʸ[j]
		cʸ[j] + ubʸ[i]*(1-zʸ[i,j]) ≥ 0.5ℓʸ[j] + ℓʸ[i]
		cʸ[j] + ubʸ[i]*(1-zʸ[i,j]) ≥ 0.5ℓʸ[j] + ℓʸ[i]
	end)
end

immutable UpperBound108 <: Cut{2}; groups; end
export UpperBound108

function getcut{U<:Union{Unary,Partition4Bit}}(prob::Problem, vars::Variables{U}, ::UpperBound108, i, j)
	ubˣ, ubʸ = prob.wub, prob.hub
	Lˣ, Lʸ = prob.W, prob.H

	(_,_,ℓˣ,ℓʸ,_,_,_,zˣ,zʸ) = unpack(vars)

	cuts = LinearConstraint[]
	if ubˣ[i] + ubˣ[j] > Lˣ
		push!(cuts, @LinearConstraint(zʸ[i,j] + zʸ[j,i] ≥ (ℓˣ[i]+ℓˣ[j]-Lˣ) / (ubˣ[i]+ubˣ[j]-Lˣ)))
	end
	if ubʸ[i] + ubʸ[j] > Lʸ
		push!(cuts, @LinearConstraint(zˣ[i,j] + zˣ[j,i] ≥ (ℓʸ[i]+ℓʸ[j]-Lʸ) / (ubʸ[i]+ubʸ[j]-Lʸ)))
	end
	return cuts
end
