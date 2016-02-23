immutable Objective1 <: Cut{2}; groups; end
immutable Objective2 <: Cut{2}; groups; end
immutable Objective3 <: Cut{2}; groups; end
immutable Objective4 <: Cut{2}; groups; end
export Objective1, Objective2, Objective3, Objective4

function getcut(prob::Problem, vars::Variables, ::Objective1, i, j)
	Lˣ, Lʸ = prob.W, prob.H
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	@LinearConstraints(begin
		dˣ[i,j] ≥ 0.5*(ℓˣ[i]+ℓˣ[j]) - Lˣ*(1-zˣ[i,j]-zˣ[j,i])
		dʸ[i,j] ≥ 0.5*(ℓʸ[i]+ℓʸ[j]) - Lʸ*(1-zʸ[i,j]-zʸ[j,i])
	end)
end

function getcut(prob::Problem, vars::Variables, ::Objective2, i, j)
	Lˣ,  Lʸ =  prob.W,   prob.H
	lbˣ, lbʸ = prob.wlb, prob.hlb
	cˣ, cʸ = vars.cˣ, vars.cʸ
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	@LinearConstraints(begin
		dˣ[i,j] ≥ cˣ[i] - cˣ[j] + ℓˣ[i] + lbˣ[j]*(zˣ[i,j]+zˣ[j,i]) - Lˣ*(1-zˣ[i,j])
		dˣ[i,j] ≥ cˣ[i] - cˣ[j] + ℓˣ[j] + lbˣ[i]*(zˣ[i,j]+zˣ[j,i]) - Lˣ*(1-zˣ[i,j])
		dˣ[i,j] ≥ cˣ[j] - cˣ[i] + ℓˣ[i] + lbˣ[j]*(zˣ[i,j]+zˣ[j,i]) - Lˣ*(1-zˣ[j,i])
		dˣ[i,j] ≥ cˣ[j] - cˣ[i] + ℓˣ[j] + lbˣ[i]*(zˣ[i,j]+zˣ[j,i]) - Lˣ*(1-zˣ[j,i])
		dʸ[i,j] ≥ cʸ[i] - cʸ[j] + ℓʸ[i] + lbʸ[j]*(zʸ[i,j]+zʸ[j,i]) - Lʸ*(1-zʸ[i,j])
		dʸ[i,j] ≥ cʸ[i] - cʸ[j] + ℓʸ[j] + lbʸ[i]*(zʸ[i,j]+zʸ[j,i]) - Lʸ*(1-zʸ[i,j])
		dʸ[i,j] ≥ cʸ[j] - cʸ[i] + ℓʸ[i] + lbʸ[j]*(zʸ[i,j]+zʸ[j,i]) - Lʸ*(1-zʸ[j,i])
		dʸ[i,j] ≥ cʸ[j] - cʸ[i] + ℓʸ[j] + lbʸ[i]*(zʸ[i,j]+zʸ[j,i]) - Lʸ*(1-zʸ[j,i])
	end)
end

function getcut(prob::Problem, vars::Variables, ::Objective3, i, j)
	lbˣ, lbʸ = prob.wlb, prob.hlb
	cˣ, cʸ = vars.cˣ, vars.cʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	@LinearConstraints(begin
		dˣ[i,j] ≥ cˣ[i] - cˣ[j] + (lbˣ[i]+lbˣ[j])*zˣ[i,j]
		dˣ[i,j] ≥ cˣ[j] - cˣ[i] + (lbˣ[i]+lbˣ[j])*zˣ[j,i]
		dʸ[i,j] ≥ cʸ[i] - cʸ[j] + (lbʸ[i]+lbʸ[j])*zʸ[i,j]
		dʸ[i,j] ≥ cʸ[j] - cʸ[i] + (lbʸ[i]+lbʸ[j])*zʸ[j,i]
	end)
end

function getcut(prob::Problem, vars::Variables, ::Objective4, i, j)
	Lˣ,  Lʸ =  prob.W,   prob.H
	lbˣ, lbʸ = prob.wlb, prob.hlb
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	@LinearConstraints(begin
		2dˣ[i,j] ≥ ℓˣ[i] - Lˣ*(1-zˣ[i,j]-zˣ[j,i]) + lbˣ[j]*(zˣ[i,j]+zˣ[j,i])
		2dˣ[i,j] ≥ ℓˣ[j] - Lˣ*(1-zˣ[i,j]-zˣ[j,i]) + lbˣ[i]*(zˣ[i,j]+zˣ[j,i])
	end)
end

immutable ThreeBoxObjective1 <: Cut{3}; groups; end
immutable ThreeBoxObjective2 <: Cut{3}; groups; end
immutable ThreeBoxObjective3 <: Cut{3}; groups; end
immutable ThreeBoxObjective4 <: Cut{3}; groups; end
export ThreeBoxObjective1, ThreeBoxObjective2, ThreeBoxObjective3, ThreeBoxObjective4

function getcut(prob::Problem, vars::Variables, ::ThreeBoxObjective1, ii, jj, kk)
	Lˣ, Lʸ = prob.W, prob.H
	lbˣ, lbʸ = prob.wlb, prob.hlb
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	cuts = LinearConstraint[]
	for (_i,_j) in combinations([ii,jj,kk],2) # the two outer ones
		for (i,j) in permutations([_i,_j]) # get both orderings
			k = setdiff([ii,jj,kk],[i,j])[1]
			append!(cuts, @LinearConstraints(begin
				dˣ[i,j] ≥ 0.5*(ℓˣ[i]+ℓˣ[j]) - Lˣ*(1-zˣ[i,j]-zˣ[j,i]) + lbˣ[k]*(zˣ[i,k]+zˣ[k,j]-1)
				dʸ[i,j] ≥ 0.5*(ℓʸ[i]+ℓʸ[j]) - Lʸ*(1-zʸ[i,j]-zʸ[j,i]) + lbʸ[k]*(zʸ[i,k]+zʸ[k,j]-1)
			end))
		end
	end
	cuts
end

function getcut(prob::Problem, vars::Variables, ::ThreeBoxObjective2, ii, jj, kk)
	Lˣ,  Lʸ =  prob.W,   prob.H
	lbˣ, lbʸ = prob.wlb, prob.hlb
	cˣ, cʸ = vars.cˣ, vars.cʸ
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	cuts = LinearConstraint[]
	for (_i,_j) in combinations([ii,jj,kk],2) # the two outer ones
		for (i,j) in permutations([_i,_j]) # get both orderings
			k = setdiff([ii,jj,kk],[i,j])[1]
			append!(cuts, @LinearConstraints(begin
				dˣ[i,j] ≥ cˣ[i] - cˣ[j] + ℓˣ[i] + lbˣ[j]*(zˣ[i,j]+zˣ[j,i]) - Lˣ*(1-zˣ[i,j]) + lbˣ[k]*(zˣ[i,k]+zˣ[k,j]-1)
				dˣ[i,j] ≥ cˣ[i] - cˣ[j] + ℓˣ[j] + lbˣ[i]*(zˣ[i,j]+zˣ[j,i]) - Lˣ*(1-zˣ[i,j]) + lbˣ[k]*(zˣ[i,k]+zˣ[k,j]-1)
				dʸ[i,j] ≥ cʸ[i] - cʸ[j] + ℓʸ[i] + lbʸ[j]*(zʸ[i,j]+zʸ[j,i]) - Lʸ*(1-zʸ[i,j]) + lbʸ[k]*(zʸ[i,k]+zʸ[k,j]-1)
				dʸ[i,j] ≥ cʸ[i] - cʸ[j] + ℓʸ[j] + lbʸ[i]*(zʸ[i,j]+zʸ[j,i]) - Lʸ*(1-zʸ[i,j]) + lbʸ[k]*(zʸ[i,k]+zʸ[k,j]-1)
			end))
		end
	end
	cuts
end

function getcut(prob::Problem, vars::Variables, ::ThreeBoxObjective3, ii, jj, kk)
	lbˣ, lbʸ = prob.wlb, prob.hlb
	cˣ, cʸ = vars.cˣ, vars.cʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	cuts = LinearConstraint[]
	for (_i,_j) in combinations([ii,jj,kk],2) # the two outer ones
		for (i,j) in permutations([_i,_j]) # get both orderings
			k = setdiff([ii,jj,kk],[i,j])[1]
			append!(cuts, @LinearConstraints(begin
				dˣ[i,j] ≥ cˣ[i] - cˣ[j] + (lbˣ[i]+lbˣ[j])*zˣ[i,j] + lbˣ[k]*(zˣ[i,k]+zˣ[k,j]-1)
				dʸ[i,j] ≥ cʸ[i] - cʸ[j] + (lbʸ[i]+lbʸ[j])*zʸ[i,j] + lbʸ[k]*(zʸ[i,k]+zʸ[k,j]-1)
			end))
		end
	end
	cuts
end

function getcut(prob::Problem, vars::Variables, ::ThreeBoxObjective4, ii, jj, kk)
	Lˣ,  Lʸ =  prob.W,   prob.H
	lbˣ, lbʸ = prob.wlb, prob.hlb
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	dˣ, dʸ = vars.dˣ, vars.dʸ
	cuts = LinearConstraint[]
	for (_i,_j) in combinations([ii,jj,kk],2) # the two outer ones
		for (i,j) in permutations([_i,_j]) # get both orderings
			k = setdiff([ii,jj,kk],[i,j])[1]
			append!(cuts, @LinearConstraints(begin
				2dˣ[i,j] ≥ ℓˣ[i] - Lˣ*(1-zˣ[i,j]-zˣ[j,i]) + lbˣ[j]*(zˣ[i,j]+zˣ[j,i]) + 2lbˣ[k]*(zˣ[i,k]+zˣ[k,j]-1)
				2dˣ[i,j] ≥ ℓˣ[j] - Lˣ*(1-zˣ[i,j]-zˣ[j,i]) + lbˣ[i]*(zˣ[i,j]+zˣ[j,i]) + 2lbʸ[k]*(zʸ[i,k]+zʸ[k,j]-1)
			end))
		end
	end
	cuts
end
