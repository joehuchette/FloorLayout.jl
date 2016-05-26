export Area, SOC, Linearize, Formulation, Unary, Binary,
	   BinaryGray, BinaryBlack, Partition4Bit,
	   Metric, Lowerbounds, Objective, CutCount,
	   base_model

abstract Metric
immutable Lowerbounds <: Metric end
immutable Objective   <: Metric end

immutable CutCount{R}
	num_rounds::Int
	metric::Metric

	CutCount() = new(1, Objective())
	CutCount(n::Int,metric::Metric) = new(n,metric)
end

abstract Area
immutable SOC <: Area end
immutable Linearize <: Area
	ε::Float64
end
Linearize() = Linearize(0.005)

abstract Formulation{A<:Area}
immutable Unary{A<:Area} <: Formulation{A}
	area::A
end
abstract Binary{A<:Area} <: Formulation{A}
immutable BinaryGray{A<:Area} <: Binary{A}
	area::A
end
immutable BinaryBlack{A<:Area} <: Binary{A}
	area::A
end
immutable Partition4Bit{A<:Area} <: Formulation{A}
	area::A
end

call{F<:Formulation}(::Type{F}) = F(SOC())

get_metrics(::Objective,   prob::Problem, i, j) = prob.c[i,j]
get_metrics(::Lowerbounds, prob::Problem, i, j) = prob.wlb[i] + prob.wlb[j] + prob.hlb[i] + prob.hlb[j]

function binary_variables(model::Model, N, ::Unary; redundant=false)
	@variable(model, v[i=1:N,j=1:N,1:2;i!=j], Bin)
	# @variable(model, v[1:N,1:N,1:2], Bin)
	@constraint(model, c1[i=1:N,j=(i+1):N], v[i,j,1] + v[j,i,1] + v[i,j,2] + v[j,i,2] == 1)

	arrx, arry = Array(Variable, N, N), Array(Variable, N, N)
	for i in 1:N, j in (i+1):N
		arrx[i,j] = v[i,j,1]
		arrx[j,i] = v[j,i,1]
		arry[i,j] = v[i,j,2]
		arry[j,i] = v[j,i,2]
	end
	return v, arrx, arry
	# return v, v[:,:,1], v[:,:,2]
end

function binary_variables(model::Model, N, ::Partition4Bit; redundant=false)
	@variable(model, v[1:N,1:N,1:2], Bin)
	@constraint(model, c1[i=1:N,j=(i+1):N], v[i,j,1] + v[i,j,2] + v[j,i,1] + v[j,i,2] ≥ 1)
	@constraint(model, c1[i=1:N,j=(i+1):N,s=1:2], v[i,j,s] + v[j,i,s] ≤ 1)

	return v, v[:,:,1], v[:,:,2]
end

##############
# Implements the binary (logarithmic) gray encoding; that is, for each pair of boxes,
# the binary code assigned to the two polyhedra in direction x (resp. y) differ by
# a single bit. This assignment is unique up to symmetry, and we choose
#		i j <= v = (0,0)
#		j i <= v = (1,0)
# 		 j
#		 i  <= v = (0,1)
# 		 i
#		 j  <= v = (1,1)
# In other words, the second bit encodes which direction the separation occurs, while
# the first encodes which box proceeds the other in that direction.
#
# We also return "helper" collections arrˣ and arrʸ such that arrˢ[i,j] = 1 iff
# i precedes j in direction s, and arrˢ[i,j] ≤ 0 otherwise. This can be encoded as an
# affine map ℒᵍ that takes {0,1}² ⟶ R⁴ in the way described above; in particular, this
# means that any valid inequality for the unary encoding of the form aᵀx + bᵀv ≤ c
# for b ≥ 0 can be transformed to a valid inequality for the binary gray encoding.

ℒ{B<:BinaryGray}(v, ::Type{B}) = begin
	[-v[1] - v[2] + 1
	  v[1] - v[2]
	 -v[1] + v[2]
	  v[1] + v[2] - 1]
end

function binary_variables{B<:BinaryGray}(model::Model, N, ::B; redundant=false)
	if redundant
		@variable(model, z⁺[1:N,1:N], Bin)
		@variable(model, z⁻[1:N,1:N], Bin)
	else
		@variable(model, z⁺[i=1:N,(i+1):N], Bin)
		@variable(model, z⁻[i=1:N,(i+1):N], Bin)
	end
	arrˣ = Array(AffExpr, N, N)
	arrʸ = Array(AffExpr, N, N)
	for i in 1:N
		for j in (i+1):N
			if redundant
				@constraint(model, z⁺[i,j] + z⁺[j,i] == 1)
				@constraint(model, z⁻[i,j] == z⁻[j,i])
			end
			val = ℒ([z⁺[i,j],z⁻[i,j]], B)
			arrˣ[i,j] = val[1]
			arrˣ[j,i] = val[2]
			arrʸ[i,j] = val[3]
			arrʸ[j,i] = val[4]
		end
	end
	return (z⁺,z⁻), arrˣ, arrʸ
end

########
#		i j <= v = (0,0)
#		j i <= v = (1,1)
# 		 j
#		 i  <= v = (1,0)
# 		 i
#		 j  <= v = (0,1)

ℒ{B<:BinaryBlack}(v, ::Type{B}) = begin
	[-v[1] - v[2] + 1
	  v[1] + v[2] - 1
	  v[1] - v[2]
	 -v[1] + v[2]]
end

function binary_variables{B<:BinaryBlack}(model::Model, N, ::B; redundant=true)
	if redundant
		@variable(model, z⁺[1:N,1:N], Bin)
		@variable(model, z⁻[1:N,1:N], Bin)
		for i in 1:N, j in (i+1):N
			@constraint(model, z⁺[i,j] + z⁺[j,i] == 1)
			@constraint(model, z⁻[i,j] + z⁻[j,i] == 1)
		end
	else
		@variable(model, tmp⁺[i=1:N,(i+1):N], Bin)
		@variable(model, tmp⁻[i=1:N,(i+1):N], Bin)
		z⁺ = Array(AffExpr, N, N)
		z⁻ = Array(AffExpr, N, N)
		for i in 1:N, j in (i+1):N
			z⁺[i,j] =     tmp⁺[i,j]
			z⁺[j,i] = 1 - tmp⁺[i,j]
			z⁻[i,j] =     tmp⁻[i,j]
			z⁻[j,i] = 1 - tmp⁻[i,j]
		end
	end
	arrˣ = Array(AffExpr, N, N)
	arrʸ = Array(AffExpr, N, N)
	for i in 1:N, j in (i+1):N
		val = ℒ([z⁺[i,j],z⁻[i,j]], B)
		arrˣ[i,j] = val[1]
		arrˣ[j,i] = val[2]
		arrʸ[i,j] = val[3]
		arrʸ[j,i] = val[4]
	end
	return (z⁺,z⁻), arrˣ, arrʸ
end

function base_model{F<:Formulation}(model::Model, prob::Problem, form::F; redundant=true, symbreak=true)
	N   = prob.N
	lbˣ = prob.wlb
	lbʸ = prob.hlb
	ubˣ = prob.wub
	ubʸ = prob.hub
	p   = prob.c
	α   = prob.area
	Lˣ  = prob.W
	Lʸ  = prob.H

	@variable(model, 0 ≤ cˣ[i=1:N] ≤ Lˣ)
	@variable(model, 0 ≤ cʸ[i=1:N] ≤ Lʸ)
	@variable(model, lbˣ[i] ≤ ℓˣ[i=1:N] ≤ ubˣ[i])
	@variable(model, lbʸ[i] ≤ ℓʸ[i=1:N] ≤ ubʸ[i])

	@variable(model, tˣ[i=1:N,(i+1):N] ≥ 0)
	@variable(model, tʸ[i=1:N,(i+1):N] ≥ 0)
	dˣ = Array(Variable, N, N)
	dʸ = Array(Variable, N, N)
	for i in 1:N
		for j in (i+1):N
			dˣ[i,j] = tˣ[i,j]
			dˣ[j,i] = tˣ[i,j]
			dʸ[i,j] = tʸ[i,j]
			dʸ[j,i] = tʸ[i,j]
		end
	end

	@constraints(model, begin
		c3[i=1:N,j=(i+1):N], dˣ[i,j] ≥ cˣ[i] - cˣ[j]
		c4[i=1:N,j=(i+1):N], dˣ[i,j] ≥ cˣ[j] - cˣ[i]
		c5[i=1:N,j=(i+1):N], dʸ[i,j] ≥ cʸ[i] - cʸ[j]
		c6[i=1:N,j=(i+1):N], dʸ[i,j] ≥ cʸ[j] - cʸ[i]
	end)
	@objective(model, Min, sum{p[i,j]*(dˣ[i,j]+dʸ[i,j]), i=1:N, j=(i+1):N})

	v, zˣ, zʸ = binary_variables(model, N, form; redundant=redundant)

	vars = Variables{F}(cˣ, cʸ, ℓˣ, ℓʸ, dˣ, dʸ, v, zˣ, zʸ)

	add_stay_in_the_box(model, prob, form, vars)
	add_area(model, prob, form, vars)

	add_nonoverlap(model, prob, form, vars)

	if symbreak
		# symmetry breaking: box 1 sits above and to the right of box 2
		maxval,linidx = findmax(p)
		ii, jj = ind2sub(size(p),linidx)
		@assert 1 ≤ ii ≤ N
		@assert 1 ≤ jj ≤ N
		@assert p[ii,jj] == maxval
		@assert ii != jj
		@constraints(model, begin
			cˣ[ii] ≤ cˣ[jj]
			cʸ[ii] ≤ cʸ[jj]
			zˣ[jj,ii] ≤ 0
			zʸ[jj,ii] ≤ 0
			(cˣ[jj]-cˣ[ii]) + (cʸ[jj]-cʸ[ii]) ≥ 0.5min(lbˣ[ii]+lbˣ[jj], lbʸ[ii]+lbʸ[jj])
		end)
	end

	return vars
end

function add_area(model::Model, prob::Problem, ::Formulation{SOC}, vars::Variables)
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	@assert (N = length(ℓˣ)) == length(ℓʸ)
	α = prob.area
	@constraint(model, c[i=1:N], ℓˣ[i]*ℓʸ[i] ≥ α[i])
end

function add_area(model::Model, prob::Problem, form::Formulation{Linearize}, vars::Variables)
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	@assert (N = length(ℓˣ)) == length(ℓʸ)
	lbˣ, lbʸ = prob.wlb, prob.hlb
	ubˣ, ubʸ = prob.wub, prob.hub
	α = prob.area
	ε = form.area.ε
    for i in 1:N
        Cᵢ = iceil(log(ubˣ[i]/lbˣ[i]) / log((1+√(ε))/(1-√(ε))))
        fac = (ubˣ[i]/lbˣ[i])^(1/Cᵢ)
        bnd = lbˣ[i]
        for k = 0:Cᵢ
            bnd *= fac
            @constraint(model, -ℓʸ[i] - α[i]/bnd^2*ℓˣ[i] <= -2*α[i]/bnd)
        end
    end

end

function add_nonoverlap(model::Model, prob::Problem, form::Formulation, vars::Variables)
	N = prob.N
	Lˣ, Lʸ = prob.W, prob.H
	cˣ, cʸ = vars.cˣ, vars.cʸ
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	for i in 1:N, j in 1:N
		i == j && continue
		@constraints(model, begin
			cˣ[i] - cˣ[j] + 0.5*(ℓˣ[i]+ℓˣ[j]) ≤ Lˣ*(1 - zˣ[i,j])
			cʸ[i] - cʸ[j] + 0.5*(ℓʸ[i]+ℓʸ[j]) ≤ Lʸ*(1 - zʸ[i,j])
		end)
	end
end

function add_nonoverlap{B<:Partition4Bit}(model::Model, prob::Problem, form::B, vars::Variables)
	N = prob.N
	Lˣ, Lʸ = prob.W, prob.H
	lbˣ, lbʸ = prob.wlb, prob.hlb
	cˣ, cʸ = vars.cˣ, vars.cʸ
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	zˣ, zʸ = vars.zˣ, vars.zʸ
	for i in 1:N, j in (i+1):N
		@constraints(model, begin
			cˣ[i] - cˣ[j] + 0.5*(ℓˣ[i]+ℓˣ[j]) ≤ Lˣ*(1-zˣ[i,j])
			cˣ[j] - cˣ[i] + 0.5*(ℓˣ[i]+ℓˣ[j]) ≤ Lˣ*(1-zˣ[j,i])
			cʸ[i] - cʸ[j] + 0.5*(ℓʸ[i]+ℓʸ[j]) ≤ Lʸ*(1-zʸ[i,j])
			cʸ[j] - cʸ[i] + 0.5*(ℓʸ[i]+ℓʸ[j]) ≤ Lʸ*(1-zʸ[j,i])

			cˣ[i] + 0.5ℓˣ[i] + Lˣ*zˣ[i,j] ≥ cˣ[j] - 0.5ℓˣ[j] + (lbˣ[i]+lbˣ[j])*(zˣ[i,j]+zˣ[j,i])
			cˣ[j] + 0.5ℓˣ[j] + Lˣ*zˣ[j,i] ≥ cˣ[i] - 0.5ℓˣ[i] + (lbˣ[i]+lbˣ[j])*(zˣ[i,j]+zˣ[j,i])
			cʸ[i] + 0.5ℓʸ[i] + Lʸ*zʸ[i,j] ≥ cʸ[j] - 0.5ℓʸ[j] + (lbʸ[i]+lbʸ[j])*(zʸ[i,j]+zʸ[j,i])
			cʸ[j] + 0.5ℓʸ[j] + Lʸ*zʸ[j,i] ≥ cʸ[i] - 0.5ℓʸ[i] + (lbʸ[i]+lbʸ[j])*(zʸ[i,j]+zʸ[j,i])
		end)
	end
end

function add_stay_in_the_box(model::Model, prob::Problem, f::Formulation, vars::Variables)
	cˣ, cʸ = vars.cˣ, vars.cʸ
	ℓˣ, ℓʸ = vars.ℓˣ, vars.ℓʸ
	@assert (N = prob.N) == length(cˣ) == length(cʸ) == length(ℓˣ) == length(ℓʸ)
	Lˣ, Lʸ = prob.W, prob.H
	for i in 1:N
		@constraints(model, begin
			cˣ[i] ≥      0.5ℓˣ[i]
			cˣ[i] ≤ Lˣ - 0.5ℓˣ[i]
			cʸ[i] ≥      0.5ℓʸ[i]
			cʸ[i] ≤ Lʸ - 0.5ℓʸ[i]
		end)
	end
end
