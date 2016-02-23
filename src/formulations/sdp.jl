import Mosek

export solveTakouda
function solveTakouda(prob::Problem)

    N   = prob.N
	lbˣ = prob.wlb
	lbʸ = prob.hlb
	ubˣ = prob.wub
	ubʸ = prob.hub
	p   = prob.c
	α   = prob.area
	Lˣ  = prob.W
	Lʸ  = prob.H

    W = W/2
    H = H/2

    β = prob.aspect

    Qˣ = Array(Float64, n, n)
    Qʸ = Array(Float64, n, n)
    Uˣ = Array(Float64, n, n)
    Uʸ = Array(Float64, n, n)
    for i in 1:n, j in 1:n
        Qˣ[i,j] = min(W, wub[i]+wub[j])
        Qʸ[i,j] = min(H, hub[i]+hub[j])
        Uˣ[i,j] = W - 0.5*(wlb[i] + wlb[j])
        Uʸ[i,j] = W - 0.5*(hlb[i] + hlb[j])
    end

    model = Model(solver=MosekSolver())

    @defVar(model, x[1:n])
    @defVar(model, y[1:n])
    @defVar(model, wlb[i] <= w[i=1:n] <= wub[i])
    @defVar(model, hlb[i] <= h[i=1:n] <= hub[i])
    @defVar(model, dˣ[1:n,1:n] >= 0)
    @defVar(model, dʸ[1:n,1:n] >= 0)
    @defVar(model, Sˣ[1:n,1:n] >= 0)
    @defVar(model, Sʸ[1:n,1:n] >= 0)
    tn  = convert(Int,  n*( n-1)/2)
    ttn = convert(Int, tn*(tn-1)/2)
    nvar = 1+2tn+2ttn
    @defVar(model, bin[1:nvar,1:nvar], SDP)
    for i in 1:nvar
        @addConstraint(model, bin[i,i] == 1)
    end

    # functions for referencing elements of σα; #assume access of upper tri
    function full_into_tn(i::Int,j::Int)
        ( 1 <= i <= n && 1 <= j <= n) || throw(BoundsError())
        i == j && error("Whoops")
        if i > j
            i,j = j,i
        end
        return sum(Int[x for x in (n-1):-1:(n-i+1)]) + (j - i)
    end

    function full_into_ttn(i::Int,j::Int)
        (1 <= i <= tn && 1 <= j <= tn) || throw(BoundsError())
        i == j && error("Whoops")
        if i > j
            i,j = j,i
        end
        return sum(Int[x for x in (tn-1):-1:(tn-i+1)]) + (j - i)
    end
    σ(i::Int,j::Int)  = bin[1,1+(full_into_tn(i,j))]
    α(i::Int,j::Int)  = bin[1,1+tn+(full_into_tn(i,j))]
    σα(i,j,k,l) = bin[1+full_into_tn(i,j),1+tn+full_into_tn(k,l)]
    σσ(i,j,k,l) = bin[1+   full_into_tn(i,j), 1+   full_into_tn(k,l)]
    αα(i,j,k,l) = bin[1+tn+full_into_tn(i,j), 1+tn+full_into_tn(k,l)]
    σσαα(i,j,k,l,p,q,r,s) = bin[1+2tn+   full_into_ttn(full_into_tn(i,j),full_into_tn(k,l)),
                                 1+2tn+ttn+full_into_ttn(full_into_tn(p,q),full_into_tn(r,s))]

    for i in 1:tn
        for j in (i+1):tn
            @addConstraint(model, bin[1,1+2tn+    full_into_ttn(i,j)] == bin[1+i,1+j])
            @addConstraint(model, bin[1,1+2tn+ttn+full_into_ttn(i,j)] == bin[1+tn+i,1+tn+j])
        end
    end
    for i in 1:n
        areaconstr = [     h[i]  √(area[i])
                      √(area[i])      w[i] ]
        @addSDPConstraint(model, areaconstr >= 0)

        aspconstr = [β[i]    w[i]
                     w[i] area[i]]
        @addSDPConstraint(model, aspconstr >= 0)
        aspconstr = [β[i]    h[i]
                     h[i] area[i]]
        @addSDPConstraint(model, aspconstr >= 0)

        sitbconstr = [0.5*(W-w[i])        x[i]
                             x[i]  0.5*(W-w[i])]
        @addSDPConstraint(model, sitbconstr >= 0)
        sitbconstr = [0.5*(H-h[i])        y[i]
                             y[i]  0.5*(H-h[i])]
        @addSDPConstraint(model, sitbconstr >= 0)

        for j in (i+1):n
            @addConstraint(model, dˣ[i,j] >= 0.5*(w[i]+w[j]) - 0.5*(1-σ(i,j))*Qˣ[i,j]) # 12
            @addConstraint(model, dʸ[i,j] >= 0.5*(h[i]+h[j]) - 0.5*(1+σ(i,j))*Qʸ[i,j]) # 19
            @addConstraint(model, dˣ[i,j] - 2Sˣ[i,j] == x[j] - x[i]) # 13
            @addConstraint(model, Sˣ[i,j] <= 0.25*(3 - σ(i,j) - α(i,j) - σα(i,j,i,j))*Uˣ[i,j]) # 15
            @addConstraint(model, Sˣ[i,j] + x[j] - x[i] >= 0) # 16
            @addConstraint(model, Sˣ[i,j] + x[j] - x[i] <= 0.25*(3 - σ(i,j) + α(i,j) + σα(i,j,i,j))*Uˣ[i,j]) # 17
            @addConstraint(model, dʸ[i,j] - 2Sʸ[i,j] == y[j] - y[i]) # 20
            @addConstraint(model, Sʸ[i,j] <= 0.25*(3 + σ(i,j) - α(i,j) + σα(i,j,i,j))*Uʸ[i,j]) # 22
            @addConstraint(model, Sʸ[i,j] + y[j] - y[i] >= 0) # 23
            @addConstraint(model, Sʸ[i,j] + y[j] - y[i] <= 0.25*(3 + σ(i,j) + α(i,j) - σα(i,j,i,j))*Uʸ[i,j]) # 24

            for k in (j+1):n
                @addConstraint(model, αα(i,j,i,k) - αα(i,j,j,k) + αα(i,k,j,k) +
                                    -σσ(i,j,i,k) - σσ(i,j,j,k) - σσ(i,k,j,k) +
                                     σσαα(i,j,i,k,i,j,i,k) - σσαα(i,j,i,k,i,j,j,k) + σσαα(i,j,i,k,i,k,j,k) +
                                     σσαα(i,j,j,k,i,j,i,k) - σσαα(i,j,j,k,i,j,j,k) + σσαα(i,j,j,k,i,k,j,k) +
                                     σσαα(i,k,j,k,i,j,i,k) - σσαα(i,k,j,k,i,j,j,k) + σσαα(i,k,j,k,i,k,j,k) == 1)
            end
            # symmetry breaking
            @addConstraint(model, dˣ[i,j] >= 0.25*(wlb[i]+wlb[j]) * (1 + σ(i,j)))
            @addConstraint(model, dʸ[i,j] >= 0.25*(hlb[i]+hlb[j]) * (1 - σ(i,j)))
            # # S valid inequalities
            @addConstraint(model, x[i] + 0.5w[i] <= x[j] - 0.5w[j] + 0.5*(3 - σ(i,j) - α(i,j) - σα(i,j,i,j))*W)
            @addConstraint(model, x[j] + 0.5w[j] <= x[i] - 0.5w[i] + 0.5*(3 - σ(i,j) + α(i,j) + σα(i,j,i,j))*W)
            @addConstraint(model, y[i] + 0.5h[i] <= y[j] - 0.5h[j] + 0.5*(3 + σ(i,j) - α(i,j) - σα(i,j,i,j))*H)
            @addConstraint(model, y[j] + 0.5h[j] <= y[i] - 0.5h[i] + 0.5*(3 + σ(i,j) + α(i,j) + σα(i,j,i,j))*H)
        end
    end

    @setObjective(model, Min, dot(c,dˣ + dʸ))
    buildInternalModel(model)
    putdouparam(MathProgBase.getrawsolver(getInternalModel(model)), MSK_DPAR_OPTIMIZER_MAX_TIME, 1.0)
    elapse = @elapsed (stat = solve(model))
    println("Objective value = $(getObjectiveValue(model)) in $elapse sec")
    return getObjectiveValue(model)
end
