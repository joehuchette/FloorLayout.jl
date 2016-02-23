module FloorLayout

using JuMP, MathProgBase, CPLEX

export Problem,
       Solution,
       yalparser,
       load_yal,
       verify_and_print

type Problem
    N::Int
    W::Real
    H::Real
    wlb
    hlb
    wub
    hub
    area
    aspect
    c
end

type Solution
    prob::Problem
    cˣ
    cʸ
    ℓˣ
    ℓʸ
    cost::Float64
    form
    cuts
    model
    time::Float64
    nvar::Int
    ncons::Int
    stat::Symbol
    nodes::Int
    gap::Float64
    lowerbound::Float64
end

type Variables{Formulation}
    cˣ
    cʸ
    ℓˣ
    ℓʸ
    dˣ
    dʸ
    v  # "raw" binary variables, sizes/interpretations depend on Formulation
    zˣ # zˢ[i,j] is 1 iff box i sits to the left of box j in direction s,
    zʸ # and 0 otherwise
end
unpack(a::Variables) = (a.cˣ,a.cʸ,a.ℓˣ,a.ℓʸ,a.dˣ,a.dʸ,a.v,a.zˣ,a.zʸ)

include("yalparser.jl")

function validate_solution(prob::Problem, sol::Solution)
    N, W, H = prob.N, prob.W, prob.H
    x, y, w, h = sol.cˣ, sol.cʸ, sol.ℓˣ, sol.ℓʸ
    tol = 1e-3*max(W,H)
    for i in 1:N
        if !(prob.wlb[i] - tol <= w[i] <= prob.wub[i] + tol) ||
           !(prob.hlb[i] - tol <= h[i] <= prob.hub[i] + tol)
            return false
        end
        if !(w[i] * h[i] >= prob.area[i] - tol)
            return false
        end
        if x[i] + 0.5w[i] >= W + tol ||
           x[i] - 0.5w[i] <= 0 - tol ||
           y[i] + 0.5h[i] >= H + tol ||
           y[i] - 0.5h[i] <= 0 - tol
            return false
        end
    end
    for i in 1:N
        for j in (i+1):N
            if abs(x[i] - x[j]) <= 0.5*(w[i] + w[j]) - tol &&
               abs(y[i] - y[j]) <= 0.5*(h[i] + h[j]) - tol
                return false
            end
        end
    end
    return true
end

function verify_and_print(fp, prob::Problem, sol::Solution)
    valid = validate_solution(prob, sol)
    println(fp, "$(sol.form): ")
    println(fp, "cuts = $(sol.cuts)")
    sol.stat == :Optimal || print(fp, "NOT OPTIMAL ")
    valid || print(fp, "NOT VALID")
    println(fp)
    println(fp, """
        $(sol.time) seconds
        $(sol.nodes) nodes
        $(sol.cost) cost
        $(sol.lowerbound) lowerbound
        $(sol.gap) gap
        $(sol.nvar) variables
        $(sol.ncons) (linear) constraints
    """)
    flush(fp)
end

export @unpack
macro unpack(prob)
    return esc(quote
        n   = $(prob).N
        wlb = $(prob).wlb
        hlb = $(prob).hlb
        wub = $(prob).wub
        hub = $(prob).hub
        c   = $(prob).c
        area = $(prob).area
        W = $(prob).W
        H = $(prob).H
    end)
end

include("get_problem_data.jl")

include("formulations/base.jl")
include("formulations/base_cuts.jl")

end
