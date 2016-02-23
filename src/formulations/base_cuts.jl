abstract Cut{R}
export Cut, addcuts, getcuts

call{T<:Cut{1}}(v::T) = T(Cut{1}, CutCount{1}())
call{T<:Cut{2}}(v::T) = T(Cut{2}, CutCount{2}())
call{T<:Cut{3}}(v::T) = T(Cut{3}, CutCount{3}())

# inelegant code to extract the "native" order for a cut;
# for example, transitivity is enforced by constraining every
# triplet of boxes, and so extract_cut_order(Transitivity) = 3
extract_cut_order{C<:Cut}(::C) = _cut_order(super(C))
_cut_order{R}(::Type{Cut{R}}) = R

function get_groups{T}(prob::Problem, cut::Cut, cutcount::CutCount{T})
    N = prob.N
    R = extract_cut_order(cut)
    if R < T
        error("Cannot add 𝒪^$R cuts when only 𝒪^$T are available")
    elseif R == T
        return combinations(1:N, R)
    else
        function _lt(p1,p2)
            v1 = mapreduce(+, combinations(p1,2)) do x
                get_metrics(cutcount.metric, prob, x...)
            end
            v2 = mapreduce(+, combinations(p2,2)) do x
                get_metrics(cutcount.metric, prob, x...)
            end
            v1 < v2
        end
        ordered_combos = sort!(collect(combinations(1:N,R)), lt=_lt, rev=true)
        numcuts = min(cutcount.num_rounds*N, length(ordered_combos))
        return ordered_combos[1:numcuts]
    end
end

get_groups(prob::Problem, cut::Cut, v::Vector{Tuple}) = v

function addcuts(model::Model, prob::Problem, vars::Variables, cut::Cut)
    groups = get_groups(prob, cut, cut.groups)
    for group in groups
        cuttingplanes = getcut(prob, vars, cut, group...)
        for c in cuttingplanes
            addConstraint(model, c)
        end
    end
end

function getcuts(model::Model, prob::Problem, vars::Variables, cut::Cut)
    groups = get_groups(prob, cut, cut.groups)
    cuts = LinearConstraint[]
    sizehint!(cuts, length(groups))
    for group in groups
        append!(cuts, getcut(prob, vars, cut, group...))
    end
    return cuts
end

include("integer_variable_inequalities.jl")
include("literature_cuts.jl")
include("objective_cuts.jl")
include("three_box_cuts.jl")
include("tightened_inequalities.jl")
include("upperbound_inequalities.jl")
