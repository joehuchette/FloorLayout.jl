immutable SequencePair <: Cut{3}; groups; end
export SequencePair
SequencePair() = SequencePair(CutCount{3}())

function getcut{B<:BinaryGray}(prob::Problem, vars::Variables{B}, ::SequencePair, i, j, k)
    N = prob.N
    z⁺, z⁻ = vars.v
    return @LinearConstraints(begin
        z⁺[i,j] + z⁺[j,k] + z⁺[k,i] ≥ 1.0
        z⁺[i,j] + z⁺[j,k] + z⁺[k,i] ≤ 2.0
    end)
end

function addcuts{B<:BinaryBlack}(model::Model, prob::Problem, vars::Variables{B}, cut::SequencePair, ranged=false)
    groups = get_groups(prob, cut, cut.groups)
    N = prob.N
    z⁺, z⁻ = vars.v
    for (i,j,k) in groups
        cuttingplanes = (if ranged
            [LinearConstraint(z⁺[i,j] + z⁺[j,k] + z⁺[k,i], 1.0, 2.0)
             LinearConstraint(z⁻[i,j] + z⁻[j,k] + z⁻[k,i], 1.0, 2.0)]
        else
            @LinearConstraints(begin
                z⁺[i,k] ≥ z⁺[i,j] + z⁺[j,k] - 1
                z⁺[k,i] ≥ z⁺[j,i] + z⁺[k,j] - 1
                z⁻[i,k] ≥ z⁻[i,j] + z⁻[j,k] - 1
                z⁻[k,i] ≥ z⁻[j,i] + z⁻[k,j] - 1
            end)
        end)
        for c in cuttingplanes
            addConstraint(model, c)
        end
    end
end
