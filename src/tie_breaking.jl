struct TieBreaker
    criteria::Vector{Pair{Symbol,Bool}}
    function TieBreaker(criteria::Vector{Pair{Symbol,Bool}})
        seen = Set{Symbol}()
        filtered = Pair{Symbol,Bool}[]
        for (col, desc) in criteria
            if col ∉ seen
                push!(filtered, col => desc)
                push!(seen, col)
            end
        end
        return new(filtered)
    end
end

Base.:+(a::TieBreaker, b::TieBreaker) = TieBreaker(vcat(a.criteria, b.criteria))

ByVdjCount() = TieBreaker([:vdj_count => true, :cdr3 => false])
ByCdr3Count() = TieBreaker([:cdr3_count => true, :cdr3 => false])
BySequenceCount() = TieBreaker([:sequence_count => true, :cdr3 => false])
ByLexicographic() = TieBreaker([:cdr3 => false])
ByFirst() = TieBreaker(Pair{Symbol,Bool}[])
ByMostNaive() = TieBreaker([
    :v_identity => true,
    :j_identity => true,
    :vdj_count => true,
    :cdr3_count => true,
    :cdr3 => false,
])

function select_hardest_candidate(candidates::DataFrame, tie_breaker::TieBreaker)::DataFrame
    if isempty(tie_breaker.criteria)
        return candidates[1:1, :]
    end

    missing_cols = [
        col for (col, _) in tie_breaker.criteria
        if col ∉ propertynames(candidates)
    ]
    if !isempty(missing_cols)
        throw(ArgumentError("Missing required columns for tie breaking: $(join(string.(missing_cols), ", "))."))
    end

    cols = [col for (col, _) in tie_breaker.criteria]
    revs = [desc for (_, desc) in tie_breaker.criteria]
    sorted = sort(candidates, cols, rev=revs)
    return sorted[1:1, :]
end
