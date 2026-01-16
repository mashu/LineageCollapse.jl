abstract type TieBreaker end

struct ByVdjCount <: TieBreaker end
struct ByCdr3Count <: TieBreaker end
struct BySequenceCount <: TieBreaker end
struct ByLexicographic <: TieBreaker end
struct ByFirst <: TieBreaker end

function select_hardest_candidate(candidates::DataFrame, ::ByVdjCount)::DataFrame
    if :vdj_count ∉ propertynames(candidates)
        throw(ArgumentError("vdj_count column is required for ByVdjCount()."))
    end
    sorted = sort(candidates, [:vdj_count, :cdr3], rev=[true, false])
    return sorted[1:1, :]
end

function select_hardest_candidate(candidates::DataFrame, ::ByCdr3Count)::DataFrame
    if :cdr3_count ∉ propertynames(candidates)
        throw(ArgumentError("cdr3_count column is required for ByCdr3Count()."))
    end
    sorted = sort(candidates, [:cdr3_count, :cdr3], rev=[true, false])
    return sorted[1:1, :]
end

function select_hardest_candidate(candidates::DataFrame, ::BySequenceCount)::DataFrame
    if :sequence_count ∉ propertynames(candidates)
        throw(ArgumentError("sequence_count column is required for BySequenceCount()."))
    end
    sorted = sort(candidates, [:sequence_count, :cdr3], rev=[true, false])
    return sorted[1:1, :]
end

function select_hardest_candidate(candidates::DataFrame, ::ByLexicographic)::DataFrame
    sorted = sort(candidates, :cdr3)
    return sorted[1:1, :]
end

function select_hardest_candidate(candidates::DataFrame, ::ByFirst)::DataFrame
    return candidates[1:1, :]
end
