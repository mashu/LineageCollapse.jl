abstract type AbstractTieBreaker end

struct TieBreaker <: AbstractTieBreaker
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

"""
    MostCommonVdjNtTieBreaker

A tie-breaker that selects the representative based on the most common VDJ_nt
sequence weighted by count, matching igdiscover's clonotypes behavior.

For each lineage, it:
1. Sums the `count` for each unique `vdj_nt` sequence
2. Selects the `vdj_nt` with the highest total count
3. Returns the first row with that `vdj_nt`
"""
struct MostCommonVdjNtTieBreaker <: AbstractTieBreaker end

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

"""
    ByMostCommonVdjNt()

Create a tie-breaker that matches igdiscover's clonotypes representative selection.

Selects the representative by finding the VDJ_nt sequence with the highest total
count across all members of the lineage, then returns the first row with that VDJ_nt.

Requires columns: `vdj_nt`, `count`
"""
ByMostCommonVdjNt() = MostCommonVdjNtTieBreaker()
