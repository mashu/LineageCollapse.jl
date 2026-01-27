"""
Abstract type for tie-breaking strategies used when selecting lineage representatives.
"""
abstract type AbstractTieBreaker end

"""
    TieBreaker <: AbstractTieBreaker

A configurable tie-breaker that sorts candidates by specified column criteria.

# Fields
- `criteria::Vector{Pair{Symbol,Bool}}`: Sorting criteria as column => descending pairs.
  Columns are checked in order; `true` means sort descending (higher is better).

# Example
```julia
# Sort by count descending, then by cdr3 ascending (lexicographic)
TieBreaker([:count => true, :cdr3 => false])
```
"""
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

"""
    ByVdjCount()

Tie-breaker that selects by highest VDJ nucleotide count, then lexicographic CDR3.
"""
ByVdjCount() = TieBreaker([:vdj_count => true, :cdr3 => false])

"""
    ByCdr3Count()

Tie-breaker that selects by highest CDR3 count, then lexicographic CDR3.
"""
ByCdr3Count() = TieBreaker([:cdr3_count => true, :cdr3 => false])

"""
    BySequenceCount()

Tie-breaker that selects by highest sequence count, then lexicographic CDR3.
"""
BySequenceCount() = TieBreaker([:sequence_count => true, :cdr3 => false])

"""
    ByLexicographic()

Tie-breaker that selects by lexicographically smallest CDR3.
"""
ByLexicographic() = TieBreaker([:cdr3 => false])

"""
    ByFirst()

Tie-breaker that selects the first candidate (no sorting).
"""
ByFirst() = TieBreaker(Pair{Symbol,Bool}[])

"""
    ByMostNaive()

Tie-breaker prioritizing sequences closest to germline (highest V/J identity),
then by VDJ count, CDR3 count, and lexicographic CDR3.
"""
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
