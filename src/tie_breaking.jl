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

Tie-breaker that picks a representative by **lineage-wide** abundance of each `vdj_nt`
(igdiscover-style): sum `count` for every distinct `vdj_nt` across **all rows in the lineage**,
take the `vdj_nt` with the largest total, then return the **first** row in the lineage that carries
that `vdj_nt`.

This is **not** the same as `ByVdjCount()`: the latter only compares rows that are already
tied at **maximum clone frequency** and uses a **per-clone** score (see `ByVdjCount()` docstring).
Because `MostCommonVdjNt` looks at the whole lineage, the chosen row need not belong to a
maximum-frequency clone if a rarer clone's `vdj_nt` wins the lineage-wide sum.
"""
struct MostCommonVdjNtTieBreaker <: AbstractTieBreaker end

Base.:+(a::TieBreaker, b::TieBreaker) = TieBreaker(vcat(a.criteria, b.criteria))

"""
    ByVdjCount()

Tie-breaker for `Hardest()` when several **clones** share the same top clone frequency:
among those **candidate rows only**, prefer the larger **`vdj_count`** (then lexicographically
smaller `cdr3`).

Here `vdj_count` is computed **per clone** `(lineage_id, d_region, j_call_first, v_call_first, cdr3)`:
within that clone, group rows by `vdj_nt`, sum `count` (or count rows if there is no `count`
column), and take the **maximum** of those per-`vdj_nt` sums. That single number is the clone's
`vdj_count` score for every row in the clone.

Contrast with `ByMostCommonVdjNt()`, which sums `vdj_nt` across the **entire lineage** and
may return a row that is **not** among the maximum-frequency clones.
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

Create a tie-breaker that matches igdiscover's clonotypes representative selection: choose the
`vdj_nt` with the highest **lineage-wide** summed `count`, then return the first row with that
`vdj_nt`. See the `MostCommonVdjNtTieBreaker` docstring and compare with `ByVdjCount()`.

Requires `vdj_nt`. If a `count` column is present, it is used as row weight; otherwise each row
counts as 1.
"""
ByMostCommonVdjNt() = MostCommonVdjNtTieBreaker()
