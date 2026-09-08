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

Among the **candidate rows already selected by `Hardest()`** (maximum clone frequency),
pick the `vdj_nt` with the largest summed abundance, then a deterministic secondary key.

This is **not** IgDiscover's representative rule. IgDiscover sums `vdj_nt` over the
**entire lineage** and, for lineages of size ≤ 2, returns the first row. Use
`ByLineageMostCommonVdjNt()` or `IgDiscover()` for that behaviour.
"""
struct MostCommonVdjNtTieBreaker <: AbstractTieBreaker end

"""
    LineageMostCommonVdjNtTieBreaker

IgDiscover-style representative: sum abundance for each `vdj_nt` across **all rows in
the lineage**, take the `vdj_nt` with the largest total, return the first row that carries
it.

If `small_lineage_first_row` is true (IgDiscover's `representative()`), lineages with
one or two rows skip the abundance calculation and return the first row in input order.
"""
struct LineageMostCommonVdjNtTieBreaker <: AbstractTieBreaker
    small_lineage_first_row::Bool
end

"""
    LowestVShmTieBreaker

Select the row with the lowest `V_SHM` (least-mutated V) among **all lineage members**.
Ties: higher abundance, then lexicographically smaller `vdj_nt`, then smaller `cdr3`.
Missing `V_SHM` is treated as worst (largest).
"""
struct LowestVShmTieBreaker <: AbstractTieBreaker end

Base.:+(a::TieBreaker, b::TieBreaker) = TieBreaker(vcat(a.criteria, b.criteria))

"""
    ByVdjCount()

Tie-breaker for `Hardest()` when several **clones** share the same top clone frequency:
among those **candidate rows only**, prefer the larger **`vdj_count`** (then lexicographically
smaller `cdr3`).

Here `vdj_count` is computed **per clone** `(lineage_id, d_region, j_call_first, v_call_first, cdr3)`:
within that clone, group rows by `vdj_nt`, sum abundance, and take the **maximum** of those
per-`vdj_nt` sums. That single number is the clone's `vdj_count` score for every row in the clone.

Contrast with `ByLineageMostCommonVdjNt()`, which sums `vdj_nt` across the **entire lineage**.
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
Deterministic and permutation-invariant.
"""
ByLexicographic() = TieBreaker([:cdr3 => false])

"""
    ByFirst()

Tie-breaker that selects the first candidate in **current row order**.

Results therefore depend on input ordering (parallel concatenations, unsorted files,
nondeterministic queries). Prefer `ByLexicographic()` when reproducibility is required.
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

Among **maximum clone-frequency candidates**, choose the `vdj_nt` with the highest summed
abundance, then the lexicographically smallest `vdj_nt`. Requires `vdj_nt`.

This operates only on `Hardest()` candidates. For lineage-wide IgDiscover-style selection
use `ByLineageMostCommonVdjNt()` or `collapse_lineages(df, IgDiscover())`.
"""
ByMostCommonVdjNt() = MostCommonVdjNtTieBreaker()

"""
    ByLineageMostCommonVdjNt(; small_lineage_first_row=false)

Choose the `vdj_nt` with the highest **lineage-wide** summed abundance, then the first row
with that `vdj_nt`.

`small_lineage_first_row=true` reproduces IgDiscover's `n <= 2` shortcut (first row,
ignoring abundance). That shortcut is order-dependent.
"""
ByLineageMostCommonVdjNt(; small_lineage_first_row::Bool=false) =
    LineageMostCommonVdjNtTieBreaker(small_lineage_first_row)

"""
    ByLowestVShm()

Select the least-mutated V sequence (`minimum V_SHM`) from the **whole lineage**.
Requires a `V_SHM` column. Not the same as IgDiscover's clonotype representative
(IgDiscover uses most-common `VDJ_nt`; `V_SHM` is used there for `_mindiffrate` only).
"""
ByLowestVShm() = LowestVShmTieBreaker()
