using CSV
using DataFrames
using LineageCollapseCLI

input = tempname() * ".tsv"
output = tempname() * ".tsv"
processed = tempname() * ".tsv"

sample = DataFrame(
    sequence_id = ["s1", "s2", "s3", "s4"],
    sequence = ["ATCGATCGAT", "ATCGATCTAT", "ATCGATGGAT", "GGGGATCGAT"],
    v_call = ["IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01", "IGHV2-1*01"],
    j_call = ["IGHJ1*01", "IGHJ1*01", "IGHJ1*01", "IGHJ2*01"],
    cdr3 = ["CGAT", "CGTT", "CGAA", "GGAT"],
    v_sequence_start = [1, 1, 1, 1],
    v_sequence_end = [3, 3, 3, 3],
    j_sequence_start = [6, 6, 6, 6],
    j_sequence_end = [10, 10, 10, 10],
    stop_codon = [false, false, false, false],
    V_errors = [0, 1, 1, 2],
    J_errors = [0, 0, 0, 1],
    count = [10, 4, 2, 1],
    VDJ_nt = ["VDJ1", "VDJ2", "VDJ2", "VDJ3"],
)

CSV.write(input, sample; delim = '\t')

LineageCollapseCLI.main([
    input,
    "--output", output,
    "--processed-output", processed,
    "--cdr3-differences", "1",
    "--collapse", "hardest",
    "--tie-breaker", "most-common-vdj-nt",
    "--quiet",
])

LineageCollapseCLI.main([
    input,
    "--output", output,
    "--no-default-filters",
    "--filter", "V_errors<=2",
    "--cdr3-similarity-percent", "75",
    "--collapse", "soft",
    "--soft-cutoff", "0.25",
    "--tie-breaker", "lexicographic",
    "--quiet",
])
