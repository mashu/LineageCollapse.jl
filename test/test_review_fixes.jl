using Test
using BioSequences
using DataFrames
using LineageCollapse

function lineage_table(;
        d_region=fill("D", 0),
        lineage_id=Int[],
        j_call_first=String[],
        v_call_first=String[],
        cdr3=String[],
        vdj_nt=String[],
        count=Int[],
        kwargs...)
    df = DataFrame(
        d_region = d_region,
        lineage_id = lineage_id,
        j_call_first = j_call_first,
        v_call_first = v_call_first,
        cdr3 = cdr3,
        vdj_nt = vdj_nt,
        count = count,
    )
    for (k, v) in kwargs
        df[!, k] = v
    end
    return df
end

@testset "Review fixes" begin
    @testset "Hamming rejects unequal lengths" begin
        @test_throws ArgumentError compute_distance(HammingDistance(), dna"ATCG", dna"ATCGATCG")
        @test_throws ArgumentError compute_distance(NormalizedHammingDistance(), dna"ATCG", dna"ATCGATCG")
        @test compute_distance(LevenshteinDistance(), dna"ATCG", dna"ATCGATCG") == 4.0
    end

    @testset "hclust upper-triangle contract" begin
        n = 8
        M = Float32[
            0.0 0.2 0.8 0.3 0.9 0.4 0.7 0.5
            0.2 0.0 0.6 0.1 0.8 0.5 0.9 0.3
            0.8 0.6 0.0 0.4 0.2 0.7 0.3 0.6
            0.3 0.1 0.4 0.0 0.5 0.2 0.8 0.4
            0.9 0.8 0.2 0.5 0.0 0.6 0.1 0.7
            0.4 0.5 0.7 0.2 0.6 0.0 0.3 0.8
            0.7 0.9 0.3 0.8 0.1 0.3 0.0 0.2
            0.5 0.3 0.6 0.4 0.7 0.8 0.2 0.0
        ]
        U = zeros(Float32, n, n)
        for i in 1:n, j in (i+1):n
            U[i, j] = M[i, j]
        end
        method = HierarchicalClustering(0.4f0)
        @test perform_clustering(method, :single, M) == perform_clustering(method, :single, U)
        @test perform_clustering(method, :complete, M) == perform_clustering(method, :complete, U)
        @test perform_clustering(method, :average, M) == perform_clustering(method, :average, U)
    end

    @testset "single-linkage transitivity" begin
        df = DataFrame(
            sequence_id = ["a", "b", "c"],
            v_call_first = fill("V1", 3),
            j_call_first = fill("J1", 3),
            cdr3 = ["AAAA", "AAAT", "AATT"],
            cdr3_length = fill(4, 3),
            d_region = fill("D", 3),
        )
        result = process_lineages(df, 1)
        @test length(unique(result.lineage_id)) == 1
        result_complete = process_lineages(df, 1; linkage=:complete)
        @test length(unique(result_complete.lineage_id)) == 2
    end

    @testset "empty public stages" begin
        empty_raw = DataFrame(
            sequence_id = String[],
            sequence = String[],
            v_call = String[],
            j_call = String[],
            cdr3 = String[],
            v_sequence_end = Int[],
            j_sequence_start = Int[],
            stop_codon = Bool[],
        )
        pre = preprocess_data(empty_raw)
        @test nrow(pre) == 0
        @test :d_region in propertynames(pre)
        @test :v_call_first in propertynames(pre)

        all_missing = DataFrame(
            sequence_id = ["s1"],
            sequence = ["ATCGATCG"],
            v_call = ["IGHV1-1*01"],
            j_call = ["IGHJ1*01"],
            cdr3 = [missing],
            v_sequence_end = [3],
            j_sequence_start = [6],
            stop_codon = [false],
        )
        pre2 = preprocess_data(all_missing)
        @test nrow(pre2) == 0

        processed = process_lineages(pre, 1)
        @test nrow(processed) == 0
        @test :lineage_id in propertynames(processed)
        @test :relative_to_max_cdr3 in propertynames(processed)

        hardest = collapse_lineages(processed, Hardest())
        @test nrow(hardest) == 0
        @test :count in propertynames(hardest)

        ig = collapse_lineages(processed, IgDiscover())
        @test nrow(ig) == 0

        soft = collapse_lineages(processed, Soft(0.2))
        @test nrow(soft) == 0
        @test :clone_frequency in propertynames(soft)
    end

    @testset "invalid AIRR coordinates are dropped" begin
        df = DataFrame(
            sequence_id = ["ok", "bad_v", "bad_j"],
            sequence = ["ATCGATCGAT", "ATCGATCGAT", "ATCGATCGAT"],
            v_call = fill("IGHV1-1*01", 3),
            j_call = fill("IGHJ1*01", 3),
            cdr3 = ["CGAT", "CGAT", "CGAT"],
            v_sequence_end = [3, 99, 3],
            j_sequence_start = [6, 6, 0],
            v_sequence_start = [1, 1, 1],
            j_sequence_end = [10, 10, 10],
            stop_codon = fill(false, 3),
        )
        result = preprocess_data(df)
        @test nrow(result) == 1
        @test result.sequence_id == ["ok"]
        @test result.d_region == ["GAT"]
    end

    @testset "missing stop_codon is kept" begin
        df = DataFrame(
            sequence_id = ["known_ok", "unknown", "stop"],
            sequence = fill("ATCGATCGAT", 3),
            v_call = fill("IGHV1-1*01", 3),
            j_call = fill("IGHJ1*01", 3),
            cdr3 = fill("CGAT", 3),
            v_sequence_end = fill(3, 3),
            j_sequence_start = fill(6, 3),
            v_sequence_start = fill(1, 3),
            j_sequence_end = fill(10, 3),
            stop_codon = [false, missing, true],
        )
        result = preprocess_data(df)
        @test Set(result.sequence_id) == Set(["known_ok", "unknown"])
    end

    @testset "min_d_region_length is inclusive" begin
        df = DataFrame(
            sequence_id = ["seq1", "seq2"],
            sequence = ["ATCGATCG", "ATCGAATCG"],
            v_call = ["IGHV1-1*01", "IGHV1-2*01"],
            j_call = ["IGHJ1*01", "IGHJ2*01"],
            cdr3 = ["CGAT", "CGAAT"],
            v_sequence_end = [3, 3],
            j_sequence_start = [6, 7],
            stop_codon = [false, false],
        )
        both = preprocess_data(df, min_d_region_length=3)
        @test nrow(both) == 2
        none = preprocess_data(df, min_d_region_length=5)
        @test nrow(none) == 0
    end

    @testset "V/J first-allele normalization" begin
        df = DataFrame(
            sequence_id = ["comma", "slash", "plain"],
            sequence = fill("ATCGATCGAT", 3),
            v_call = ["IGHV1-2*01,IGHV1-2*02", "IGHV3-23*01/IGHV3-23D*01", "IGHV1-2*01"],
            j_call = ["IGHJ4*02,IGHJ4*01", "IGHJ4*02", "IGHJ4*02"],
            cdr3 = fill("CGAT", 3),
            v_sequence_end = fill(3, 3),
            j_sequence_start = fill(6, 3),
            v_sequence_start = fill(1, 3),
            j_sequence_end = fill(10, 3),
            stop_codon = fill(false, 3),
        )
        result = preprocess_data(df)
        @test result.v_call_first[1] == "IGHV1-2*01"
        @test result.v_call_first[2] == "IGHV3-23*01/IGHV3-23D*01"
        @test result.v_call_first[3] == "IGHV1-2*01"
        @test result.j_call_first[1] == "IGHJ4*02"
    end

    @testset "Hardest uses max-frequency candidates only" begin
        df = lineage_table(
            d_region = fill("D", 5),
            lineage_id = fill(1, 5),
            j_call_first = fill("J", 5),
            v_call_first = fill("V", 5),
            cdr3 = ["AAA", "AAA", "AAA", "BBB", "BBB"],
            vdj_nt = ["AAA", "AAA", "AAA", "BBB", "BBB"],
            count = [1, 1, 1, 50, 50],
        )
        hardest = collapse_lineages(df, Hardest(); tie_breaker=ByMostCommonVdjNt())
        @test hardest.vdj_nt[1] == "AAA"

        lineage_wide = collapse_lineages(df, Hardest(); tie_breaker=ByLineageMostCommonVdjNt())
        @test lineage_wide.vdj_nt[1] == "BBB"
    end

    @testset "n<=2 shortcut is IgDiscover-only" begin
        df = lineage_table(
            d_region = ["D", "D"],
            lineage_id = [1, 1],
            j_call_first = ["J", "J"],
            v_call_first = ["V", "V"],
            cdr3 = ["AAA", "BBB"],
            vdj_nt = ["LOW", "HIGH"],
            count = [1, 100],
        )
        scientific = collapse_lineages(df, Hardest(); tie_breaker=ByMostCommonVdjNt())
        @test scientific.vdj_nt[1] == "HIGH"

        ig = collapse_lineages(df, IgDiscover())
        @test ig.vdj_nt[1] == "LOW"
        @test ig.count[1] == 101
    end

    @testset "ByLowestVShm is primary and deterministic" begin
        df = lineage_table(
            d_region = fill("D", 3),
            lineage_id = fill(1, 3),
            j_call_first = fill("J", 3),
            v_call_first = fill("V", 3),
            cdr3 = ["AAA", "BBB", "CCC"],
            vdj_nt = ["A100", "B000", "C001"],
            count = [100, 30, 20],
            V_SHM = [0.08, 0.00, 0.00],
        )
        result = collapse_lineages(df, Hardest(); tie_breaker=ByLowestVShm())
        @test result.vdj_nt[1] == "B000"
        @test result.count[1] == 150
    end

    @testset "permutation invariance of scientific representative" begin
        df = lineage_table(
            d_region = fill("D", 4),
            lineage_id = fill(1, 4),
            j_call_first = fill("J", 4),
            v_call_first = fill("V", 4),
            cdr3 = ["AAA", "AAA", "BBB", "CCC"],
            vdj_nt = ["WIN", "AAA", "BBB", "CCC"],
            count = [10, 1, 3, 2],
            sequence_id = ["s1", "s2", "s3", "s4"],
        )
        original = collapse_lineages(df, Hardest(); tie_breaker=ByMostCommonVdjNt())
        for perm in (df[[4, 1, 3, 2], :], df[[3, 4, 2, 1], :], df[[2, 3, 4, 1], :])
            got = collapse_lineages(perm, Hardest(); tie_breaker=ByMostCommonVdjNt())
            @test got.vdj_nt[1] == original.vdj_nt[1]
            @test got.sequence_id[1] == original.sequence_id[1]
        end
    end

    @testset "abundance :rows ignores count" begin
        df = lineage_table(
            d_region = fill("D", 3),
            lineage_id = fill(1, 3),
            j_call_first = fill("J", 3),
            v_call_first = fill("V", 3),
            cdr3 = ["AAA", "BBB", "BBB"],
            vdj_nt = ["A", "B", "B"],
            count = [100, 1, 1],
        )
        by_count = collapse_lineages(df, IgDiscover(); abundance=:count)
        @test by_count.vdj_nt[1] == "A"
        @test by_count.count[1] == 102

        by_rows = collapse_lineages(df, IgDiscover(); abundance=:rows)
        @test by_rows.vdj_nt[1] == "B"
        @test by_rows.count[1] == 3
    end

    @testset "single-linkage connected components match hclust" begin
        seqs = [dna"AAAA", dna"AAAT", dna"AATT", dna"TTTT"]
        clusters, _ = cluster_single_linkage(HammingDistance(), seqs, 1.0f0)
        dist = compute_pairwise_distance(HammingDistance(), seqs)
        h = perform_clustering(HierarchicalClustering(1.0f0), :single, dist)
        function same_parts(a, b)
            n = length(a)
            for i in 1:n, j in (i+1):n
                (a[i] == a[j]) == (b[i] == b[j]) || return false
            end
            true
        end
        @test same_parts(clusters, h)
    end

    @testset "synthetic ground-truth repertoire" begin
        df = DataFrame(
            sequence_id = ["c1a", "c1b", "c2", "c3"],
            v_call_first = fill("V1", 4),
            j_call_first = fill("J1", 4),
            cdr3 = ["AAAAAAAA", "AAAATAAA", "GGGGGGGG", "CCCCCCCC"],
            cdr3_length = fill(8, 4),
            d_region = fill("D", 4),
            vdj_nt = ["vdj1", "vdj1b", "vdj2", "vdj3"],
            count = [5, 1, 2, 1],
        )
        lineages = process_lineages(df, 1)
        @test length(unique(lineages.lineage_id)) == 3
        @test only(unique(lineages[lineages.cdr3 .== "AAAAAAAA", :lineage_id])) ==
              only(unique(lineages[lineages.cdr3 .== "AAAATAAA", :lineage_id]))

        collapsed = collapse_lineages(lineages, IgDiscover())
        clone1 = collapsed[collapsed.cdr3 .== "AAAAAAAA", :]
        @test nrow(clone1) == 1
        @test clone1.count[1] == 6
    end
end
