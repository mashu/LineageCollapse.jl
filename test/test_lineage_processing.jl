using Test
using BioSequences
using DataFrames
using LineageCollapse
using LinearAlgebra

@testset "Lineage Processing" begin
    @testset "Distance Metrics" begin
        seq1 = dna"ATCG"
        seq2 = dna"ATTG"
        seq3 = dna"ATCGATCG"

        @test compute_distance(HammingDistance(), seq1, seq2) == 1.0
        @test compute_distance(HammingDistance(), seq1, seq1) == 0.0
        @test compute_distance(LevenshteinDistance(), seq1, seq2) == 1.0
        @test compute_distance(LevenshteinDistance(), seq1, seq3) == 4.0
    end

    @testset "Pairwise Distance Computation" begin
        sequences = [dna"ATCT", dna"ATCG", dna"ATTT"]
        expected_hamming = [0.0 1.0 1.0; 0.0 0.0 2.0; 0.0 0.0 0.0]
        expected_levenshtein =  [0.0 1.0 1.0; 0.0 0.0 2.0; 0.0 0.0 0.0]

        @test compute_pairwise_distance(HammingDistance(), sequences) ≈ expected_hamming
        @test compute_pairwise_distance(LevenshteinDistance(), sequences) ≈ expected_levenshtein
    end

    @testset "Hierarchical Clustering" begin
        dist_matrix = [0.0 1.0 1.0; 1.0 0.0 5.0; 1.0 5.0 0.0]
        clustering = HierarchicalClustering(2)  # cutoff of 6

        result = perform_clustering(clustering, :average, dist_matrix)
        @test length(unique(result)) == 2  # Should form 2 clusters
        @test result[1] != result[3]  # First and last sequences should be in different clusters
        @test result[1] == result[2]  # First two sequences should be in the same cluster
    end

    @testset "Process Lineages" begin
        # Create a sample DataFrame
        df = DataFrame(
            sequence_id = ["seq1", "seq2", "seq3", "seq4"],
            sequence = ["ATCGATCG", "ATCGATCT", "ATCGATGG", "ATCTATCG"],
            v_call_first = ["V1", "V1", "V2", "V2"],
            j_call_first = ["J1", "J1", "J2", "J2"],
            cdr3 = ["ATCGATCG", "ATCGATCT", "ATCGATGG", "ATCTATCG"],
            cdr3_length = [8, 8, 8, 8],
            d_region = ["CGAT", "CGAT", "CGAT", "CTAT"],
            v_sequence_end = [3, 3, 3, 3],
            j_sequence_start = [6, 6, 6, 6]
        )

        result = process_lineages(df, distance_metric=HammingDistance(), clustering_method=HierarchicalClustering(0.5))

        @test nrow(result) == 4
        @test :lineage_id in propertynames(result)
        @test :cluster in propertynames(result)
        @test :cluster_size in propertynames(result)
        @test :cdr3_count in propertynames(result)
        @test :cdr3_frequency in propertynames(result)

        # Check if sequences with the same V and J genes are grouped together
        @test result[result.v_call_first .== "V1", :lineage_id] == result[result.v_call_first .== "V1", :lineage_id]
        @test result[result.v_call_first .== "V2", :lineage_id] == result[result.v_call_first .== "V2", :lineage_id]
        @test result[result.v_call_first .== "V1", :lineage_id] != result[result.v_call_first .== "V2", :lineage_id]
    end

    @testset "Process Lineages with Levenshtein Distance" begin
        df = DataFrame(
            sequence_id = ["seq1", "seq2", "seq3", "seq4"],
            sequence = ["ATCGATCG", "ATCGATCT", "ATCGATGG", "ATCTATCG"],
            v_call_first = ["V1", "V1", "V2", "V2"],
            j_call_first = ["J1", "J1", "J2", "J2"],
            cdr3 = ["ATCGATCG", "ATCGATCT", "ATCGATGG", "ATCTATCG"],
            cdr3_length = [8, 8, 8, 8],
            d_region = ["CGAT", "CGAT", "CGAT", "CTAT"],
            v_sequence_end = [3, 3, 3, 3],
            j_sequence_start = [6, 6, 6, 6]
        )

        result = process_lineages(df, distance_metric=LevenshteinDistance(), clustering_method=HierarchicalClustering(0.5))

        @test nrow(result) == 4
        @test :lineage_id in propertynames(result)
    end

    @testset "Process Lineages with Normalized Hamming Distance and check grouping" begin
        df = DataFrame(
            sequence_id = 1:20,
            lineage_id_mismatches_20=[1,1,1,2,3,3,3,3,3,3,4,4,4,5,6,6,6,6,6,6],
            cdr3 = ["CCCCCCCCCCCCCCCCCCCC", "CCCAACCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCGGCCC", "GGGGGGGGGGGGGGGGGGGG", "AAAAAAAAAAAAAAAAAAAT", "AAAAAAAAACAAAAAAAAAA", "AAAAGAAAAAAAAAAAAAAA", "AAAAAAAAAAGGAAAAAAAA", "AAAAACCAAAAAAAAAAAAA", "AATTAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC", "CCCAACCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCGGCCC", "GGGGGGGGGGGGGGGGGGGG", "AAAAAAAAAAAAAAAAAAAT", "AAAAAAAAACAAAAAAAAAA", "AAAAGAAAAAAAAAAAAAAA", "AAAAAAAAAAGGAAAAAAAA", "AAAAACCAAAAAAAAAAAAA", "AATTAAAAAAAAAAAAAAAA"],
            d_region = fill("CGAT", 20),
            cdr3_length = fill(20, 20),
            v_call_first = vcat(fill("V1",10),  fill("V2",10)),
            j_call_first = fill("J1", 20),
            v_sequence_end = fill(3,20),
            j_sequence_start = fill(3,20)
        )
        function rand_index(clustering1, clustering2)
            n = length(clustering1)
            if n != length(clustering2)
                error("Clusterings must have the same length")
            end

            a, b = 0, 0
            for i in 1:n-1
                for j in i+1:n
                    same1 = clustering1[i] == clustering1[j]
                    same2 = clustering2[i] == clustering2[j]
                    if same1 == same2
                        a += 1
                    else
                        b += 1
                    end
                end
            end

            return a / (a + b)
        end

        result = process_lineages(df, distance_metric=NormalizedHammingDistance(), clustering_method=HierarchicalClustering(0.2))
        @test rand_index(result.lineage_id, df.lineage_id_mismatches_20) == 1.0
    end

    @testset "Process Lineages with mismatch threshold" begin
        df = DataFrame(
            sequence_id = ["seq1", "seq2", "seq3"],
            sequence = ["ATCG", "ATGG", "GGGG"],
            v_call_first = ["V1", "V1", "V1"],
            j_call_first = ["J1", "J1", "J1"],
            cdr3 = ["ATCG", "ATGG", "GGGG"],
            cdr3_length = [4, 4, 4],
            d_region = ["CG", "GG", "GG"],
            v_sequence_end = [1, 1, 1],
            j_sequence_start = [3, 3, 3]
        )

        result_abs = process_lineages(df, 1)
        @test result_abs[result_abs.cdr3 .== "ATCG", :lineage_id] ==
              result_abs[result_abs.cdr3 .== "ATGG", :lineage_id]
        @test result_abs[result_abs.cdr3 .== "ATCG", :lineage_id] !=
              result_abs[result_abs.cdr3 .== "GGGG", :lineage_id]

        result_frac = process_lineages(df, 0.25)
        @test result_frac[result_frac.cdr3 .== "ATCG", :lineage_id] ==
              result_frac[result_frac.cdr3 .== "ATGG", :lineage_id]
        @test result_frac[result_frac.cdr3 .== "ATCG", :lineage_id] !=
              result_frac[result_frac.cdr3 .== "GGGG", :lineage_id]
    end

    @testset "collapse_lineages function" begin
        # Sample data for testing
        test_df = DataFrame(
            d_region = ["CGAT", "CGAT", "CGAT", "CGAT", "CGAT", "CGAT"],
            lineage_id = [1, 1, 1, 2, 2, 2],
            j_call_first = ["J1", "J1", "J1", "J2", "J2", "J2"],
            v_call_first = ["V1", "V1", "V1", "V2", "V2", "V2"],
            cdr3 = ["ATCG", "ATCG", "ATTG", "GCTA", "GCTA", "GCTT"],
            vdj_nt = ["VDJ1", "VDJ1", "VDJ2", "VDJ3", "VDJ3", "VDJ4"],
            count = [5, 3, 2, 4, 4, 2]
        )

        @testset "Hardest collapse strategy" begin
            result = collapse_lineages(test_df, Hardest())

            @test nrow(result) == 2  # Should have one row per lineage
            @test result[result.lineage_id .== 1, :cdr3][1] == "ATCG"  # Most frequent for lineage 1
            @test result[result.lineage_id .== 2, :cdr3][1] == "GCTA"  # Most frequent for lineage 2
        end

        @testset "Soft collapse strategy" begin
            result = collapse_lineages(test_df, Soft(0.3))

            @test nrow(result) == 6  # Should keep sequences above 30% frequency per lineage
            @test "ATCG" in result.cdr3  # Should keep ATCG in lineage 1
            @test "ATTG" in result.cdr3  # Should keep ATTG in lineage 1 (20%, but rounded to 30%)
            @test "GCTA" in result.cdr3  # Should keep GCTA in lineage 2
            @test "GCTT" in result.cdr3  # Should keep GCTT in lineage 2 (20%, but rounded to 30%)

            # Check frequencies
            @test result[result.cdr3 .== "ATCG", :clone_frequency][1] ≈ 0.6666666666666666 atol=1e-6
            @test result[result.cdr3 .== "ATTG", :clone_frequency][1] ≈ 0.3333333333333333 atol=1e-6
            @test result[result.cdr3 .== "GCTA", :clone_frequency][1] ≈ 0.6666666666666666 atol=1e-6
            @test result[result.cdr3 .== "GCTT", :clone_frequency][1] ≈ 0.3333333333333333 atol=1e-6
        end

        @testset "Edge cases" begin
            # Single sequence per lineage
            single_seq_df = DataFrame(
                d_region = ["CGAT", "CGAT"],
                lineage_id = [1, 2],
                j_call_first = ["J1", "J2"],
                v_call_first = ["V1", "V2"],
                cdr3 = ["ATCG", "GCTA"],
                vdj_nt = ["VDJ1", "VDJ2"],
                count = [1, 1]
            )
            result = collapse_lineages(single_seq_df, Hardest())
            @test nrow(result) == 2
            @test Set(result.cdr3) == Set(["ATCG", "GCTA"])

            # All sequences with equal frequency
            equal_freq_df = DataFrame(
                d_region = ["CGAT", "CGAT", "CGAT"],
                lineage_id = [1, 1, 1],
                j_call_first = ["J1", "J1", "J1"],
                v_call_first = ["V1", "V1", "V1"],
                cdr3 = ["ATCG", "ATTG", "ATAG"],
                vdj_nt = ["VDJ1", "VDJ2", "VDJ3"],
                count = [1, 1, 1]
            )
            result = collapse_lineages(equal_freq_df, Hardest())
            @test nrow(result) == 1
            @test result.cdr3[1] in ["ATCG", "ATTG", "ATAG"]
        end

        @testset "Hardest tie-breaking" begin
            tie_df = DataFrame(
                d_region = ["D1", "D1", "D2", "D2"],
                lineage_id = [1, 1, 1, 1],
                j_call_first = ["J1", "J1", "J1", "J1"],
                v_call_first = ["V1", "V1", "V1", "V1"],
                cdr3 = ["AAA", "AAA", "BBB", "BBB"],
                vdj_nt = ["VDJ1", "VDJ1", "VDJ2", "VDJ3"],
                cdr3_count = [2, 2, 2, 2],
                v_identity = [0.95, 0.95, 0.99, 0.99],
                j_identity = [0.96, 0.96, 0.98, 0.98],
            )

            result_default = collapse_lineages(tie_df, Hardest())
            @test nrow(result_default) == 1
            @test result_default.cdr3[1] == "AAA"

            result_lex = collapse_lineages(tie_df, Hardest(); tie_breaker=ByLexicographic())
            @test nrow(result_lex) == 1
            @test result_lex.cdr3[1] == "AAA"

            result_cdr3 = collapse_lineages(tie_df, Hardest(); tie_breaker=ByCdr3Count())
            @test nrow(result_cdr3) == 1
            @test result_cdr3.cdr3[1] == "AAA"

            result_naive = collapse_lineages(tie_df, Hardest(); tie_breaker=ByMostNaive())
            @test nrow(result_naive) == 1
            @test result_naive.cdr3[1] == "BBB"

            result_composite_naive = collapse_lineages(tie_df, Hardest();
                tie_breaker=ByMostNaive() + ByVdjCount(),
            )
            @test nrow(result_composite_naive) == 1
            @test result_composite_naive.cdr3[1] == "BBB"

            result_composite_vdj = collapse_lineages(tie_df, Hardest();
                tie_breaker=ByVdjCount() + ByMostNaive(),
            )
            @test nrow(result_composite_vdj) == 1
            @test result_composite_vdj.cdr3[1] == "AAA"

            atol_df = DataFrame(
                d_region = ["D1", "D1", "D1", "D2", "D2"],
                lineage_id = [1, 1, 1, 1, 1],
                j_call_first = ["J1", "J1", "J1", "J1", "J1"],
                v_call_first = ["V1", "V1", "V1", "V1", "V1"],
                cdr3 = ["AAA", "AAA", "AAA", "BBB", "BBB"],
                vdj_nt = ["VDJ1", "VDJ2", "VDJ2", "VDJ3", "VDJ3"],
            )
            result_count = collapse_lineages(atol_df, Hardest();
                tie_breaker=BySequenceCount(),
                tie_atol=0.25,
            )
            @test nrow(result_count) == 1
            @test result_count.cdr3[1] == "AAA"
        end

        @testset "Invalid inputs" begin
            @test_throws ArgumentError Soft(-0.1)
            @test_throws ArgumentError Soft(1.1)
        end

        @testset "ByMostCommonVdjNt tie-breaker (igdiscover compatible)" begin
            # Test case: multiple vdj_nt sequences with different counts
            # igdiscover selects based on total count for each vdj_nt
            vdj_df = DataFrame(
                d_region = ["D1", "D1", "D1", "D1", "D1"],
                lineage_id = [1, 1, 1, 1, 1],
                j_call_first = ["J1", "J1", "J1", "J1", "J1"],
                v_call_first = ["V1", "V1", "V1", "V1", "V1"],
                cdr3 = ["AAA", "AAA", "AAA", "BBB", "BBB"],
                vdj_nt = ["VDJ_A", "VDJ_A", "VDJ_B", "VDJ_C", "VDJ_C"],
                count = [5, 3, 10, 4, 4],  # VDJ_A total=8, VDJ_B=10, VDJ_C=8
                extra_col = ["a1", "a2", "b1", "c1", "c2"]
            )

            result = collapse_lineages(vdj_df, Hardest(); tie_breaker=ByMostCommonVdjNt())
            @test nrow(result) == 1
            # Should select VDJ_B which has the highest total count (10)
            @test result.vdj_nt[1] == "VDJ_B"
            @test result.extra_col[1] == "b1"  # Preserves the row's columns
        end

        @testset "Hardest strategy automatically aggregates count and nVDJ_nt" begin
            agg_df = DataFrame(
                d_region = ["D1", "D1", "D1", "D2", "D2"],
                lineage_id = [1, 1, 1, 2, 2],
                j_call_first = ["J1", "J1", "J1", "J2", "J2"],
                v_call_first = ["V1", "V1", "V1", "V2", "V2"],
                cdr3 = ["AAA", "AAA", "BBB", "CCC", "CCC"],
                vdj_nt = ["VDJ1", "VDJ2", "VDJ3", "VDJ4", "VDJ5"],
                count = [5, 3, 2, 4, 6]
            )

            result = collapse_lineages(agg_df, Hardest(); tie_breaker=ByMostCommonVdjNt())

            @test nrow(result) == 2
            # Lineage 1: count should be 5+3+2=10, nVDJ_nt should be 3
            lin1 = result[result.lineage_id .== 1, :]
            @test lin1.count[1] == 10
            @test lin1.nVDJ_nt[1] == 3

            # Lineage 2: count should be 4+6=10, nVDJ_nt should be 2
            lin2 = result[result.lineage_id .== 2, :]
            @test lin2.count[1] == 10
            @test lin2.nVDJ_nt[1] == 2
        end

        @testset "ByMostCommonVdjNt matches igdiscover representative selection" begin
            # This test mimics igdiscover's representative() function behavior:
            # When n > 2: pick the row with the most common VDJ_nt (count-weighted)
            igdiscover_df = DataFrame(
                d_region = ["D1", "D1", "D1", "D1"],
                lineage_id = [1, 1, 1, 1],
                j_call_first = ["J1", "J1", "J1", "J1"],
                v_call_first = ["V1", "V1", "V1", "V1"],
                cdr3 = ["AAA", "AAA", "BBB", "BBB"],
                vdj_nt = ["MOST_COMMON", "MOST_COMMON", "LESS_COMMON", "LESS_COMMON"],
                count = [10, 5, 3, 2],  # MOST_COMMON=15, LESS_COMMON=5
                V_errors = [1, 2, 3, 4],  # Should preserve value from selected row
                J_errors = [10, 20, 30, 40]
            )

            result = collapse_lineages(igdiscover_df, Hardest(); tie_breaker=ByMostCommonVdjNt())

            @test nrow(result) == 1
            @test result.vdj_nt[1] == "MOST_COMMON"
            @test result.V_errors[1] == 1  # First row with MOST_COMMON vdj_nt
            @test result.J_errors[1] == 10
            @test result.count[1] == 20  # Sum of all counts
            @test result.nVDJ_nt[1] == 2  # Two unique VDJ_nt sequences
        end
    end
end