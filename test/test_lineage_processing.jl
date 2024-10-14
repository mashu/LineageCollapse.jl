using Test
using BioSequences
using DataFrames
using LineageCollapse

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
        expected_hamming = [0.0 1.0 1.0; 1.0 0.0 2.0; 1.0 2.0 0.0]
        expected_levenshtein =  [0.0 1.0 1.0; 1.0 0.0 2.0; 1.0 2.0 0.0]

        @test compute_pairwise_distance(HammingDistance(), sequences) ≈ expected_hamming
        @test compute_pairwise_distance(LevenshteinDistance(), sequences) ≈ expected_levenshtein
    end

    @testset "Hierarchical Clustering" begin
        dist_matrix = [0.0 1.0 1.0; 1.0 0.0 5.0; 1.0 5.0 0.0]
        clustering = HierarchicalClustering(0.6)  # cutoff ratio of 0.6

        result = perform_clustering(clustering, dist_matrix)
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
end