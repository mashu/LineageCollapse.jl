using Test
using DataFrames
using BioSequences
using LineageCollapse
using CSV
using CairoMakie

@testset "LineageCollapse.jl" begin
    @testset "Data Loading" begin
        # Create a temporary TSV file for testing
        test_data = """
        sequence\tv_sequence_end\tj_sequence_start\tcdr3\tv_call\tj_call\tstop_codon
        ACGT\t2\t3\tCDR3\tV1\tJ1\tfalse
        TGCA\t3\t4\tCDR3_2\tV2\tJ2\ttrue
        """
        test_file = tempname() * ".tsv"
        write(test_file, test_data)    
        df = load_data(test_file)

        @test size(df) == (2, 7)
        @test names(df) == ["sequence", "v_sequence_end", "j_sequence_start", "cdr3", "v_call", "j_call", "stop_codon"]
        @test df.sequence == ["ACGT", "TGCA"]
    
        # Test with custom columns
        custom_data = """
        seq\tcdr3\tv_call
        ACGT\tCDR3\tV1
        TGCA\tCDR3_2\tV2
        """
        custom_file = tempname() * ".tsv"
        write(custom_file, custom_data)
    
        custom_df = load_data(custom_file, required_columns=[:seq, :cdr3, :v_call])
        @test size(custom_df) == (2, 3)
        @test names(custom_df) == ["seq", "cdr3", "v_call"]
    end

    @testset "Preprocessing" begin
        test_df = DataFrame(
            sequence = ["ACGT", "TGCA", "GTAC", "CATG"],
            v_sequence_end = [2, 3, 2, 3],
            j_sequence_start = [3, 4, 3, 4],
            cdr3 = ["CDR3", "CDR3_2", missing, "CDR3_4"],
            v_call = ["V1,V2", "V2", "V3", "V4"],
            j_call = ["J1", "J2,J3", "J3", "J4"],
            stop_codon = [false, true, false, false]
        )

        processed_df = preprocess_data(test_df, min_d_region_length=1)
        @test nrow(processed_df) == 2  # After dropping missing CDR3 and stop codons
        @test :d_region in propertynames(processed_df)
        @test :v_call_first in propertynames(processed_df)
        @test :j_call_first in propertynames(processed_df)
        @test :cdr3_length in propertynames(processed_df)
        @test processed_df.v_call_first == ["V1", "V4"]
        @test processed_df.j_call_first == ["J1", "J4"]
    end

    @testset "Lineage Processing" begin
        @testset "Hamming Distance" begin
            seq1 = LongDNA{4}("ATCG")
            seq2 = LongDNA{4}("ATTG")
            @test LineageCollapse.hamming(seq1, seq2) == 1.0
        end

        @testset "Pairwise Hamming" begin
            seqs = [LongDNA{4}("ATCG"), LongDNA{4}("ATTG"), LongDNA{4}("ATCG")]
            dist_matrix = LineageCollapse.pairwise_hamming(seqs)
            @test size(dist_matrix) == (3, 3)
            @test dist_matrix[1,2] == dist_matrix[2,1] == 1.0
            @test dist_matrix[1,3] == dist_matrix[3,1] == 0.0
        end

        @testset "Process Lineages" begin
            test_df = DataFrame(
                v_call_first = ["V1", "V1", "V2", "V2"],
                j_call_first = ["J1", "J1", "J2", "J2"],
                cdr3_length = [4, 4, 5, 5],
                cdr3 = ["ATCG", "ATTG", "ATCGA", "ATCGA"],
                d_region = ["TC", "TT", "TCG", "TCG"],
                cluster_size = [1, 1, 1, 1]
            )
            
            processed_df = process_lineages(test_df)
            @test :cluster in propertynames(processed_df)
            @test :cdr3_frequency in propertynames(processed_df)
            @test nrow(processed_df) <= nrow(test_df)  # Some rows might be filtered out
        end
    end

    @testset "Visualization" begin
        test_df = DataFrame(
            cluster_size = [1, 2, 3, 4, 5],
            cdr3_length = [10, 11, 12, 13, 14],
            cdr3_frequency = [0.1, 0.2, 0.3, 0.4, 0.5],
            v_call_first = ["V1", "V2", "V3", "V4", "V5"]
        )

        fig = plot_diagnostics(test_df)
        @test fig isa CairoMakie.Figure
    end
end