using Test
using DataFrames
using LineageCollapse

@testset "Preprocessing" begin
    @testset "preprocess_data function" begin
        # Create a sample DataFrame for testing
        df = DataFrame(
            sequence_id = ["seq1", "seq2", "seq3", "seq4", "seq5"],
            sequence = ["ATCGATCG", "ATCGATCT", "ATCGATGG", "ATCTATCG", "ATCGATCG"],
            v_call = ["IGHV1-1*01", "IGHV1-2*01", "IGHV1-3*01", "IGHV1-4*01", "IGHV1-1*01"],
            j_call = ["IGHJ1*01", "IGHJ2*01", "IGHJ3*01", "IGHJ4*01", "IGHJ1*01"],
            cdr3 = ["CGAT", "CGAT", missing, "CTAT", "CGAT"],
            v_sequence_end = [3, 3, 3, 3, 3],
            j_sequence_start = [6, 6, 6, 6, 6],
            stop_codon = [false, false, false, true, false]
        )

        result = preprocess_data(df)

        @test result isa DataFrame
        @test nrow(result) == 2  # 5 original rows - 1 missing CDR3 - 1 stop codon - 1 duplicate sequence
        @test :d_region in propertynames(result)
        @test :v_call_first in propertynames(result)
        @test :j_call_first in propertynames(result)
        @test :cdr3_length in propertynames(result)

        @test all(result.d_region .== "GAT")
        @test result.v_call_first == ["IGHV1-1*01", "IGHV1-2*01"]
        @test result.j_call_first == ["IGHJ1*01", "IGHJ2*01"]
        @test all(result.cdr3_length .== 4)
    end

    @testset "preprocess_data with custom min_d_region_length" begin
        df = DataFrame(
            sequence_id = ["seq1", "seq2"],
            sequence = ["ATCGATCG", "ATCGAATCG"],
            v_call = ["IGHV1-1*01", "IGHV1-2*01"],
            j_call = ["IGHJ1*01", "IGHJ2*01"],
            cdr3 = ["CGAT", "CGAAT"],
            v_sequence_end = [3, 3],
            j_sequence_start = [6, 7],
            stop_codon = [false, false]
        )

        result = preprocess_data(df, min_d_region_length=3)

        @test nrow(result) == 1
        @test result.sequence_id == ["seq2"]
        @test result.d_region == ["GAAT"]
    end

    @testset "preprocess_data with all invalid data" begin
        df = DataFrame(
            sequence_id = ["seq1", "seq2"],
            sequence = ["ATCGATCG", "ATCGATCG"],
            v_call = ["IGHV1-1*01", "IGHV1-1*01"],
            j_call = ["IGHJ1*01", "IGHJ1*01"],
            cdr3 = [missing, missing],
            v_sequence_end = [3, 3],
            j_sequence_start = [6, 6],
            stop_codon = [true, true]
        )

        result = preprocess_data(df)

        @test nrow(result) == 0
    end
end