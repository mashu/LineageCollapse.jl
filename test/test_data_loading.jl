using Test
using DataFrames
using LineageCollapse

@testset "Data Loading" begin
    @testset "load_data function" begin
        # Create a temporary TSV file for testing
        test_data = """
        sequence_id\tsequence\tv_sequence_end\tj_sequence_start\tcdr3\tv_call\tj_call\tstop_codon
        seq1\tATCGATCG\t3\t6\tCGAT\tIGHV1-1*01\tIGHJ1*01\tfalse
        seq2\tATCGATCT\t3\t6\tCGAT\tIGHV1-2*01\tIGHJ2*01\tfalse
        """
        test_file = tempname() * ".tsv"
        write(test_file, test_data)

        # Test loading data
        df = load_data(test_file)

        @test df isa DataFrame
        @test nrow(df) == 2
        @test ncol(df) == 8
        @test names(df) == ["sequence_id", "sequence", "v_sequence_end", "j_sequence_start", "cdr3", "v_call", "j_call", "stop_codon"]
        @test df.sequence_id == ["seq1", "seq2"]
        @test df.sequence == ["ATCGATCG", "ATCGATCT"]
        @test df.v_sequence_end == [3, 3]
        @test df.j_sequence_start == [6, 6]
        @test df.cdr3 == ["CGAT", "CGAT"]
        @test df.v_call == ["IGHV1-1*01", "IGHV1-2*01"]
        @test df.j_call == ["IGHJ1*01", "IGHJ2*01"]
        @test df.stop_codon == [false, false]

        # Clean up
        rm(test_file)
    end

    @testset "load_data with custom delimiter" begin
        # Create a temporary CSV file for testing
        test_data = """
        sequence_id,sequence,v_sequence_end,j_sequence_start,cdr3,v_call,j_call,stop_codon
        seq1,ATCGATCG,3,6,CGAT,IGHV1-1*01,IGHJ1*01,false
        """
        test_file = tempname() * ".csv"
        write(test_file, test_data)

        # Test loading data with custom delimiter
        df = load_data(test_file, delimiter=',')

        @test df isa DataFrame
        @test nrow(df) == 1
        @test ncol(df) == 8
        @test names(df) == ["sequence_id", "sequence", "v_sequence_end", "j_sequence_start", "cdr3", "v_call", "j_call", "stop_codon"]

        # Clean up
        rm(test_file)
    end

    @testset "load_data with missing required columns" begin
        # Create a temporary TSV file with missing columns
        test_data = """
        sequence_id\tsequence\tv_sequence_end\tj_sequence_start
        seq1\tATCGATCG\t3\t6
        """
        test_file = tempname() * ".tsv"
        write(test_file, test_data)

        # Test that loading data with missing columns throws an error
        @test_throws Exception load_data(test_file)

        # Clean up
        rm(test_file)
    end
end