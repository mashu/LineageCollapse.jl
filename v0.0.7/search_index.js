var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = LineageCollapse","category":"page"},{"location":"#LineageCollapse","page":"Home","title":"LineageCollapse","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LineageCollapse, a Julia package for collapsing lineages in AIRR data.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"LineageCollapse provides tools for processing and analyzing adaptive immune receptor repertoire (AIRR) data. It offers functions for data loading, preprocessing, lineage collapsing, and visualization.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install LineageCollapse using Julia's package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"LineageCollapse\")","category":"page"},{"location":"#Basic-Usage","page":"Home","title":"Basic Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Here's a quick example of how to use LineageCollapse:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using LineageCollapse\n\n# Load data\ndf = load_data(\"path/to/your/airr_data.tsv.gz\")\n\n# Preprocess data\npreprocessed_df = preprocess_data(df, min_d_region_length=3)\n\n# Assign lineages using default length Normalized Hamming distance and Hierarchical clustering\nresult1 = process_lineages(preprocessed_df)\n\n# Assign lineages using length Normalized Hamming distance\n# with Hierarchical clustering but different CDR3 similarity cutoff\nresult2 = process_lineages(preprocessed_df, \n                                    distance_metric=NormalizedHammingDistance(), \n                                    clustering_method=HierarchicalClustering(0.1))\n# Collapse identical CDR3s but only within each cluster\ncollapsed = combine(groupby(result, [:d_region, :lineage_id, :j_call_first, :v_call_first, :cdr3]), :cdr3 => length => :count)\n\n# Generate diagnostic plots (requires CairoMakie)\n# using CairoMakie\n# plot_diagnostics(collapsed_df)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For more detailed information on each function and its options, please refer to the API documentation below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LineageCollapse]","category":"page"},{"location":"#LineageCollapse.HierarchicalClustering","page":"Home","title":"LineageCollapse.HierarchicalClustering","text":"HierarchicalClustering(cutoff::Float64)\n\nA type representing hierarchical clustering with a cutoff.\n\nArguments\n\ncutoff::Float64: The cutoff value for the clustering, below which clusters are merged. Higher values result in fewer clusters.\n\n\n\n\n\n","category":"type"},{"location":"#LineageCollapse.compute_distance-Tuple{HammingDistance, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}","page":"Home","title":"LineageCollapse.compute_distance","text":"compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64\n\nCompute the distance between two LongDNA{4} sequences using the specified distance metric.\n\n\n\n\n\n","category":"method"},{"location":"#LineageCollapse.compute_distance-Tuple{LevenshteinDistance, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}","page":"Home","title":"LineageCollapse.compute_distance","text":"compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64\n\nCompute the distance between two LongDNA{4} sequences using the specified distance metric.\n\n\n\n\n\n","category":"method"},{"location":"#LineageCollapse.compute_distance-Tuple{NormalizedHammingDistance, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}, BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}","page":"Home","title":"LineageCollapse.compute_distance","text":"compute_distance(metric::DistanceMetric, x::LongDNA{4}, y::LongDNA{4})::Float64\n\nCompute the distance between two LongDNA{4} sequences using the specified distance metric.\n\n\n\n\n\n","category":"method"},{"location":"#LineageCollapse.compute_pairwise_distance-Tuple{Union{DistanceMetric, LineageCollapse.NormalizedDistanceMetric}, Vector{BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}}","page":"Home","title":"LineageCollapse.compute_pairwise_distance","text":"compute_pairwise_distance(metric::Union{DistanceMetric, NormalizedDistanceMetric}, sequences::Vector{LongDNA{4}})::Matrix{Float64}\n\nCompute pairwise distances between sequences using the specified distance metric.\n\n\n\n\n\n","category":"method"},{"location":"#LineageCollapse.deduplicate_data","page":"Home","title":"LineageCollapse.deduplicate_data","text":"deduplicate_data(df::DataFrame, use_barcode::Bool=false)::DataFrame\n\nDeduplicate the input DataFrame based on sequence or sequence+barcode.\n\nArguments\n\ndf::DataFrame: Input DataFrame.\nuse_barcode::Bool=false: Whether to use barcode for deduplication.\n\nReturns\n\nDataFrame: Deduplicated DataFrame.\n\n\n\n\n\n","category":"function"},{"location":"#LineageCollapse.load_data-Tuple{String}","page":"Home","title":"LineageCollapse.load_data","text":"load_data(filepath::String; \n          delimiter::Char='\t', \n          required_columns=[:sequence_id, :sequence, :v_sequence_end, :j_sequence_start, :cdr3, :v_call, :j_call, :stop_codon])::DataFrame\n\nLoad data from a file (compressed or uncompressed) and return a DataFrame.\n\nArguments\n\nfilepath::String: Path to the data file.\ndelimiter::Char='\t': Delimiter used in the data file (default: tab).\nrequired_columns::Vector{Symbol}: Required columns to select from the data file.\n\nReturns\n\nDataFrame: DataFrame containing the loaded data.\n\nThrows\n\nArgumentError: If any of the required columns are missing in the data file.\n\n\n\n\n\n","category":"method"},{"location":"#LineageCollapse.perform_clustering-Tuple{HierarchicalClustering, Symbol, Matrix{Float64}}","page":"Home","title":"LineageCollapse.perform_clustering","text":"perform_clustering(method::HierarchicalClustering, linkage::Symbol, dist_matrix::Matrix{Float64})::Vector{Int}\n\nPerform hierarchical clustering on the distance matrix using the specified method and linkage.\n\n\n\n\n\n","category":"method"},{"location":"#LineageCollapse.preprocess_data-Tuple{DataFrames.DataFrame}","page":"Home","title":"LineageCollapse.preprocess_data","text":"preprocess_data(df::DataFrame; min_d_region_length::Union{Int,Nothing}=nothing, deduplicate::Bool=false, use_barcode::Bool=false)::DataFrame\n\nPreprocess the input DataFrame by performing data cleaning and transformation.\n\nArguments\n\ndf::DataFrame: Input DataFrame.\nmin_d_region_length::Union{Int,Nothing}=nothing: Minimum length of the D region to keep. If nothing, no filtering is applied.\ndeduplicate::Bool=false: Whether to deduplicate the DataFrame.\nuse_barcode::Bool=false: Whether to use barcode for deduplication (only applicable if deduplicate is true).\n\nReturns\n\nDataFrame: Preprocessed DataFrame.\n\n\n\n\n\n","category":"method"},{"location":"#LineageCollapse.process_lineages-Tuple{DataFrames.DataFrame}","page":"Home","title":"LineageCollapse.process_lineages","text":"process_lineages(df::DataFrame; \n                distance_metric::Union{DistanceMetric, NormalizedDistanceMetric} = NormalizedHammingDistance(),\n                clustering_method::ClusteringMethod = HierarchicalClustering(0.1),\n                linkage::Symbol = :single)::DataFrame\n\nProcess lineages from a DataFrame of CDR3 sequences.\n\n\n\n\n\n","category":"method"}]
}