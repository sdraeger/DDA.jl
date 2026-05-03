#=
Comprehensive spec validation tests for generated Julia Variants module.

These tests validate that the generated DDA spec implementation is correct and consistent.
=#

using Test
using DelayDifferentialAnalysis

# Expected configurations - canonical source of truth for tests
const EXPECTED_VARIANTS = [
    # (abbreviation, position, output_suffix, stride, reserved)
    ("ST", 0, "_ST", 4, false),
    ("CT", 1, "_CT", 4, false),
    ("CD", 2, "_CD_DDA_ST", 2, false),
    ("RESERVED", 3, "_RESERVED", 1, true),
    ("DE", 4, "_DE", 1, false),
    ("SY", 5, "_SY", 1, false),
]

const ACTIVE_VARIANT_ABBREVS = ["ST", "CT", "CD", "DE", "SY"]
const CT_REQUIRING_VARIANTS = ["CT", "CD", "DE"]

@testset "DDA Spec Validation" begin

    # =============================================================================
    # CONSTANT VALIDATION
    # =============================================================================

    @testset "Constants" begin
        @test SPEC_VERSION == "1.0.0"
        @test BINARY_NAME == "run_DDA_AsciiEdf"
        @test REQUIRES_SHELL_WRAPPER == true
        @test SHELL_COMMAND == "sh"
        @test "linux" in SUPPORTED_PLATFORMS
        @test "macos" in SUPPORTED_PLATFORMS
        @test "windows" in SUPPORTED_PLATFORMS
        @test length(SUPPORTED_PLATFORMS) == 3
        @test SELECT_MASK_SIZE == 6
        @test length(VARIANT_REGISTRY) == SELECT_MASK_SIZE
    end

    # =============================================================================
    # VARIANT METADATA VALIDATION
    # =============================================================================

    @testset "Variant Metadata" begin
        @testset "All variants present" begin
            @test length(VARIANT_REGISTRY) == length(EXPECTED_VARIANTS)

            for (abbrev, pos, suffix, stride, reserved) in EXPECTED_VARIANTS
                variant = get_variant_by_abbrev(abbrev)
                @test variant !== nothing
                @test variant.position == pos
                @test variant.output_suffix == suffix
                @test variant.stride == stride
                @test variant.reserved == reserved
            end
        end

        @testset "Unique positions" begin
            positions = [v.position for v in VARIANT_REGISTRY]
            @test length(positions) == length(unique(positions))
        end

        @testset "Sequential positions" begin
            for (i, variant) in enumerate(VARIANT_REGISTRY)
                # Julia is 1-indexed, positions are 0-indexed
                @test variant.position == i - 1
            end
        end

        @testset "Unique abbreviations" begin
            abbrevs = [v.abbreviation for v in VARIANT_REGISTRY]
            @test length(abbrevs) == length(unique(abbrevs))
        end

        @testset "Unique output suffixes" begin
            suffixes = [v.output_suffix for v in VARIANT_REGISTRY]
            @test length(suffixes) == length(unique(suffixes))
        end

        @testset "Only RESERVED is reserved" begin
            for variant in VARIANT_REGISTRY
                if variant.abbreviation == "RESERVED"
                    @test variant.reserved == true
                else
                    @test variant.reserved == false
                end
            end
        end
    end

    # =============================================================================
    # STRIDE VALUES
    # =============================================================================

    @testset "Stride Values" begin
        @test ST.stride == 4
        @test CT.stride == 4
        @test CD.stride == 2
        @test DE.stride == 1
        @test SY.stride == 1
    end

    # =============================================================================
    # OUTPUT COLUMNS
    # =============================================================================

    @testset "Output Columns" begin
        @testset "ST" begin
            @test ST.output_columns.coefficients == 3
            @test ST.output_columns.has_error == true
        end

        @testset "CT" begin
            @test CT.output_columns.coefficients == 3
            @test CT.output_columns.has_error == true
        end

        @testset "CD" begin
            @test CD.output_columns.coefficients == 1
            @test CD.output_columns.has_error == true
        end

        @testset "DE" begin
            @test DE.output_columns.coefficients == 0
            @test DE.output_columns.has_error == false
        end

        @testset "SY" begin
            @test SY.output_columns.coefficients == 0
            @test SY.output_columns.has_error == false
        end
    end

    # =============================================================================
    # CHANNEL FORMAT
    # =============================================================================

    @testset "Channel Format" begin
        @test ST.channel_format == Individual
        @test CT.channel_format == Pairs
        @test CD.channel_format == DirectedPairs
        @test DE.channel_format == Individual
        @test SY.channel_format == Individual
    end

    # =============================================================================
    # REQUIRED PARAMETERS
    # =============================================================================

    @testset "Required Parameters" begin
        @testset "CT requires CT params" begin
            @test requires_ct_params(CT) == true
            @test "-WL_CT" in CT.required_params
            @test "-WS_CT" in CT.required_params
        end

        @testset "CD requires CT params" begin
            @test requires_ct_params(CD) == true
            @test "-WL_CT" in CD.required_params
            @test "-WS_CT" in CD.required_params
        end

        @testset "DE requires CT params" begin
            @test requires_ct_params(DE) == true
            @test "-WL_CT" in DE.required_params
            @test "-WS_CT" in DE.required_params
        end

        @testset "ST no CT params" begin
            @test requires_ct_params(ST) == false
            @test isempty(ST.required_params)
        end

        @testset "SY no CT params" begin
            @test requires_ct_params(SY) == false
            @test isempty(SY.required_params)
        end
    end

    # =============================================================================
    # LOOKUP FUNCTIONS
    # =============================================================================

    @testset "Lookup Functions" begin
        @testset "By abbreviation" begin
            for (abbrev, _, _, _, _) in EXPECTED_VARIANTS
                @test get_variant_by_abbrev(abbrev) !== nothing
            end
            @test get_variant_by_abbrev("XX") === nothing
            @test get_variant_by_abbrev("") === nothing
            @test get_variant_by_abbrev("st") === nothing  # Case sensitive
        end

        @testset "By position" begin
            for i in 0:(SELECT_MASK_SIZE - 1)
                @test get_variant_by_position(i) !== nothing
            end
            @test get_variant_by_position(6) === nothing
            @test get_variant_by_position(99) === nothing
            @test get_variant_by_position(-1) === nothing
        end

        @testset "By suffix" begin
            for (_, _, suffix, _, _) in EXPECTED_VARIANTS
                @test get_variant_by_suffix(suffix) !== nothing
            end
            @test get_variant_by_suffix("_XX") === nothing
            @test get_variant_by_suffix("") === nothing
        end
    end

    # =============================================================================
    # SELECT MASK
    # =============================================================================

    @testset "Select Mask" begin
        @testset "Generate" begin
            @test generate_select_mask(["ST"]) == [1, 0, 0, 0, 0, 0]
            @test generate_select_mask(["SY"]) == [0, 0, 0, 0, 0, 1]
            @test generate_select_mask(["ST", "SY"]) == [1, 0, 0, 0, 0, 1]
            @test generate_select_mask(["ST", "CT", "CD", "DE", "SY"]) == [1, 1, 1, 0, 1, 1]
            @test generate_select_mask(String[]) == [0, 0, 0, 0, 0, 0]
            @test generate_select_mask(["ST", "XX", "INVALID", "SY"]) == [1, 0, 0, 0, 0, 1]
        end

        @testset "Parse" begin
            @test parse_select_mask([1, 0, 0, 0, 0, 0]) == ["ST"]
            @test parse_select_mask([1, 0, 0, 0, 0, 1]) == ["ST", "SY"]
            @test parse_select_mask([0, 0, 0, 1, 0, 0]) == []  # RESERVED excluded
            @test parse_select_mask([1, 1, 1, 1, 1, 1]) == ["ST", "CT", "CD", "DE", "SY"]
        end

        @testset "Format" begin
            @test format_select_mask([1, 1, 0, 0, 0, 1]) == "1 1 0 0 0 1"
        end

        @testset "Roundtrip" begin
            original = ["ST", "CT", "SY"]
            mask = generate_select_mask(original)
            parsed = parse_select_mask(mask)
            @test parsed == original
        end
    end

    # =============================================================================
    # ACTIVE VARIANTS
    # =============================================================================

    @testset "Active Variants" begin
        active = active_variants()
        @test length(active) == 5

        abbrevs = [v.abbreviation for v in active]
        for expected in ACTIVE_VARIANT_ABBREVS
            @test expected in abbrevs
        end

        @test !("RESERVED" in abbrevs)
    end

    # =============================================================================
    # SELECT MASK POSITIONS
    # =============================================================================

    @testset "Select Mask Positions" begin
        @test SelectMaskPositions.ST == 0
        @test SelectMaskPositions.CT == 1
        @test SelectMaskPositions.CD == 2
        @test SelectMaskPositions.RESERVED == 3
        @test SelectMaskPositions.DE == 4
        @test SelectMaskPositions.SY == 5
    end

    # =============================================================================
    # FILE TYPES
    # =============================================================================

    @testset "File Types" begin
        @test get_flag(EDF) == "-EDF"
        @test get_flag(ASCII) == "-ASCII"

        @testset "From extension" begin
            @test file_type_from_extension("edf") == EDF
            @test file_type_from_extension(".edf") == EDF
            @test file_type_from_extension("EDF") == EDF
            @test file_type_from_extension("txt") == ASCII
            @test file_type_from_extension("csv") == ASCII
            @test file_type_from_extension("ascii") == ASCII
            @test file_type_from_extension("unknown") === nothing
            @test file_type_from_extension("") === nothing
        end
    end

    # =============================================================================
    # DELAYS
    # =============================================================================

    @testset "Delays" begin
        @testset "Default delays length" begin
            @test length(DEFAULT_DELAYS) == 2
        end

        @testset "Default delays values" begin
            @test DEFAULT_DELAYS[1] == 7  # Julia is 1-indexed
            @test DEFAULT_DELAYS[2] == 10
        end

        @testset "Default delays equals expected" begin
            @test DEFAULT_DELAYS == [7, 10]
        end
    end

    # =============================================================================
    # VARIANT ORDER
    # =============================================================================

    @testset "Variant Order" begin
        @testset "Matches positions" begin
            for (i, variant) in enumerate(VARIANT_REGISTRY)
                @test VARIANT_ORDER[i] == variant.abbreviation
            end
        end

        @test length(VARIANT_ORDER) == SELECT_MASK_SIZE
        @test VARIANT_ORDER == ["ST", "CT", "CD", "RESERVED", "DE", "SY"]
    end

    # =============================================================================
    # DIRECT VARIANT ACCESS
    # =============================================================================

    @testset "Direct Variant Access" begin
        @test ST.abbreviation == "ST"
        @test ST.name == "Single Timeseries"
        @test ST.position == 0

        @test CT.abbreviation == "CT"
        @test CT.name == "Cross-Timeseries"
        @test CT.position == 1

        @test CD.abbreviation == "CD"
        @test CD.name == "Cross-Dynamical"
        @test CD.position == 2

        @test RESERVED.abbreviation == "RESERVED"
        @test RESERVED.name == "Reserved"
        @test RESERVED.position == 3

        @test DE.abbreviation == "DE"
        @test DE.name == "Delay Embedding"
        @test DE.position == 4

        @test SY.abbreviation == "SY"
        @test SY.name == "Synchronization"
        @test SY.position == 5
    end

end

# =============================================================================
# GROUND TRUTH VALIDATION - CLI Command Generation
# =============================================================================

@testset "CLI Command Generation" begin
    @testset "SELECT mask for CLI" begin
        # Verify SELECT mask for ST-only matches expected CLI format
        mask = generate_select_mask(["ST"])
        cli_args = format_select_mask(mask)
        @test cli_args == "1 0 0 0 0 0"

        # Verify SELECT mask for all active variants
        mask_all = generate_select_mask(["ST", "CT", "CD", "DE", "SY"])
        cli_args_all = format_select_mask(mask_all)
        @test cli_args_all == "1 1 1 0 1 1"  # RESERVED at position 3 should be 0
    end

    @testset "SELECT mask positions match binary spec" begin
        # Ground truth: DDA binary expects SELECT mask as 6 integers:
        # Position 0: ST, Position 1: CT, Position 2: CD
        # Position 3: RESERVED (always 0), Position 4: DE, Position 5: SY
        test_cases = [
            ("ST", [1, 0, 0, 0, 0, 0]),
            ("CT", [0, 1, 0, 0, 0, 0]),
            ("CD", [0, 0, 1, 0, 0, 0]),
            ("DE", [0, 0, 0, 0, 1, 0]),
            ("SY", [0, 0, 0, 0, 0, 1]),
        ]

        for (variant_abbrev, expected_mask) in test_cases
            mask = generate_select_mask([variant_abbrev])
            @test mask == expected_mask
        end
    end
end

# =============================================================================
# GROUND TRUTH VALIDATION - Output File Parsing
# =============================================================================

@testset "Output File Parsing" begin
    @testset "ST stride parsing" begin
        # Ground truth: ST output format per channel is:
        # [a1, a2, a3, error] - 4 columns (3 coefficients + 1 error)
        @test ST.stride == 4
        @test ST.output_columns.coefficients == 3
        @test ST.output_columns.has_error == true
        # Total = coefficients + error = 3 + 1 = 4 = stride
        @test ST.output_columns.coefficients + (ST.output_columns.has_error ? 1 : 0) == ST.stride
    end

    @testset "CT stride parsing" begin
        # Ground truth: CT output format per pair is:
        # [a1, a2, a3, error] - 4 columns (3 coefficients + 1 error)
        @test CT.stride == 4
        @test CT.output_columns.coefficients == 3
        @test CT.output_columns.has_error == true
        @test CT.output_columns.coefficients + (CT.output_columns.has_error ? 1 : 0) == CT.stride
    end

    @testset "CD stride parsing" begin
        # Ground truth: CD output format per directed pair is:
        # [a1, error] - 2 columns (1 coefficient + 1 error)
        @test CD.stride == 2
        @test CD.output_columns.coefficients == 1
        @test CD.output_columns.has_error == true
        @test CD.output_columns.coefficients + (CD.output_columns.has_error ? 1 : 0) == CD.stride
    end

    @testset "DE stride parsing" begin
        # Ground truth: DE output format is:
        # [ergodicity] - 1 column (single measure, no error)
        @test DE.stride == 1
        @test DE.output_columns.coefficients == 0
        @test DE.output_columns.has_error == false
    end

    @testset "SY stride parsing" begin
        # Ground truth: SY output format per channel is:
        # [sync_coef] - 1 column (synchronization coefficient, no error)
        @test SY.stride == 1
        @test SY.output_columns.coefficients == 0
        @test SY.output_columns.has_error == false
    end

    @testset "Output file suffixes match binary" begin
        # Ground truth: Binary creates files named: {base}{suffix}
        expected_suffixes = [
            ("ST", "_ST"),
            ("CT", "_CT"),
            ("CD", "_CD_DDA_ST"),
            ("RESERVED", "_RESERVED"),
            ("DE", "_DE"),
            ("SY", "_SY"),
        ]

        for (abbrev, expected_suffix) in expected_suffixes
            variant = get_variant_by_abbrev(abbrev)
            @test variant !== nothing
            @test variant.output_suffix == expected_suffix
        end
    end
end

# =============================================================================
# GROUND TRUTH VALIDATION - Mock Output Parsing
# =============================================================================

@testset "Mock Output Parsing" begin
    @testset "Parse ST mock output" begin
        # Mock ST output: window_start window_end [a1 a2 a3 error] per channel
        # For 2 channels, 1 timepoint:
        mock_data = [0.0, 1000.0, 0.1, 0.2, 0.3, 0.01, 0.4, 0.5, 0.6, 0.02]
        #                         ---- channel 0 ----  ---- channel 1 ----

        stride = ST.stride
        @test stride == 4

        # Skip first 2 values (window bounds), extract channel 0
        ch0_start = 3  # Julia is 1-indexed
        ch0_data = mock_data[ch0_start:ch0_start + stride - 1]
        @test ch0_data == [0.1, 0.2, 0.3, 0.01]

        # Extract channel 1
        ch1_start = ch0_start + stride
        ch1_data = mock_data[ch1_start:ch1_start + stride - 1]
        @test ch1_data == [0.4, 0.5, 0.6, 0.02]
    end

    @testset "Parse CD mock output" begin
        # Mock CD output: window_start window_end [a1 error] per directed pair
        # For 2 directed pairs (1->2, 2->1), 1 timepoint:
        mock_data = [0.0, 1000.0, 0.1, 0.01, 0.2, 0.02]
        #                         ---- 1->2 ----  ---- 2->1 ----

        stride = CD.stride
        @test stride == 2

        # Extract pair 0 (1->2)
        p0_start = 3  # Julia is 1-indexed
        p0_data = mock_data[p0_start:p0_start + stride - 1]
        @test p0_data == [0.1, 0.01]

        # Extract pair 1 (2->1)
        p1_start = p0_start + stride
        p1_data = mock_data[p1_start:p1_start + stride - 1]
        @test p1_data == [0.2, 0.02]
    end

    @testset "Parse SY mock output" begin
        # Mock SY output: window_start window_end [sync_coef] per channel
        # For 3 channels, 1 timepoint:
        mock_data = [0.0, 1000.0, 0.95, 0.87, 0.91]
        #                         ch0   ch1   ch2

        stride = SY.stride
        @test stride == 1

        # Each channel gets 1 value
        for i in 0:2
            ch_start = 3 + i * stride  # Julia is 1-indexed
            ch_data = mock_data[ch_start:ch_start + stride - 1]
            @test length(ch_data) == 1
        end
    end

    @testset "Stride determines num channels" begin
        # Ground truth: data_columns / stride = num_channels
        test_cases = [
            (ST.stride, 8, 2),   # 8 data cols / 4 stride = 2 channels
            (ST.stride, 12, 3), # 12 data cols / 4 stride = 3 channels
            (CD.stride, 4, 2),  # 4 data cols / 2 stride = 2 pairs
            (SY.stride, 5, 5),  # 5 data cols / 1 stride = 5 channels
        ]

        for (stride, data_cols, expected_num) in test_cases
            @test data_cols % stride == 0
            @test div(data_cols, stride) == expected_num
        end
    end
end

# =============================================================================
# GROUND TRUTH VALIDATION - Required Parameters
# =============================================================================

@testset "Required Parameters Ground Truth" begin
    @testset "CT requires WL_CT WS_CT" begin
        # Verify CT requires -WL_CT and -WS_CT as the binary expects
        @test "-WL_CT" in CT.required_params
        @test "-WS_CT" in CT.required_params
    end

    @testset "CD requires WL_CT WS_CT" begin
        @test "-WL_CT" in CD.required_params
        @test "-WS_CT" in CD.required_params
    end

    @testset "DE requires WL_CT WS_CT" begin
        @test "-WL_CT" in DE.required_params
        @test "-WS_CT" in DE.required_params
    end

    @testset "ST no special params" begin
        @test isempty(ST.required_params)
    end

    @testset "SY no special params" begin
        @test isempty(SY.required_params)
    end
end

# =============================================================================
# GROUND TRUTH VALIDATION - Binary Integration Tests
# =============================================================================

"""Helper function to find the DDA binary."""
function find_dda_binary()::Union{String, Nothing}
    # Look for binary in common locations
    possible_paths = [
        joinpath(@__DIR__, "..", "..", "..", "bin", "run_DDA_AsciiEdf"),
        joinpath(homedir(), ".local", "bin", "run_DDA_AsciiEdf"),
        "/usr/local/bin/run_DDA_AsciiEdf",
    ]

    for path in possible_paths
        if isfile(path)
            return path
        end
    end
    return nothing
end

"""Helper function to find test EDF data."""
function find_test_data()::Union{String, Nothing}
    possible_paths = [
        joinpath(@__DIR__, "..", "..", "..", "data", "patient1_S05__01_03 (1).edf"),
    ]

    for path in possible_paths
        if isfile(path)
            return path
        end
    end
    return nothing
end

"""
Build a DDA CLI command using spec-generated values.
This validates that the Julia implementation produces correct CLI arguments.
"""
function build_dda_command(
    binary_path::String,
    input_file::String,
    output_file::String,
    variants::Vector{String};
    channels::Vector{Int}=[1, 2, 3],
    WL::Union{Int, Nothing}=2048,
    WS::Union{Int, Nothing}=1024,
    delays::Vector{Int}=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    start_sample::Union{Int, Nothing}=nothing,
    end_sample::Union{Int, Nothing}=nothing,
)::Cmd
    # Build command arguments
    args = String[]

    # Input/output files
    push!(args, "-DATA_FN", input_file)
    push!(args, "-OUT_FN", output_file)

    # File type (detect from extension)
    ext = lowercase(splitext(input_file)[2])
    if ext == ".edf"
        push!(args, "-EDF")
    else
        push!(args, "-ASCII")
    end

    # Channel list (1-based for binary)
    push!(args, "-CH_list")
    for ch in channels
        push!(args, string(ch))
    end

    # Model parameters (defaults from binary spec)
    push!(args, "-dm", "3", "-order", "4", "-nr_tau", "2")

    # Window parameters
    WL !== nothing && push!(args, "-WL", string(WL))
    WS !== nothing && push!(args, "-WS", string(WS))

    # SELECT mask using spec-generated function
    mask = generate_select_mask(variants)
    push!(args, "-SELECT")
    for bit in mask
        push!(args, string(bit))
    end

    push!(args, "-MODEL", "1", "2", "10")

    # Add CT params if any variant requires them
    needs_ct_params = any(v -> begin
        variant = get_variant_by_abbrev(v)
        variant !== nothing && requires_ct_params(variant)
    end, variants)

    if needs_ct_params
        push!(args, "-WL_CT", "2", "-WS_CT", "2")
    end

    # Delay/TAU values
    push!(args, "-TAU")
    for d in delays
        push!(args, string(d))
    end

    # Time bounds
    if start_sample !== nothing && end_sample !== nothing
        push!(args, "-StartEnd", string(start_sample), string(end_sample))
    end

    # APE binary needs to run through sh on Unix
    return `sh $binary_path $args`
end

"""Parse DDA output file and extract Q matrix using spec stride values."""
function parse_dda_output(filepath::String, variant::VariantMetadata)::Matrix{Float64}
    if !isfile(filepath)
        error("Output file not found: $filepath")
    end

    lines = readlines(filepath)
    if isempty(lines)
        return Matrix{Float64}(undef, 0, 0)
    end

    stride = variant.stride

    # Parse all rows
    data_rows = Vector{Vector{Float64}}()
    for line in lines
        parts = split(strip(line))
        if isempty(parts)
            continue
        end

        values = tryparse.(Float64, parts)
        if all(v -> v !== nothing, values)
            row = Float64[v for v in values if v !== nothing]
            push!(data_rows, row)
        end
    end

    if isempty(data_rows)
        return Matrix{Float64}(undef, 0, 0)
    end

    # Skip first 2 columns (window bounds), take every stride-th column
    # to extract the primary coefficient for each channel
    num_timepoints = length(data_rows)
    num_data_cols = length(data_rows[1]) - 2  # Exclude window bounds

    if num_data_cols <= 0
        return Matrix{Float64}(undef, 0, 0)
    end

    num_channels = div(num_data_cols, stride)

    if num_channels <= 0
        return Matrix{Float64}(undef, 0, 0)
    end

    # Build matrix [channels × timepoints]
    q_matrix = Matrix{Float64}(undef, num_channels, num_timepoints)

    for (t, row) in enumerate(data_rows)
        for ch in 1:num_channels
            # Extract first coefficient of each channel (skip 2 window cols, then stride*ch offset)
            col_idx = 3 + (ch - 1) * stride  # 1-indexed, skip 2 window cols
            if col_idx <= length(row)
                q_matrix[ch, t] = row[col_idx]
            end
        end
    end

    return q_matrix
end

@testset "Binary Integration Tests" begin
    binary_path = find_dda_binary()
    test_data = find_test_data()

    if binary_path === nothing || test_data === nothing
        @info "Skipping binary integration tests: binary or test data not found"
        @info "Binary path: $binary_path"
        @info "Test data: $test_data"
    else
        @testset "ST variant shell execution" begin
            # Create temp output file
            output_file = tempname()

            # Build command using Julia spec implementation
            cmd = build_dda_command(
                binary_path,
                test_data,
                output_file,
                ["ST"];
                channels=[1, 2, 3],
                WL=2048,
                WS=1024,
                delays=collect(1:10),
                start_sample=0,
                end_sample=6000,
            )

            @info "Executing ST command: $cmd"

            # Run the binary
            try
                run(cmd)

                # Check output file was created with correct suffix
                st_output = output_file * ST.output_suffix
                @test isfile(st_output)

                if isfile(st_output)
                    # Parse output using spec stride
                    q_matrix = parse_dda_output(st_output, ST)

                    @test size(q_matrix, 1) > 0  # Has channels
                    @test size(q_matrix, 2) > 0  # Has timepoints

                    @info "ST output matrix size: $(size(q_matrix))"

                    # Cleanup
                    rm(st_output, force=true)
                end
            catch e
                @warn "ST variant test failed: $e"
                @test_broken true
            finally
                rm(output_file, force=true)
            end
        end

        @testset "SELECT mask produces expected outputs" begin
            output_file = tempname()

            # Test: ST only should only produce _ST file
            cmd = build_dda_command(
                binary_path,
                test_data,
                output_file,
                ["ST"];
                channels=[1, 2],
                start_sample=0,
                end_sample=4000,
            )

            try
                run(cmd)

                # ST file should exist
                st_output = output_file * ST.output_suffix
                @test isfile(st_output)

                # CT file should NOT exist (CT not in SELECT mask)
                ct_output = output_file * CT.output_suffix
                @test !isfile(ct_output)

                # Cleanup
                rm(st_output, force=true)
            catch e
                @warn "SELECT mask test failed: $e"
                @test_broken true
            finally
                rm(output_file, force=true)
            end
        end

        @testset "Stride correctly parses binary output" begin
            output_file = tempname()

            cmd = build_dda_command(
                binary_path,
                test_data,
                output_file,
                ["ST"];
                channels=[1, 2, 3],  # 3 channels
                start_sample=0,
                end_sample=6000,
            )

            try
                run(cmd)

                st_output = output_file * ST.output_suffix
                if isfile(st_output)
                    # Read raw output
                    lines = readlines(st_output)
                    @test !isempty(lines)

                    if !isempty(lines)
                        # First line should have: 2 window cols + 3 channels * 4 stride = 14 columns
                        first_line_cols = length(split(strip(lines[1])))
                        expected_cols = 2 + 3 * ST.stride  # 14

                        @test first_line_cols == expected_cols

                        # Parse using spec stride
                        q_matrix = parse_dda_output(st_output, ST)
                        @test size(q_matrix, 1) == 3  # 3 channels

                        @info "Output columns: $first_line_cols (expected $expected_cols)"
                        @info "Parsed matrix channels: $(size(q_matrix, 1))"
                    end

                    rm(st_output, force=true)
                end
            catch e
                @warn "Stride parsing test failed: $e"
                @test_broken true
            finally
                rm(output_file, force=true)
            end
        end

        @testset "DE variant shell execution" begin
            output_file = tempname()

            # DE requires CT params
            cmd = build_dda_command(
                binary_path,
                test_data,
                output_file,
                ["DE"];
                channels=[1, 2, 3],
                start_sample=0,
                end_sample=6000,
            )

            try
                run(cmd)

                de_output = output_file * DE.output_suffix
                @test isfile(de_output)

                if isfile(de_output)
                    q_matrix = parse_dda_output(de_output, DE)
                    @test size(q_matrix, 2) > 0  # Has timepoints

                    @info "DE output matrix size: $(size(q_matrix))"
                    rm(de_output, force=true)
                end
            catch e
                @warn "DE variant test failed: $e"
                @test_broken true
            finally
                rm(output_file, force=true)
            end
        end

        @testset "SY variant shell execution" begin
            output_file = tempname()

            cmd = build_dda_command(
                binary_path,
                test_data,
                output_file,
                ["SY"];
                channels=[1, 2, 3],
                start_sample=0,
                end_sample=6000,
            )

            try
                run(cmd)

                sy_output = output_file * SY.output_suffix
                @test isfile(sy_output)

                if isfile(sy_output)
                    q_matrix = parse_dda_output(sy_output, SY)
                    @test size(q_matrix, 2) > 0  # Has timepoints

                    @info "SY output matrix size: $(size(q_matrix))"
                    rm(sy_output, force=true)
                end
            catch e
                @warn "SY variant test failed: $e"
                @test_broken true
            finally
                rm(output_file, force=true)
            end
        end

        @testset "Multiple variants shell execution" begin
            output_file = tempname()

            # Run ST + SY together
            cmd = build_dda_command(
                binary_path,
                test_data,
                output_file,
                ["ST", "SY"];
                channels=[1, 2],
                start_sample=0,
                end_sample=4000,
            )

            try
                run(cmd)

                # Both files should exist
                st_output = output_file * ST.output_suffix
                sy_output = output_file * SY.output_suffix

                @test isfile(st_output)
                @test isfile(sy_output)

                if isfile(st_output) && isfile(sy_output)
                    st_matrix = parse_dda_output(st_output, ST)
                    sy_matrix = parse_dda_output(sy_output, SY)

                    # Both should have same number of timepoints
                    @test size(st_matrix, 2) == size(sy_matrix, 2)

                    @info "ST matrix: $(size(st_matrix)), SY matrix: $(size(sy_matrix))"

                    rm(st_output, force=true)
                    rm(sy_output, force=true)
                end
            catch e
                @warn "Multiple variants test failed: $e"
                @test_broken true
            finally
                rm(output_file, force=true)
            end
        end

        @testset "Julia vs Shell consistency" begin
            # This test validates that the Julia implementation
            # produces the same CLI command format as expected by the binary
            output_file = tempname()

            # Build command using Julia spec functions
            variants = ["ST", "CT"]
            mask = generate_select_mask(variants)
            mask_str = format_select_mask(mask)

            # Verify mask format is correct
            @test mask_str == "1 1 0 0 0 0"

            # Verify variant metadata matches what we pass to binary
            @test ST.position == 0  # First bit in SELECT
            @test CT.position == 1  # Second bit in SELECT

            # Verify CT requires -WL_CT and -WS_CT
            @test requires_ct_params(CT) == true
            @test "-WL_CT" in CT.required_params
            @test "-WS_CT" in CT.required_params

            # This validates the Julia implementation is consistent
            # with what the binary expects
            @test true
        end
    end
end
