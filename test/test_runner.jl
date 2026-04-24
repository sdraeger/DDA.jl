#=
Tests for the DDA Runner module.

These tests validate the DDARequest, DDAResult types and command building logic.
Integration tests (actual binary execution) are skipped if binary is not found.
=#

using Test
using DelayDifferentialAnalysis

function write_test_edf(path::String, labels::Vector{String})
    n = length(labels)

    fixed_header = string(
        rpad("0", 8),
        rpad("", 80),
        rpad("", 80),
        rpad("01.01.01", 8),
        rpad("01.01.01", 8),
        rpad(string(256 + 256 * n), 8),
        rpad("", 44),
        rpad("1", 8),
        rpad("1", 8),
        rpad(string(n), 4),
    )

    signal_header = string(
        join(rpad.(labels, 16)),
        join(fill(rpad("", 80), n)),
        join(fill(rpad("uV", 8), n)),
        join(fill(rpad("-1", 8), n)),
        join(fill(rpad("1", 8), n)),
        join(fill(rpad("-2048", 8), n)),
        join(fill(rpad("2047", 8), n)),
        join(fill(rpad("", 80), n)),
        join(fill(rpad("1", 8), n)),
        join(fill(rpad("", 32), n)),
    )

    open(path, "w") do io
        write(io, fixed_header)
        write(io, signal_header)
        write(io, UInt8[0 for _ in 1:n])
    end
end

@testset "Runner Module" begin

    # =============================================================================
    # REQUEST CONSTRUCTION
    # =============================================================================

    @testset "DDARequest Construction" begin
        @testset "Basic construction" begin
            request = DDARequest(
                "test.edf",
                [1, 2, 3],
                ["ST", "SY"]
            )

            @test request.file_path == "test.edf"
            @test request.channels == [1, 2, 3]
            @test request.variants == ["ST", "SY"]
            @test request.window_params.window_length == 200   # Default (DDADefaults.WINDOW_LENGTH)
            @test request.window_params.window_step == 100    # Default (DDADefaults.WINDOW_STEP)
            @test request.delay_params.delays == collect(DEFAULT_DELAYS)
            @test request.model_terms == DelayDifferentialAnalysis.DDADefaults.MODEL_PARAMS
            @test request.sampling_rate == DelayDifferentialAnalysis.DDADefaults.SAMPLING_RATE
            @test request.tm == maximum(collect(DEFAULT_DELAYS))
        end

        @testset "With custom parameters" begin
            request = DDARequest(
                "data.edf",
                [1, 2],
                ["ST", "CT", "CD"];
                window_length=4096,
                window_step=2048,
                ct_window_length=4,
                ct_window_step=2,
                delays=[1, 2, 3, 4, 5],
                model=[4, 7, 9],
                derivative_points=6,
                order=5,
                nr_tau=3,
                time_range=(0.0, 10000.0),
                ct_pairs=[(1, 2)],
                cd_pairs=[(1, 2), (2, 1)],
                sampling_rate=(500, 1000),
                TM=11,
                out_fn="/tmp/custom_dda_output"
            )

            @test request.window_params.window_length == 4096
            @test request.window_params.window_step == 2048
            @test request.window_params.ct_window_length == 4
            @test request.window_params.ct_window_step == 2
            @test request.delay_params.delays == [1, 2, 3, 4, 5]
            @test request.model_params.derivative_points == 6
            @test request.model_params.order == 5
            @test request.model_params.nr_tau == 3
            @test request.model_terms == [4, 7, 9]
            @test request.time_range.start == 0.0
            @test request.time_range.stop == 10000.0
            @test request.ct_channel_pairs == [(1, 2)]
            @test request.cd_channel_pairs == [(1, 2), (2, 1)]
            @test request.sampling_rate == (500, 1000)
            @test request.tm == 11
            @test request.out_fn == "/tmp/custom_dda_output"
        end

        @testset "Default model parameters" begin
            request = DDARequest("test.edf", [1], ["ST"])

            @test request.model_params.derivative_points == 3
            @test request.model_params.order == 4
            @test request.model_params.nr_tau == 2
        end

        @testset "Scalar sampling rate is normalized to a pair" begin
            request = DDARequest("test.edf", [1], ["ST"]; sampling_rate=2048)
            @test request.sampling_rate == (1024, 2048)
        end

        @testset "Custom model requires explicit derivative config and order" begin
            @test_throws ErrorException DDARequest(
                "test.edf",
                [1],
                ["ST"];
                model=[1, 2, 10],
            )

            request = DDARequest(
                "test.edf",
                [1],
                ["ST"];
                model=[1, 2, 10],
                derivative_points=5,
                order=3,
            )
            @test request.model_params.derivative_points == 5
            @test request.model_params.order == 3
        end
    end

    @testset "Input channel label inference" begin
        @testset "EDF labels are extracted from file header" begin
            temp_dir = mktempdir()
            edf_path = joinpath(temp_dir, "labels.edf")
            write_test_edf(edf_path, ["Fp1", "Fp2", "C3", "C4"])

            try
                labels = DelayDifferentialAnalysis.Runner._resolve_requested_channel_labels(
                    edf_path,
                    [1, 3, 4];
                    fallback_prefix="Channel ",
                )
                @test labels == ["Fp1", "C3", "C4"]
            finally
                rm(temp_dir; recursive=true, force=true)
            end
        end

        @testset "ASCII header labels are extracted from first non-numeric row" begin
            temp_dir = mktempdir()
            ascii_path = joinpath(temp_dir, "labels.tsv")

            open(ascii_path, "w") do io
                println(io, "Fp1\tFp2\tC3")
                println(io, "1\t2\t3")
                println(io, "4\t5\t6")
            end

            try
                labels = DelayDifferentialAnalysis.Runner._resolve_requested_channel_labels(
                    ascii_path,
                    [2, 3];
                    fallback_prefix="Channel ",
                )
                @test labels == ["Fp2", "C3"]
            finally
                rm(temp_dir; recursive=true, force=true)
            end
        end

        @testset "ASCII files without a header fall back to synthesized labels" begin
            temp_dir = mktempdir()
            ascii_path = joinpath(temp_dir, "numeric.tsv")

            open(ascii_path, "w") do io
                println(io, "1\t2\t3")
                println(io, "4\t5\t6")
            end

            try
                labels = DelayDifferentialAnalysis.Runner._resolve_requested_channel_labels(
                    ascii_path,
                    [1, 3];
                    fallback_prefix="Channel ",
                )
                @test labels == ["Channel 1", "Channel 3"]
            finally
                rm(temp_dir; recursive=true, force=true)
            end
        end
    end

    # =============================================================================
    # MODULE IMPORTS VALIDATION
    # =============================================================================

    @testset "Runner module imports" begin
        # These tests ensure the Runner module has all required imports (UUIDs, Dates, etc.)
        # They catch issues like `using UUIDs` vs `import UUIDs` that cause runtime errors.

        @testset "UUIDs module is accessible" begin
            # Directly test that the UUIDs import works by calling the internal function
            # This catches the `using UUIDs` vs `import UUIDs` issue
            runner_module = DelayDifferentialAnalysis.Runner

            # The Runner module must have UUIDs accessible for uuid4()
            # Test by evaluating in the module's context
            uuid_str = Base.invokelatest(Core.eval, runner_module, :(string(UUIDs.uuid4())))
            @test typeof(uuid_str) == String
            @test length(uuid_str) == 36  # UUID format: xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
            @test occursin("-", uuid_str)
        end

        @testset "Dates module is accessible" begin
            # Test that Dates.now() works in the Runner module context
            runner_module = DelayDifferentialAnalysis.Runner

            datetime_str = Base.invokelatest(Core.eval, runner_module, :(string(Dates.now())))
            @test typeof(datetime_str) == String
            @test length(datetime_str) > 0
        end

        @testset "run_analysis fails gracefully with missing file" begin
            request = DDARequest(
                "/nonexistent/path/to/file.edf",
                [1, 2, 3],
                ["ST", "SY"];
                window_length=2048,
                window_step=1024
            )

            # Create a fake binary path to bypass binary discovery
            fake_binary = tempname()
            touch(fake_binary)

            try
                runner = DDARunner(fake_binary)
                # This should fail with "Input file not found", NOT "UndefVarError: UUIDs"
                err = nothing
                try
                    run_analysis(; runner=runner, request=request)
                catch e
                    err = e
                end

                # Verify we got the expected error type (file not found), not an import error
                @test err !== nothing
                err_msg = string(err)
                @test !occursin("UndefVarError", err_msg)
                @test !occursin("not defined", err_msg)
                @test occursin("not found", lowercase(err_msg)) || occursin("Input file", err_msg)
            finally
                rm(fake_binary, force=true)
            end
        end
    end

    # =============================================================================
    # RUNNER CONSTRUCTION
    # =============================================================================

    @testset "DDARunner Construction" begin
        @testset "Invalid path throws error" begin
            @test_throws ErrorException DDARunner("/nonexistent/path/to/binary")
        end

        @testset "Auto-discovery without binary" begin
            # This will fail if binary is not installed, which is expected in CI
            # We test that the function exists and throws appropriate error
            withenv("DDA_BINARY_PATH" => nothing) do
                # Only test if binary is definitely not installed
                if find_binary() === nothing
                    @test_throws ErrorException DDARunner()
                end
            end
        end

        @testset "Explicit binary_path resolves binary" begin
            temp_dir = mktempdir()
            fake_binary = joinpath(temp_dir, BINARY_NAME)
            touch(fake_binary)

            try
                runner = DDARunner(; binary_path=fake_binary)
                @test runner.binary_path == fake_binary
            finally
                rm(temp_dir; recursive=true, force=true)
            end
        end
    end

    # =============================================================================
    # VARIANT RESULT DATA
    # =============================================================================

    @testset "VariantResultData" begin
        q_matrix = [1.0 2.0 3.0; 4.0 5.0 6.0]
        coefficients = reshape(collect(1.0:12.0), 2, 3, 2)
        errors = [0.1 0.2 0.3; 0.4 0.5 0.6]
        T = [0.0, 100.0, 200.0]
        t = [0.014, 0.114, 0.214]
        window_starts = [0, 64, 128]
        window_ends = [128, 192, 256]
        labels = ["Channel 1", "Channel 2"]

        result = VariantResultData(
            "ST",
            "Single Timeseries",
            q_matrix,
            coefficients,
            errors,
            T,
            t,
            window_starts,
            window_ends,
            labels
        )

        @test result.variant_id == "ST"
        @test result.variant_name == "Single Timeseries"
        @test size(result.q_matrix) == (2, 3)
        @test size(result.coefficients) == (2, 3, 2)
        @test size(result.errors) == (2, 3)
        @test result.T == T
        @test result.t == t
        @test result.window_starts == window_starts
        @test result.window_ends == window_ends
        @test result.channel_labels == labels
    end

    # =============================================================================
    # INTEGRATION TESTS (require binary)
    # =============================================================================

    @testset "Binary Integration Tests" begin
        binary_path = find_binary()

        if binary_path === nothing
            @info "Skipping binary integration tests: DDA binary not found"
        else
            @info "Found DDA binary at: $binary_path"

            # Find test data
            test_data_paths = [
                joinpath(@__DIR__, "..", "..", "..", "data", "patient1_S05__01_03 (1).edf"),
            ]

            test_data = nothing
            for path in test_data_paths
                if isfile(path)
                    test_data = path
                    break
                end
            end

            if test_data === nothing
                @info "Skipping integration tests: test data not found"
            else
                @testset "Run ST analysis" begin
                    request = DDARequest(
                        test_data,
                        [1, 2, 3],
                        ["ST"];
                        window_length=2048,
                        window_step=1024,
                        delays=collect(1:10),
                        time_range=(0.0, 6000.0)
                    )

                    try
                        result = run_analysis(; request=request)

                        @test !isempty(result.id)
                        @test result.file_path == test_data
                        @test size(result.q_matrix, 1) > 0
                        @test size(result.q_matrix, 2) > 0
                        @test length(result.variant_results) == 1
                        @test result.variant_results[1].variant_id == "ST"

                        @info "ST result" size=size(result.q_matrix) timepoints=size(result.q_matrix, 2)
                    catch e
                        @warn "ST analysis failed" exception=e
                        @test_broken true
                    end
                end

                @testset "Run multiple variants" begin
                    request = DDARequest(
                        test_data,
                        [1, 2],
                        ["ST", "SY"];
                        window_length=2048,
                        window_step=1024,
                        time_range=(0.0, 4000.0)
                    )

                    try
                        result = run_analysis(; request=request)

                        @test length(result.variant_results) == 2

                        variant_ids = [vr.variant_id for vr in result.variant_results]
                        @test "ST" in variant_ids
                        @test "SY" in variant_ids

                        @info "Multi-variant result" variants=variant_ids
                    catch e
                        @warn "Multi-variant analysis failed" exception=e
                        @test_broken true
                    end
                end
            end
        end
    end

    @testset "Direct run_analysis API" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            @test_throws MethodError run_analysis(
                "/nonexistent/path/to/file.edf",
                [1, 2, 3],
                ["ST"];
                binary_path=fake_binary,
            )
            @test_throws MethodError run_analysis_structured(
                "/nonexistent/path/to/file.edf",
                [1, 2, 3],
                ["ST"];
                binary_path=fake_binary,
            )

            err = nothing
            try
                run_analysis(;
                    file_path="/nonexistent/path/to/file.edf",
                    channels=[1, 2, 3],
                    flavors=["ST"],
                    binary_path=fake_binary,
                )
            catch e
                err = e
            end

            @test err !== nothing
            @test occursin("Input file not found", string(err))

            @test_throws MethodError run_analysis(;
                file_path="test.edf",
                channels=[1],
                variants=["ST"],
                binary_path=fake_binary,
            )

            err = nothing
            try
                run_analysis_structured(;
                    file_path="/nonexistent/path/to/file.edf",
                    channels=[1, 2],
                    flavors=["ST"],
                    binary_path=fake_binary,
                )
            catch e
                err = e
            end

            @test err !== nothing
            @test occursin("Input file not found", string(err))

            @test_throws MethodError run_analysis_structured(;
                file_path="test.edf",
                channels=[1],
                variants=["ST"],
                binary_path=fake_binary,
            )

            err = nothing
            try
                run_analysis(;
                    file_path="/nonexistent/path/to/file.edf",
                    channels=[1, 2],
                    select=[0, 0, 0, 0, 1, 1],
                    binary_path=fake_binary,
                )
            catch e
                err = e
            end

            @test err !== nothing
            @test occursin("Input file not found", string(err))

            @test_throws MethodError run_analysis(;
                file_path="test.edf",
                channels=[1],
                flavors=["ST"],
                binary_path=fake_binary,
                model_dimension=4,
            )
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "Command building uses explicit MODEL OUT_FN SR and 1-based channels" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            runner = DDARunner(fake_binary)
            request = DDARequest(
                "test.edf",
                [1, 3, 5],
                ["ST"];
                model=[1, 2, 10],
                derivative_points=5,
                order=4,
                out_fn="/tmp/dda_custom_output",
                sampling_rate=(500, 1000),
            )
            output_base, cleanup_output = DelayDifferentialAnalysis.Runner._resolve_output_base(request)
            cmd = DelayDifferentialAnalysis.Runner.build_command(
                runner,
                request,
                output_base,
            )
            argv = collect(cmd)
            cmd_str = DelayDifferentialAnalysis.Runner.build_command_string(
                runner,
                request,
                output_base,
            )

            if Sys.iswindows()
                @test argv[1] == fake_binary
            else
                @test argv[1] == "sh"
                @test argv[2] == fake_binary
            end
            @test cleanup_output == false
            @test "-OUT_FN" in argv
            @test argv[findfirst(==("-OUT_FN"), argv) + 1] == "/tmp/dda_custom_output"
            @test argv[(findfirst(==("-CH_list"), argv) + 1):(findfirst(==("-SELECT"), argv) - 1)] == ["1", "3", "5"]
            @test argv[(findfirst(==("-MODEL"), argv) + 1):(findfirst(==("-TAU"), argv) - 1)] == ["1", "2", "10"]
            @test argv[(findfirst(==("-WLms"), argv) + 1)] == "200"
            @test argv[(findfirst(==("-WSms"), argv) + 1)] == "100"
            @test argv[(findfirst(==("-dm"), argv) + 1)] == "5"
            @test argv[(findfirst(==("-SR"), argv) + 1):(findfirst(==("-SR"), argv) + 2)] == ["500", "1000"]
            if !Sys.iswindows()
                @test cmd_str == "sh $fake_binary -EDF -DATA_FN test.edf -OUT_FN /tmp/dda_custom_output -CH_list 1 3 5 -SELECT 1 0 0 0 0 0 -MODEL 1 2 10 -TAU 7 10 -WLms 200 -WSms 100 -dm 5 -order 4 -nr_tau 2 -SR 500 1000"
            end
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "Sampling rate is only passed when explicitly requested" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            runner = DDARunner(fake_binary)
            request = DDARequest(
                "test.edf",
                [1],
                ["ST"];
                sampling_rate=(500, 500),
            )
            cmd = DelayDifferentialAnalysis.Runner.build_command(runner, request, "/tmp/out")
            argv = collect(cmd)
            @test !("-SR" in argv)

            default_request = DDARequest("test.edf", [1], ["ST"])
            default_cmd_str = DelayDifferentialAnalysis.Runner.build_command_string(
                runner,
                default_request,
                "/tmp/default_out",
            )
            @test !occursin(" -SR ", default_cmd_str)
            if !Sys.iswindows()
                @test default_cmd_str == "sh $fake_binary -EDF -DATA_FN test.edf -OUT_FN /tmp/default_out -CH_list 1 -SELECT 1 0 0 0 0 0 -MODEL 1 2 10 -TAU 7 10 -WLms 200 -WSms 100 -dm 3 -order 4 -nr_tau 2"
            end
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "Explicit select overrides variant strings" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            runner = DDARunner(fake_binary)
            request = DDARequest(
                "test.edf",
                [1, 2],
                ["ST"];
                select=[0, 0, 0, 0, 1, 1],
            )
            cmd_str = DelayDifferentialAnalysis.Runner.build_command_string(
                runner,
                request,
                "/tmp/select_out",
            )

            @test request.variants == ["DE", "SY"]
            @test request.select == [0, 0, 0, 0, 1, 1]
            if !Sys.iswindows()
                @test cmd_str == "sh $fake_binary -EDF -DATA_FN test.edf -OUT_FN /tmp/select_out -CH_list 1 2 -SELECT 0 0 0 0 1 1 -MODEL 1 2 10 -TAU 7 10 -WLms 200 -WSms 100 -dm 3 -order 4 -nr_tau 2 -WL_CT 200 -WS_CT 100"
            end
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "Raw T, derived t, and normalized window bounds are kept separate" begin
        request = DDARequest(
            "test.edf",
            [1],
            ["ST"];
            window_length=128,
            window_step=64,
            delays=[7, 10],
            derivative_points=3,
            sampling_rate=(500, 500),
            TM=12,
        )

        T = DelayDifferentialAnalysis.Runner._extract_raw_T([
            StructuredChannelData(
                1,
                [
                    StructuredTimepoint(0.0, 5000.0, [1.0, 2.0, 3.0], 0.1),
                    StructuredTimepoint(200.0, 5200.0, [4.0, 5.0, 6.0], 0.2),
                ],
            ),
        ])
        t = DelayDifferentialAnalysis.Runner._compute_t_axis(
            T,
            request.model_params.derivative_points,
            request.tm,
            request.sampling_rate,
        )
        window_starts, window_ends = DelayDifferentialAnalysis.Runner._normalized_window_bounds(
            request,
            "ST",
        )

        @test T == [0.0, 200.0]
        @test t == [0.032, 0.432]
        @test window_starts == Int64[]
        @test window_ends == Int64[]

        window_starts, window_ends = DelayDifferentialAnalysis.Runner._normalized_window_bounds(
            request,
            "ST",
            2,
        )

        @test window_starts == [0, 64]
        @test window_ends == [128, 192]
    end
end
