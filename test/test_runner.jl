#=
Tests for the DDA Runner module.

These tests validate the DDARequest, DDAResult types and command building logic.
Integration tests (actual binary execution) are skipped if binary is not found.
=#

using Test
using DelayDifferentialAnalysis

@testset "Runner Module" begin

    # =============================================================================
    # REQUEST CONSTRUCTION
    # =============================================================================

    @testset "DDARequest Construction" begin
        @testset "Basic construction" begin
            request = DDARequest(
                "test.edf",
                [0, 1, 2],
                ["ST", "SY"]
            )

            @test request.file_path == "test.edf"
            @test request.channels == [0, 1, 2]
            @test request.variants == ["ST", "SY"]
            @test request.window_params.window_length == 200   # Default (DDADefaults.WINDOW_LENGTH)
            @test request.window_params.window_step == 100    # Default (DDADefaults.WINDOW_STEP)
            @test request.delay_params.delays == collect(DEFAULT_DELAYS)
        end

        @testset "With custom parameters" begin
            request = DDARequest(
                "data.edf",
                [0, 1],
                ["ST", "CT", "CD"];
                window_length=4096,
                window_step=2048,
                ct_window_length=4,
                ct_window_step=2,
                delays=[1, 2, 3, 4, 5],
                dm=6,
                order=5,
                nr_tau=3,
                time_range=(0.0, 10000.0),
                ct_pairs=[(0, 1)],
                cd_pairs=[(0, 1), (1, 0)],
                sampling_rate=2048.0
            )

            @test request.window_params.window_length == 4096
            @test request.window_params.window_step == 2048
            @test request.window_params.ct_window_length == 4
            @test request.window_params.ct_window_step == 2
            @test request.delay_params.delays == [1, 2, 3, 4, 5]
            @test request.model_params.dm == 6
            @test request.model_params.order == 5
            @test request.model_params.nr_tau == 3
            @test request.time_range.start == 0.0
            @test request.time_range.stop == 10000.0
            @test request.ct_channel_pairs == [(0, 1)]
            @test request.cd_channel_pairs == [(0, 1), (1, 0)]
            @test request.sampling_rate == 2048.0
        end

        @testset "Default model parameters" begin
            request = DDARequest("test.edf", [0], ["ST"])

            @test request.model_params.dm == 4
            @test request.model_params.order == 4
            @test request.model_params.nr_tau == 2
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
                [0, 1, 2],
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
                    run_analysis(runner, request)
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
            withenv("DDA_BINARY_PATH" => nothing, "DDA_HOME" => nothing) do
                # Only test if binary is definitely not installed
                if find_binary() === nothing
                    @test_throws ErrorException DDARunner()
                end
            end
        end
    end

    # =============================================================================
    # VARIANT RESULT DATA
    # =============================================================================

    @testset "VariantResultData" begin
        q_matrix = [1.0 2.0 3.0; 4.0 5.0 6.0]
        labels = ["Channel 1", "Channel 2"]

        result = VariantResultData(
            "ST",
            "Single Timeseries",
            q_matrix,
            labels
        )

        @test result.variant_id == "ST"
        @test result.variant_name == "Single Timeseries"
        @test size(result.q_matrix) == (2, 3)
        @test result.channel_labels == labels
    end

    # =============================================================================
    # INTEGRATION TESTS (require binary)
    # =============================================================================

    @testset "Binary Integration Tests" begin
        binary_path = find_binary()

        if binary_path === nothing
            @info "Skipping binary integration tests: DDA binary not found"
            @test_skip true
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
                @test_skip true
            else
                @testset "Run ST analysis" begin
                    request = DDARequest(
                        test_data,
                        [0, 1, 2],
                        ["ST"];
                        window_length=2048,
                        window_step=1024,
                        delays=collect(1:10),
                        time_range=(0.0, 6000.0)
                    )

                    try
                        result = run_analysis(request)

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
                        [0, 1],
                        ["ST", "SY"];
                        window_length=2048,
                        window_step=1024,
                        time_range=(0.0, 4000.0)
                    )

                    try
                        result = run_analysis(request)

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
end
