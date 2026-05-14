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
            @test request.window_params.WL === nothing
            @test request.window_params.WS === nothing
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
                WL=4096,
                WS=2048,
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

            @test request.window_params.WL == 4096
            @test request.window_params.WS == 2048
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

        @testset "Scalar sampling rate is preserved" begin
            request = DDARequest("test.edf", [1], ["ST"]; sampling_rate=2048)
            @test request.sampling_rate == 2048
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

        @testset "run_DDA fails gracefully with missing file" begin
            request = DDARequest(
                "/nonexistent/path/to/file.edf",
                [1, 2, 3],
                ["ST", "SY"];
                WL=2048,
                WS=1024
            )

            # Create a fake binary path to bypass binary discovery
            fake_binary = tempname()
            touch(fake_binary)

            try
                runner = DDARunner(fake_binary)
                # This should fail with "Input file not found", NOT "UndefVarError: UUIDs"
                err = nothing
                try
                    run_DDA(; runner=runner, request=request)
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
        T = [
            0.0 128.0
            100.0 228.0
            200.0 328.0
        ]
        A = [
            1.0 2.0 0.1 7.0 8.0 0.4
            3.0 4.0 0.2 9.0 10.0 0.5
            5.0 6.0 0.3 11.0 12.0 0.6
        ]
        coefficients = Array{Float64,3}(undef, 2, 3, 2)
        coefficients[1, 1, :] = [1.0, 2.0]
        coefficients[1, 2, :] = [3.0, 4.0]
        coefficients[1, 3, :] = [5.0, 6.0]
        coefficients[2, 1, :] = [7.0, 8.0]
        coefficients[2, 2, :] = [9.0, 10.0]
        coefficients[2, 3, :] = [11.0, 12.0]
        errors = [0.1 0.2 0.3; 0.4 0.5 0.6]
        t = [0.014, 0.114, 0.214]
        window_starts = [0, 64, 128]
        window_ends = [128, 192, 256]
        labels = ["Channel 1", "Channel 2"]

        result = VariantResultData(
            "ST",
            "Single Timeseries",
            A,
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
        @test result.A == A
        @test !(:q_matrix in fieldnames(typeof(result)))
        @test eltype(result.T) == Int64
        @test size(result.coefficients) == (2, 3, 2)
        @test size(result.errors) == (2, 3)
        @test result.T == T
        @test result.t == t
        @test result.window_starts == window_starts
        @test result.window_ends == window_ends
        @test result.channel_labels == labels
    end

    @testset "DDAResult exposes primary derived t axis" begin
        request = DDARequest(
            "test.edf",
            [1, 2],
            ["ST"];
            derivative_points=3,
            TM=12,
            sampling_rate=(500, 500),
        )
        T = Int64[
            0 128
            200 328
        ]
        t = (T[:, 1] .+ 1 .+ request.tm .+ request.model_params.derivative_points) ./ 500
        A = [
            1.0 2.0 0.1 7.0 8.0 0.4
            3.0 4.0 0.2 9.0 10.0 0.5
        ]
        coefficients = Array{Float64,3}(undef, 2, 2, 2)
        coefficients[1, 1, :] = [1.0, 2.0]
        coefficients[1, 2, :] = [3.0, 4.0]
        coefficients[2, 1, :] = [7.0, 8.0]
        coefficients[2, 2, :] = [9.0, 10.0]
        errors = [0.1 0.2; 0.4 0.5]
        labels = ["Channel 1", "Channel 2"]
        variant = VariantResultData(
            "ST",
            "Single Timeseries",
            A,
            coefficients,
            errors,
            T,
            t,
            [0, 200],
            [128, 328],
            labels,
        )

        result = DDAResult(
            "analysis",
            "test.edf",
            labels,
            T,
            t,
            A,
            [variant],
            request.window_params,
            request.delay_params,
            "2026-05-14T00:00:00",
        )
        legacy_result = DDAResult(
            "analysis",
            "test.edf",
            labels,
            T,
            A,
            [variant],
            request.window_params,
            request.delay_params,
            "2026-05-14T00:00:00",
        )

        @test :t in fieldnames(typeof(result))
        @test result.T == T
        @test result.t == [0.032, 0.432]
        @test result.t == t
        @test legacy_result.t == t
        @test result.ST === variant
        @test result.ST.A == A
        @test :ST in propertynames(result)
        @test !(:CT in propertynames(result))
        @test_throws ErrorException result.CT
    end

    @testset "VariantResultData partitions raw binary output into T and A" begin
        channels = [
            StructuredChannelData(
                1,
                [
                    StructuredTimepoint(10.0, 138.0, [1.0, 2.0, 3.0], 0.1),
                    StructuredTimepoint(110.0, 238.0, [4.0, 5.0, 6.0], 0.2),
                ],
            ),
            StructuredChannelData(
                2,
                [
                    StructuredTimepoint(10.0, 138.0, [7.0, 8.0, 9.0], 0.3),
                    StructuredTimepoint(110.0, 238.0, [10.0, 11.0, 12.0], 0.4),
                ],
            ),
        ]
        request = DDARequest("test.edf", [1, 2], ["ST"]; sampling_rate=nothing)
        variant = get_variant_by_abbrev("ST")

        result = DelayDifferentialAnalysis.Runner._pack_variant_result(
            "ST",
            variant,
            channels,
            request,
            ["Channel 1", "Channel 2"],
        )

        @test result.T == [
            10.0 138.0
            110.0 238.0
        ]
        @test eltype(result.T) == Int64
        @test result.A == [
            1.0 2.0 3.0 0.1 7.0 8.0 9.0 0.3
            4.0 5.0 6.0 0.2 10.0 11.0 12.0 0.4
        ]
        @test result.coefficients[1, 2, :] == [4.0, 5.0, 6.0]
        @test result.coefficients[2, 1, :] == [7.0, 8.0, 9.0]
        @test result.errors == [0.1 0.2; 0.3 0.4]
        @test result.t == [24.0, 124.0]
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
                        WL=2048,
                        WS=1024,
                        delays=collect(1:10),
                        time_range=(0.0, 6000.0)
                    )

                    try
                        result = run_DDA(; request=request)

                        @test !isempty(result.id)
                        @test result.file_path == test_data
                        @test size(result.T, 1) > 0
                        @test size(result.T, 2) == 2
                        @test eltype(result.T) == Int64
                        @test length(result.t) == size(result.T, 1)
                        @test size(result.A, 1) > 0
                        @test size(result.A, 2) > 1
                        @test length(result.variant_results) == 1
                        @test result.variant_results[1].variant_id == "ST"
                        @test result.t == result.variant_results[1].t

                        @info "ST result" size=size(result.A) windows=size(result.A, 1)
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
                        WL=2048,
                        WS=1024,
                        time_range=(0.0, 4000.0)
                    )

                    try
                        result = run_DDA(; request=request)

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

    @testset "Direct run_DDA API" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            @test :run_DDA in names(DelayDifferentialAnalysis)
            @test !(:run_analysis in names(DelayDifferentialAnalysis))

            @test_throws UndefVarError run_analysis

            @test_throws MethodError run_DDA(
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
                run_DDA(;
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

            @test_throws MethodError run_DDA(;
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
                run_DDA(;
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

            @test_throws MethodError run_DDA(;
                file_path="test.edf",
                channels=[1],
                flavors=["ST"],
                binary_path=fake_binary,
                model_dimension=4,
            )

            @test_throws MethodError DDARequest(
                "test.edf",
                [1],
                ["ST"];
                window_length=200,
                window_step=100,
            )
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "Command building uses explicit MODEL OUT_FN SR WL WS and 1-based channels" begin
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
                WL=200,
                WS=100,
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
            @test !("-WLms" in argv)
            @test !("-WSms" in argv)
            @test argv[(findfirst(==("-WL"), argv) + 1)] == "200"
            @test argv[(findfirst(==("-WS"), argv) + 1)] == "100"
            @test argv[(findfirst(==("-dm"), argv) + 1)] == "5"
            @test argv[(findfirst(==("-SR"), argv) + 1):(findfirst(==("-SR"), argv) + 2)] == ["500", "1000"]
            if !Sys.iswindows()
                @test cmd_str == "sh $fake_binary -EDF -DATA_FN test.edf -OUT_FN /tmp/dda_custom_output -CH_list 1 3 5 -SELECT 1 0 0 0 0 0 -MODEL 1 2 10 -TAU 7 10 -WL 200 -WS 100 -dm 5 -order 4 -nr_tau 2 -SR 500 1000"
            end
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "Sampling rate is emitted as scalar, tuple, or omitted by default" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            runner = DDARunner(fake_binary)
            scalar_request = DDARequest(
                "test.edf",
                [1],
                ["ST"];
                sampling_rate=500,
            )
            scalar_argv = collect(DelayDifferentialAnalysis.Runner.build_command(runner, scalar_request, "/tmp/scalar_out"))
            scalar_sr_index = findfirst(==("-SR"), scalar_argv)
            @test scalar_sr_index !== nothing
            @test scalar_argv[scalar_sr_index + 1] == "500"
            @test scalar_sr_index + 1 == length(scalar_argv) ||
                  scalar_argv[scalar_sr_index + 2] != "500"

            tuple_request = DDARequest(
                "test.edf",
                [1],
                ["ST"];
                sampling_rate=(500, 500),
            )
            tuple_argv = collect(DelayDifferentialAnalysis.Runner.build_command(runner, tuple_request, "/tmp/tuple_out"))
            tuple_sr_index = findfirst(==("-SR"), tuple_argv)
            @test tuple_sr_index !== nothing
            @test tuple_argv[(tuple_sr_index + 1):(tuple_sr_index + 2)] == ["500", "500"]

            @test_throws ErrorException DDARequest("test.edf", [1], ["ST"]; sampling_rate=500.5)
            @test_throws ErrorException DDARequest("test.edf", [1], ["ST"]; sampling_rate=(500, 1000.5))

            default_request = DDARequest("test.edf", [1], ["ST"])
            default_cmd_str = DelayDifferentialAnalysis.Runner.build_command_string(
                runner,
                default_request,
                "/tmp/default_out",
            )
            @test !occursin(" -SR ", default_cmd_str)
            @test !occursin(" -WL ", default_cmd_str)
            @test !occursin(" -WS ", default_cmd_str)
            @test !occursin(" -WLms ", default_cmd_str)
            @test !occursin(" -WSms ", default_cmd_str)
            if !Sys.iswindows()
                @test default_cmd_str == "sh $fake_binary -EDF -DATA_FN test.edf -OUT_FN /tmp/default_out -CH_list 1 -SELECT 1 0 0 0 0 0 -MODEL 1 2 10 -TAU 7 10 -dm 3 -order 4 -nr_tau 2"
            end
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "run_DDA passthrough kwargs map directly to binary flags" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            runner = DDARunner(fake_binary)
            request = DDARequest(
                "test.edf",
                [1, 2],
                ["ST"];
                tau_file="/tmp/tau_values.txt",
                tau2=[11, 12],
                model2=[2, 5, 9],
                WL_CT=33,
                WS_CT=22,
                no_norm=true,
                WN_list=[4, 8, 12],
            )
            cmd_str = DelayDifferentialAnalysis.Runner.build_command_string(
                runner,
                request,
                "/tmp/passthrough_out",
            )
            argv = collect(DelayDifferentialAnalysis.Runner.build_command(
                runner,
                request,
                "/tmp/passthrough_out",
            ))

            @test "-TAU_file" in argv
            @test argv[findfirst(==("-TAU_file"), argv) + 1] == "/tmp/tau_values.txt"
            @test argv[(findfirst(==("-TAU2"), argv) + 1):(findfirst(==("-MODEL2"), argv) - 1)] == ["11", "12"]
            @test argv[(findfirst(==("-MODEL2"), argv) + 1):(findfirst(==("-WL_CT"), argv) - 1)] == ["2", "5", "9"]
            @test argv[findfirst(==("-WL_CT"), argv) + 1] == "33"
            @test "-WS_CT" in argv
            @test !("-WS_ct" in argv)
            @test argv[findfirst(==("-WS_CT"), argv) + 1] == "22"
            @test "-NoNorm" in argv
            @test "-WN_list" in argv
            @test !("WN_list" in argv)
            @test argv[(findfirst(==("-WN_list"), argv) + 1):end] == ["4", "8", "12"]
            if !Sys.iswindows()
                @test cmd_str == "sh $fake_binary -EDF -DATA_FN test.edf -OUT_FN /tmp/passthrough_out -CH_list 1 2 -SELECT 1 0 0 0 0 0 -MODEL 1 2 10 -TAU 7 10 -dm 3 -order 4 -nr_tau 2 -TAU_file /tmp/tau_values.txt -TAU2 11 12 -MODEL2 2 5 9 -WL_CT 33 -WS_CT 22 -NoNorm -WN_list 4 8 12"
            end
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "run_DDA passthrough defaults are absent and list types are validated" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            runner = DDARunner(fake_binary)
            request = DDARequest("test.edf", [1], ["ST"])
            argv = collect(DelayDifferentialAnalysis.Runner.build_command(
                runner,
                request,
                "/tmp/default_passthrough_out",
            ))

            @test !("-TAU_file" in argv)
            @test !("-TAU2" in argv)
            @test !("-MODEL2" in argv)
            @test !("-WL_CT" in argv)
            @test !("-WS_CT" in argv)
            @test !("-NoNorm" in argv)
            @test !("-WN_list" in argv)

            @test_throws ErrorException DDARequest("test.edf", [1], ["ST"]; tau2=7)
            @test_throws ErrorException DDARequest("test.edf", [1], ["ST"]; model2=[1, "bad"])
            @test_throws ErrorException DDARequest("test.edf", [1], ["ST"]; WN_list=3)
            @test_throws MethodError DDARequest("test.edf", [1], ["ST"]; WL_ct=33)
            @test_throws MethodError DDARequest("test.edf", [1], ["ST"]; WS_ct=22)
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "CT window flags are explicit and not inferred from WL WS" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            runner = DDARunner(fake_binary)
            inferred_request = DDARequest(
                "test.edf",
                [1, 2],
                ["DE"];
                WL=200,
                WS=100,
            )
            inferred_argv = collect(DelayDifferentialAnalysis.Runner.build_command(
                runner,
                inferred_request,
                "/tmp/no_inferred_ct_out",
            ))

            @test "-WL" in inferred_argv
            @test "-WS" in inferred_argv
            @test !("-WL_CT" in inferred_argv)
            @test !("-WS_CT" in inferred_argv)

            explicit_request = DDARequest(
                "test.edf",
                [1, 2],
                ["DE"];
                WL=200,
                WS=100,
                ct_window_length=33,
                ct_window_step=22,
            )
            explicit_argv = collect(DelayDifferentialAnalysis.Runner.build_command(
                runner,
                explicit_request,
                "/tmp/explicit_ct_out",
            ))

            @test explicit_argv[findfirst(==("-WL_CT"), explicit_argv) + 1] == "33"
            @test explicit_argv[findfirst(==("-WS_CT"), explicit_argv) + 1] == "22"
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
                WL=200,
                WS=100,
            )
            cmd_str = DelayDifferentialAnalysis.Runner.build_command_string(
                runner,
                request,
                "/tmp/select_out",
            )

            @test request.variants == ["DE", "SY"]
            @test request.select == [0, 0, 0, 0, 1, 1]
            if !Sys.iswindows()
                @test cmd_str == "sh $fake_binary -EDF -DATA_FN test.edf -OUT_FN /tmp/select_out -CH_list 1 2 -SELECT 0 0 0 0 1 1 -MODEL 1 2 10 -TAU 7 10 -WL 200 -WS 100 -dm 3 -order 4 -nr_tau 2"
            end
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "run_DDA executes CT once per pair when mixed with ST" begin
        if Sys.iswindows()
            @info "Skipping CT pair execution test on Windows"
        else
            temp_dir = mktempdir()
            fake_binary = joinpath(temp_dir, "fake_dda.sh")
            input_file = joinpath(temp_dir, "input.txt")
            output_base = joinpath(temp_dir, "mixed")

            open(input_file, "w") do io
                println(io, "1 2 3")
                println(io, "4 5 6")
            end

            open(fake_binary, "w") do io
                write(io, raw"""#!/usr/bin/env sh
out=''
st=0
ct=0
channels=''
while [ "$#" -gt 0 ]; do
  case "$1" in
    -OUT_FN) out="$2"; shift 2 ;;
    -SELECT) st="$2"; ct="$3"; shift 7 ;;
    -CH_list)
      shift
      channels=''
      while [ "$#" -gt 0 ] && [ "${1#-}" = "$1" ]; do
        channels="$channels $1"
        shift
      done
      ;;
    *) shift ;;
  esac
done
: > "$out.info"
if [ "$st" = "1" ]; then
  line='0 10'
  for _ch in $channels; do line="$line 1.0 2.0 3.0 0.1"; done
  printf '%s\n' "$line" > "${out}_ST"
fi
if [ "$ct" = "1" ]; then
  printf '%s\n' '0 10 4.0 5.0 6.0 0.2' > "${out}_CT"
fi
""")
            end
            chmod(fake_binary, 0o755)

            try
                @test_throws ErrorException run_DDA(;
                    file_path=input_file,
                    channels=[1, 2, 3],
                    flavors=["ST", "CT"],
                    binary_path=fake_binary,
                    WL=128,
                    WS=100,
                    WL_CT=128,
                    WS_CT=100,
                    out_fn=joinpath(temp_dir, "invalid_ct"),
                )

                result = run_DDA(;
                    file_path=input_file,
                    channels=[1, 2, 3],
                    flavors=["ST", "CT"],
                    binary_path=fake_binary,
                    WL=128,
                    WS=100,
                    out_fn=output_base,
                )

                variant_ids = [vr.variant_id for vr in result.variant_results]
                @test variant_ids == ["ST", "CT"]
                @test result.ST === result.variant_results[1]
                @test result.CT === result.variant_results[2]
                @test result.ST.A == result.variant_results[1].A
                @test result.CT.A == result.variant_results[2].A
                @test :ST in propertynames(result)
                @test :CT in propertynames(result)
                @test !(:DE in propertynames(result))
                ct_result = result.variant_results[2]
                @test size(ct_result.A) == (1, 12)
                @test ct_result.channel_labels == [
                    "Channel 1-Channel 2",
                    "Channel 1-Channel 3",
                    "Channel 2-Channel 3",
                ]
            finally
                rm(temp_dir; recursive=true, force=true)
            end
        end
    end

    @testset "Raw T, derived t, and normalized window bounds are kept separate" begin
        request = DDARequest(
            "test.edf",
            [1],
            ["ST"];
            WL=128,
            WS=64,
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

        ct_request = DDARequest(
            "test.edf",
            [1, 2],
            ["CT"];
            WL=128,
            WS=64,
            ct_window_length=2,
            ct_window_step=2,
        )
        ct_result = DelayDifferentialAnalysis.Runner._pack_variant_result(
            "CT",
            DelayDifferentialAnalysis.CT,
            [
                StructuredChannelData(
                    1,
                    [
                        StructuredTimepoint(0.0, 144.0, [1.0, 2.0, 3.0], 0.1),
                        StructuredTimepoint(100.0, 244.0, [4.0, 5.0, 6.0], 0.2),
                    ],
                ),
            ],
            ct_request,
            ["Channel 1-Channel 2"],
        )
        @test ct_result.window_starts == [0, 100]
        @test ct_result.window_ends == [144, 244]
    end
end
