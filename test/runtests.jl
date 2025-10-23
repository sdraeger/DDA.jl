using Test
using DDAJl


@testset "DDAJl.jl" begin
    @testset "Error Types" begin
        @test BinaryNotFoundError("/nonexistent") isa DDAError
        @test FileNotFoundError("/nonexistent") isa DDAError
        @test ExecutionFailedError("test") isa DDAError
        @test ParseError("test") isa DDAError
    end

    @testset "DDARunner Creation" begin
        @test_throws BinaryNotFoundError DDARunner("/nonexistent/binary")

        if isfile("./bin/run_DDA_AsciiEdf")
            runner = DDARunner("./bin/run_DDA_AsciiEdf")
            @test binary_path(runner) == "./bin/run_DDA_AsciiEdf"
        end
    end

    @testset "Type Construction" begin
        time_range = TimeRange(0.0, 100.0)
        @test time_range.start == 0.0
        @test time_range.stop == 100.0

        window_params = WindowParameters(1024, 512, nothing, nothing)
        @test window_params.window_length == 1024
        @test window_params.window_step == 512

        scale_params = ScaleParameters(1.0, 10.0, 10)
        @test scale_params.scale_min == 1.0
        @test scale_params.scale_max == 10.0
        @test scale_params.scale_num == 10
    end

    @testset "Parser" begin
        test_content = """
        # Comment line
        1.0 2.0 3.0 4.0 5.0 6.0
        7.0 8.0 9.0 10.0 11.0 12.0
        """

        result = DDAJl.parse_dda_output(test_content)
        @test !isempty(result)
        @test size(result, 1) > 0
        @test size(result, 2) > 0

        empty_content = """
        # Only comments
        # More comments
        """

        @test_throws ParseError DDAJl.parse_dda_output(empty_content)
    end

    @testset "DDAResult Construction" begin
        window_params = WindowParameters(1024, 512, nothing, nothing)
        scale_params = ScaleParameters(1.0, 10.0, 10)

        result = DDAResult(
            "test-id",
            "/path/to/file.edf",
            ["Channel 1"],
            zeros(Float64, 10, 100),
            window_params,
            scale_params
        )

        @test result.id == "test-id"
        @test result.file_path == "/path/to/file.edf"
        @test length(result.channels) == 1
        @test size(result.q_matrix) == (10, 100)
        @test isnothing(result.variant_results)
        @test isnothing(result.raw_output)
    end

    @testset "End-to-End Execution" begin
        if isfile("./bin/run_DDA_AsciiEdf") && isfile("./data/test.edf")
            runner = DDARunner("./bin/run_DDA_AsciiEdf")

            time_range = TimeRange(0.0, 100.0)
            window_params = WindowParameters(1024, 512, nothing, nothing)
            scale_params = ScaleParameters(1.0, 10.0, 10)

            # Create preprocessing options
            preprocessing_opts = DDAJl.PreprocessingOptions(nothing, nothing, nothing)

            # Create algorithm selection (enable ST variant)
            algo_selection = DDAJl.AlgorithmSelection(["ST"], "1 0 0 0")

            # Create DDA request
            request = DDAJl.DDARequest(
                "./data/test.edf",
                nothing,  # channels
                time_range,
                preprocessing_opts,
                algo_selection,
                window_params,
                scale_params,
                nothing  # ct_channel_pairs
            )

            result = DDAJl.run(runner, request)

            println("DDA Result: ", result)

            @test result isa DDAResult
            @test result.file_path == "./data/test.edf"
            @test size(result.q_matrix, 1) > 0
            @test size(result.q_matrix, 2) > 0
            @test !isnothing(result.variant_results)
            @test length(result.variant_results) == 1
            @test result.variant_results[1].variant_id == "ST"
        end
    end
end
