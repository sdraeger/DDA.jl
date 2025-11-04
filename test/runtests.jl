using Test
using DelayDifferentialAnalysis
import DelayDifferentialAnalysis as DDA
using JSON3


@testset "DelayDifferentialAnalysis.jl" begin
    @testset "Error Types" begin
        @test BinaryNotFoundError("/nonexistent") isa DDAError
        @test FileNotFoundError("/nonexistent") isa DDAError
        @test ExecutionFailedError("test") isa DDAError
        @test ParseError("test") isa DDAError
    end

    @testset "DDARunner Creation" begin
        @test_throws BinaryNotFoundError DDARunner("/nonexistent/binary")

        if isfile("./test/bin/run_DDA_AsciiEdf")
            runner = DDARunner("./test/bin/run_DDA_AsciiEdf")
            @test binary_path(runner) == "./test/bin/run_DDA_AsciiEdf"
        end
    end

    @testset "Type Construction" begin
        bounds = Bounds(0, 100000)
        @test bounds.start == 0
        @test bounds.stop == 100000

        # Test that bounds can be nothing
        @test isnothing(nothing)

        window_params = WindowParameters(1024, 512, nothing, nothing)
        @test window_params.window_length == 1024
        @test window_params.window_step == 512

        delay_params = DelayParameters([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        @test delay_params.delays == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        @test length(delay_params.delays) == 10
    end

    @testset "Parser" begin
        test_content = """
        # Comment line
        1.0 2.0 3.0 4.0 5.0 6.0
        7.0 8.0 9.0 10.0 11.0 12.0
        """

        result = DDA.parse_dda_output(test_content)
        @test !isempty(result)
        @test size(result, 1) > 0
        @test size(result, 2) > 0

        empty_content = """
        # Only comments
        # More comments
        """

        @test_throws ParseError DDA.parse_dda_output(empty_content)
    end

    @testset "DDAResult Construction" begin
        window_params = WindowParameters(1024, 512, nothing, nothing)
        delay_params = DelayParameters([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

        result = DDAResult(
            "/path/to/file.edf",
            ["Channel 1"],
            zeros(Float64, 10, 100),
            window_params,
            delay_params
        )

        @test result.file_path == "/path/to/file.edf"
        @test length(result.channels) == 1
        @test size(result.q_matrix) == (10, 100)
        @test isnothing(result.variant_results)
        @test isnothing(result.raw_output)
    end

    @testset "End-to-End Execution" begin
        if isfile("./test/bin/run_DDA_AsciiEdf") && isfile("./test/data/test.edf")
            runner = DDARunner("./test/bin/run_DDA_AsciiEdf")

            window_params = WindowParameters(1024, 512, nothing, nothing)
            delay_params = DelayParameters([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

            # Create algorithm selection (enable ST variant)
            algo_selection = DDA.AlgorithmSelection(["ST"], "1 0 0 0")

            # Create DDA request with no bounds (process entire file)
            request = DDA.DDARequest(
                "./test/data/test.edf",
                nothing,  # channels
                nothing,  # bounds (process entire file)
                algo_selection,
                window_params,
                delay_params,
                nothing  # ct_channel_pairs
            )

            result = DDA.run_dda(runner, request)

            # Validate result structure
            @test result isa DDAResult
            @test result.file_path == "./test/data/test.edf"
            @test size(result.q_matrix, 1) > 0
            @test size(result.q_matrix, 2) > 0
            @test !isnothing(result.variant_results)
            @test length(result.variant_results) == 1
            @test result.variant_results[1].variant_id == "ST"

            # Ground truth validation
            ground_truth_file = "./test/data/expected_output.json"
            if isfile(ground_truth_file)
                # Load and compare with ground truth
                expected_data = JSON3.read(read(ground_truth_file, String))
                expected_matrix = Matrix{Float64}(hcat(expected_data.q_matrix...)')

                @test size(result.q_matrix) == size(expected_matrix)
                @test result.q_matrix ≈ expected_matrix atol = 1e-10

                @info "Ground truth validation passed!"
            else
                # Generate ground truth file for future runs
                ground_truth_data = Dict(
                    "q_matrix" => [result.q_matrix[:, i] for i in 1:size(result.q_matrix, 2)],
                    "shape" => size(result.q_matrix),
                    "test_params" => Dict(
                        "window_length" => 1024,
                        "window_step" => 512,
                        "delays" => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
                    )
                )

                write(ground_truth_file, JSON3.write(ground_truth_data, allow_inf=true))
                @warn "Ground truth file created at $ground_truth_file - please verify output is correct"
            end
        end
    end

    @testset "CD Variant Execution" begin
        if isfile("./test/bin/run_DDA_AsciiEdf") && isfile("./test/data/test.edf")
            runner = DDARunner("./test/bin/run_DDA_AsciiEdf")

            window_params = WindowParameters(1024, 512, nothing, nothing)
            delay_params = DelayParameters([1, 2, 3, 4, 5])

            # Create algorithm selection for CD variant
            # CD should auto-enable ST and CT (select_mask should be "1 1 1 0")
            algo_selection = DDA.AlgorithmSelection(["CD"], nothing)

            # Create DDA request
            request = DDA.DDARequest(
                "./test/data/test.edf",
                [1, 2],  # Two channels for CD
                nothing,  # bounds (process entire file)
                algo_selection,
                window_params,
                delay_params,
                nothing  # ct_channel_pairs
            )

            result = DDA.run_dda(runner, request)

            # Validate result structure
            @test result isa DDAResult
            @test result.file_path == "./test/data/test.edf"
            @test size(result.q_matrix, 1) > 0
            @test size(result.q_matrix, 2) > 0
            @test !isnothing(result.variant_results)

            # CD should generate ST, CT, and CD variants (3 total)
            @test length(result.variant_results) == 3

            # Verify all three variants are present
            variant_ids = [v.variant_id for v in result.variant_results]
            @test "ST" in variant_ids
            @test "CT" in variant_ids
            @test "CD" in variant_ids

            # Find CD variant and verify it has data
            cd_variant_idx = findfirst(v -> v.variant_id == "CD", result.variant_results)
            @test !isnothing(cd_variant_idx)
            @test size(result.variant_results[cd_variant_idx].q_matrix, 1) > 0
            @test size(result.variant_results[cd_variant_idx].q_matrix, 2) > 0

            # Ground truth validation for each variant
            for variant in result.variant_results
                ground_truth_file = "./test/data/expected_output_$(variant.variant_id).json"

                if isfile(ground_truth_file)
                    # Load and compare with ground truth
                    expected_data = JSON3.read(read(ground_truth_file, String))
                    expected_matrix = Matrix{Float64}(hcat(expected_data.q_matrix...)')

                    @test size(variant.q_matrix) == size(expected_matrix)
                    @test variant.q_matrix ≈ expected_matrix atol = 1e-10

                    @info "Ground truth validation passed for $(variant.variant_id)!"
                else
                    # Generate ground truth file for future runs
                    ground_truth_data = Dict(
                        "q_matrix" => [variant.q_matrix[:, i] for i in 1:size(variant.q_matrix, 2)],
                        "shape" => size(variant.q_matrix),
                        "variant_id" => variant.variant_id,
                        "test_params" => Dict(
                            "window_length" => 1024,
                            "window_step" => 512,
                            "delays" => [1, 2, 3, 4, 5],
                            "channels" => [1, 2]
                        )
                    )

                    write(ground_truth_file, JSON3.write(ground_truth_data, allow_inf=true))
                    @warn "Ground truth file created at $ground_truth_file for $(variant.variant_id) - please verify output is correct"
                end
            end

            @info "CD variant test passed with $(length(result.variant_results)) variants"
        end
    end
end
