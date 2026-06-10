using Test
using DelayDifferentialAnalysis

@testset "Structure selection" begin
    function fake_result(error)
        return (
            variant_results=[
                (
                    variant_id="ST",
                    errors=reshape([Float64(error)], 1, 1),
                ),
            ],
        )
    end

    @testset "selects the lowest-error model and delay candidate" begin
        calls = []
        scores = Dict(
            ([1, 2, 6], [7, 10]) => 0.5,
            ([1, 2, 6], [10, 20]) => 0.3,
            ([1, 2, 10], [7, 10]) => 0.7,
            ([1, 2, 10], [10, 20]) => 0.1,
        )

        run_once = (; model, delays, kwargs...) -> begin
            model_vector = Int[model...]
            delay_vector = Int[delays...]
            push!(calls, (model=model_vector, delays=delay_vector, kwargs=kwargs))
            return fake_result(scores[(model_vector, delay_vector)])
        end

        result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
            run_once;
            file_path="data.ascii",
            channels=[1, 2],
            binary_path="/tmp/run_DDA_AsciiEdf",
            candidate_models=[[1, 2, 6], [1, 2, 10]],
            candidate_delays=[[7, 10], [10, 20]],
            derivative_points=4,
            order=3,
            WL=3000,
            WS=200,
            input_format=:ascii,
        )

        @test result.best_model == [1, 2, 10]
        @test result.best_delays == [10, 20]
        @test result.best_score == 0.1
        @test length(result.trials) == 4
        @test calls[1].model == [1, 2, 6]
        @test calls[1].delays == [7, 10]
        @test calls[1].kwargs[:flavors] == ["ST"]
        @test calls[1].kwargs[:WL] == 3000
        @test calls[1].kwargs[:WS] == 200
        @test calls[1].kwargs[:input_format] == :ascii
        @test calls[1].kwargs[:nr_tau] == 2
    end

    @testset "supports matrix model candidates and stable output paths" begin
        out_dir = mktempdir()
        matrix_model = [0 0 1; 0 0 2; 1 1 1]
        seen = []

        run_once = (; model, out_fn, kwargs...) -> begin
            push!(seen, (model=model, out_fn=out_fn))
            return fake_result(length(seen))
        end

        try
            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                binary_path="/tmp/run_DDA_AsciiEdf",
                candidate_models=[matrix_model],
                candidate_delays=[[32, 9]],
                derivative_points=4,
                order=3,
                out_dir=out_dir,
            )

            @test result.best_model == matrix_model
            @test seen[1].model == matrix_model
            @test startswith(seen[1].out_fn, out_dir)
            @test basename(seen[1].out_fn) == "structure_selection_m1_d1"
            @test result.trials[1].out_fn == seen[1].out_fn
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "scores supported error metrics" begin
        result = (
            variant_results=[
                (variant_id="ST", errors=[1.0 2.0 9.0]),
            ],
        )

        @test DelayDifferentialAnalysis.StructureSelection._score_result(result, :mean_error) == 4.0
        @test DelayDifferentialAnalysis.StructureSelection._score_result(result, :median_error) == 2.0
        @test DelayDifferentialAnalysis.StructureSelection._score_result(result, :minimum_error) == 1.0
        @test_throws ErrorException DelayDifferentialAnalysis.StructureSelection._score_result(
            (variant_results=[(variant_id="CT", errors=[1.0])],),
            :mean_error,
        )
        @test_throws ErrorException DelayDifferentialAnalysis.StructureSelection._score_result(result, :unknown)
    end
end
