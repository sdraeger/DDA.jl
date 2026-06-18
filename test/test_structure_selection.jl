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

    function fake_errors(errors)
        return (
            variant_results=[
                (
                    variant_id="ST",
                    errors=Float64.(errors),
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

    @testset "delays keyword accepts explicit nested delay candidates" begin
        calls = []
        scores = Dict(
            ([1, 2, 6], [7, 10]) => 0.5,
            ([1, 2, 6], [10, 20]) => 0.2,
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
            candidate_models=[[1, 2, 6]],
            delays=[[7, 10], [10, 20]],
            derivative_points=4,
            order=3,
        )

        @test result.best_delays == [10, 20]
        @test length(calls) == 2
        @test calls[1].delays == [7, 10]
        @test !haskey(calls[1].kwargs, :tau_file)
    end

    @testset "flat delay pools generate model-specific tau files" begin
        calls = []

        run_once = (; model, tau_file, delays, nr_tau, out_fn, kwargs...) -> begin
            model_vector = Int[model...]
            push!(calls, (
                model=model_vector,
                tau_file=tau_file,
                delays=Int[delays...],
                nr_tau=nr_tau,
                out_fn=out_fn,
                kwargs=kwargs,
            ))
            if model_vector == [1]
                return fake_errors([0.7, 0.4, 0.6])
            elseif model_vector == [4]
                return fake_errors([0.1, 0.8, 0.5])
            end
            return fake_errors([0.9])
        end

        out_dir = mktempdir()
        try
            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                binary_path="/tmp/run_DDA_AsciiEdf",
                candidate_models=[[1], [4]],
                delays=10:20:50,
                derivative_points=4,
                order=2,
                nr_delays=2,
                out_dir=out_dir,
                tau_file_suffix="_run1",
            )

            @test result.best_model == [4]
            @test result.best_delays == [10, 30]
            @test result.best_score == 0.1
            @test result.artifacts_dir !== nothing
            @test dirname(calls[1].tau_file) == result.artifacts_dir
            @test basename(result.artifacts_dir) |> startswith("structure_selection_")
            @test sort(basename.(calls[i].tau_file for i in eachindex(calls))) == [
                "TAU_ALL__1_0_run1",
                "TAU_ALL__2_1_run1",
            ]
            @test read(calls[1].tau_file, String) == "10\n30\n50\n"
            @test read(calls[2].tau_file, String) == "10 30\n10 50\n30 50\n"
            @test calls[1].delays == [10]
            @test calls[1].nr_tau == 1
            @test calls[2].delays == [10, 30]
            @test calls[2].nr_tau == 2
            @test calls[1].kwargs[:flavors] == ["ST"]
            @test result.trials[1].tau_file == calls[1].tau_file
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "generated tau-file selection maps error rows back to delays" begin
        run_once = (; kwargs...) -> fake_errors([0.5; 0.2; 0.9])

        out_dir = mktempdir()
        try
            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                candidate_models=[[4]],
                delays=[10, 20, 30],
                derivative_points=4,
                order=2,
                out_dir=out_dir,
            )

            @test result.best_delays == [10, 30]
            @test result.best_score == 0.2
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "generated tau files honor user prefix directory" begin
        calls = []

        run_once = (; tau_file, out_fn, kwargs...) -> begin
            push!(calls, (tau_file=tau_file, out_fn=out_fn))
            n_rows = length(readlines(tau_file))
            return fake_errors(fill(1.0, n_rows))
        end

        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                candidate_models=[[1], [4]],
                delays=10:20:50,
                derivative_points=4,
                order=2,
                nr_delays=2,
                tau_file_prefix=prefix,
                tau_file_suffix="_run1",
            )

            @test result.artifacts_dir == prefix
            @test isdir(prefix)
            @test sort(basename.(call.tau_file for call in calls)) == [
                "TAU_ALL__1_0_run1",
                "TAU_ALL__2_1_run1",
            ]
            @test all(call -> dirname(call.tau_file) == prefix, calls)
            @test all(call -> dirname(call.out_fn) == prefix, calls)
            @test isfile(joinpath(prefix, "TAU_ALL__1_0_run1"))
            @test isfile(joinpath(prefix, "TAU_ALL__2_1_run1"))
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "one-delay tau-file models are re-encoded for nr_tau one" begin
        calls = []

        run_once = (; model, nr_tau, tau_file, out_fn, kwargs...) -> begin
            push!(calls, (
                model=Int[model...],
                nr_tau=nr_tau,
                tau_file=tau_file,
                out_fn=out_fn,
            ))
            n_rows = length(readlines(tau_file))
            return fake_errors(fill(1.0, n_rows))
        end

        out_dir = mktempdir()
        try
            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                N_MOD=3,
                DDAorder=3,
                delays=6:8,
                derivative_points=5,
                out_dir=out_dir,
            )

            reported_call = only(filter(call -> basename(call.out_fn) == "structure_selection_m7_d1", calls))
            @test reported_call.nr_tau == 1
            @test reported_call.model == [1, 2, 3]
            @test basename(reported_call.tau_file) == "TAU_ALL__1_0"
            @test result.trials[7].model == [1, 3, 6]
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "generated tau-file artifacts can be cleaned up on error" begin
        parent_dir = mktempdir()

        try
            @test_throws ErrorException DelayDifferentialAnalysis.StructureSelection._structure_selection(
                (; kwargs...) -> error("DDA failed");
                file_path="data.ascii",
                channels=[1],
                candidate_models=[[1]],
                delays=10:20,
                derivative_points=4,
                order=2,
                out_dir=parent_dir,
                cleanup_on_error=true,
            )

            @test isempty(readdir(parent_dir))
        finally
            rm(parent_dir; recursive=true, force=true)
        end
    end

    @testset "joint model scope uses all channels together" begin
        calls = []

        run_once = (; channels, model, delays, kwargs...) -> begin
            push!(calls, (channels=Int[channels...], model=Int[model...], delays=Int[delays...], kwargs=kwargs))
            return fake_result(model == [1, 2, 10] ? 0.1 : 0.5)
        end

        result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
            run_once;
            file_path="data.ascii",
            channels=[1, 2, 3],
            binary_path="/tmp/run_DDA_AsciiEdf",
            candidate_models=[[1, 2, 6], [1, 2, 10]],
            delays=[[7, 10]],
            derivative_points=4,
            order=3,
            model_scope=:joint,
        )

        @test result isa StructureSelectionResult
        @test result.best_model == [1, 2, 10]
        @test length(calls) == 2
        @test all(call -> call.channels == [1, 2, 3], calls)
        @test !haskey(calls[1].kwargs, :model_scope)
    end

    @testset "per-channel model scope selects one model per channel" begin
        calls = []
        scores = Dict(
            (1, [1, 2, 6]) => 0.1,
            (1, [1, 2, 10]) => 0.9,
            (2, [1, 2, 6]) => 0.8,
            (2, [1, 2, 10]) => 0.2,
        )

        run_once = (; channels, model, delays, out_fn, kwargs...) -> begin
            channel = only(channels)
            model_vector = Int[model...]
            push!(calls, (
                channels=Int[channels...],
                model=model_vector,
                delays=Int[delays...],
                out_fn=out_fn,
                kwargs=kwargs,
            ))
            return fake_result(scores[(channel, model_vector)])
        end

        out_dir = mktempdir()
        try
            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=1:2,
                binary_path="/tmp/run_DDA_AsciiEdf",
                candidate_models=[[1, 2, 6], [1, 2, 10]],
                delays=[[7, 10]],
                derivative_points=4,
                order=3,
                model_scope="per_channel",
                out_dir=out_dir,
            )

            @test result isa PerChannelStructureSelectionResult
            @test length(result.results) == 2
            @test result.results[1].channel_index == 1
            @test result.results[1].channel == 1
            @test result.results[1].selection.best_model == [1, 2, 6]
            @test result.results[1].selection.best_score == 0.1
            @test result.results[2].channel_index == 2
            @test result.results[2].channel == 2
            @test result.results[2].selection.best_model == [1, 2, 10]
            @test result.results[2].selection.best_score == 0.2
            @test length(calls) == 4
            @test calls[1].channels == [1]
            @test calls[3].channels == [2]
            @test basename(calls[1].out_fn) == "structure_selection_ch1_m1_d1"
            @test basename(calls[3].out_fn) == "structure_selection_ch2_m1_d1"
            @test !haskey(calls[1].kwargs, :model_scope)
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "rejects invalid structure-selection model scope" begin
        @test_throws ErrorException DelayDifferentialAnalysis.StructureSelection._structure_selection(
            (; kwargs...) -> fake_result(1.0);
            file_path="data.ascii",
            channels=[1],
            candidate_models=[[1, 2, 6]],
            delays=[[7, 10]],
            derivative_points=4,
            order=3,
            model_scope=:unknown,
        )
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

    @testset "make_MOD generates Claudia-style MOD libraries" begin
        MOD = make_MOD(1, 2)

        @test MOD == [
            1 0 0 0 0
            0 0 1 0 0
            0 0 0 1 0
        ]

        @test DelayDifferentialAnalysis.StructureSelection._p_dda(2; nr_delays=2) == [
            0 1
            0 2
            1 1
            1 2
            2 2
        ]

        @test_throws ErrorException make_MOD(1, 2; nr_delays=3)
        @test_throws MethodError structure_selection(1, 2)
    end

    @testset "prints MOD and P_DDA as terminal Unicode" begin
        MOD = make_MOD(1, 2)
        io = IOBuffer()

        print_structure_selection(io, MOD, 2)
        text = String(take!(io))

        @test occursin("P_DDA", text)
        @test occursin("1: [0, 1]", text)
        @test occursin("MOD rows × P_DDA terms", text)
        @test occursin("model | 1 | 2 | 3 | 4 | 5", text)
        @test occursin("1 | ✓ |   |   |   |", text)
        @test occursin("2 |   |   | ✓ |   |", text)
        @test occursin("Models", text)
        @test occursin("ẋ = a₁·x₁", text)
        @test occursin("ẋ = a₁·x₁²", text)
        @test occursin("ẋ = a₁·x₁·x₂", text)
        @test !occursin("\\dot{x}", text)
    end

    @testset "write_model_LaTeX remains available explicitly" begin
        MOD = make_MOD(1, 2)
        P_DDA = DelayDifferentialAnalysis.StructureSelection._p_dda(2)
        io = IOBuffer()

        write_model_LaTeX(io, MOD, nothing, P_DDA, 3)
        text = String(take!(io))

        @test occursin("\\dot{x}", text)
        @test occursin("a_1 x_1 x_2", text)
    end

    @testset "structure_selection can generate candidate models from make_MOD" begin
        calls = []
        scores = Dict(
            ([1], [7, 10]) => 0.7,
            ([1], [10, 20]) => 0.6,
            ([3], [7, 10]) => 0.5,
            ([3], [10, 20]) => 0.4,
            ([4], [7, 10]) => 0.3,
            ([4], [10, 20]) => 0.1,
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
            channels=[1],
            binary_path="/tmp/run_DDA_AsciiEdf",
            N_MOD=1,
            DDAorder=2,
            delays=[[7, 10], [10, 20]],
            derivative_points=4,
            input_format=:ascii,
        )

        @test result.best_model == [4]
        @test result.best_delays == [10, 20]
        @test result.best_score == 0.1
        @test length(result.trials) == 6
        @test calls[1].model == [1]
        @test calls[1].kwargs[:order] == 2
        @test calls[1].kwargs[:nr_tau] == 2
    end

    @testset "candidate_delays remains a deprecated alias" begin
        calls = []

        run_once = (; model, delays, kwargs...) -> begin
            push!(calls, (model=Int[model...], delays=Int[delays...], kwargs=kwargs))
            return fake_result(delays == [32, 9] ? 0.2 : 0.9)
        end

        result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
            run_once;
            file_path="data.ascii",
            channels=[1],
            binary_path="/tmp/run_DDA_AsciiEdf",
            candidate_models=[[1, 2, 6]],
            candidate_delays=[[7, 10], [32, 9]],
            derivative_points=4,
            order=3,
        )

        @test result.best_delays == [32, 9]
        @test length(calls) == 2
    end

    @testset "structure_selection rejects tau_file row expansion" begin
        tau_path = tempname()
        write(tau_path, "7 10\n32 9\n")

        try
            @test_throws ErrorException DelayDifferentialAnalysis.StructureSelection._structure_selection(
                (; kwargs...) -> fake_result(0.1);
                file_path="data.ascii",
                channels=[1],
                binary_path="/tmp/run_DDA_AsciiEdf",
                candidate_models=[[1, 2, 6]],
                tau_file=tau_path,
                derivative_points=4,
                order=3,
            )
        finally
            rm(tau_path; force=true)
        end
    end
end
