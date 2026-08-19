using Random
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

    function write_structure_st_series(
        path,
        errors::AbstractArray{<:Real,3},
        n_terms::Integer;
        T=collect(0:(size(errors, 3) - 1)),
    )
        n_channels, n_tau, _ = size(errors)
        open(path, "w") do io
            for window_idx in axes(errors, 3)
                fields = zeros(Float64, Int(n_terms) + 1, n_channels, n_tau)
                fields[end, :, :] .= Float64.(errors[:, :, window_idx])
                println(io, "$(T[window_idx]) $(T[window_idx] + 1) ", join(vec(fields), " "))
            end
        end
        return nothing
    end

    function write_structure_st(path, errors::AbstractMatrix{<:Real}, n_terms::Integer)
        values = reshape(Float64.(errors), size(errors, 1), size(errors, 2), 1)
        return write_structure_st_series(path, values, n_terms)
    end

    @testset "structure_selection_compute computes outputs without selecting" begin
        calls = []
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            run_once = (; kwargs...) -> begin
                push!(calls, kwargs)
                return nothing
            end

            run = DelayDifferentialAnalysis.StructureSelection._structure_selection_compute(
                run_once;
                file_path="data.ascii",
                binary_path="/tmp/run_DDA_AsciiEdf",
                input_format=:ascii,
                channels=nothing,
                N_MOD=1,
                DDAorder=2,
                delays=10:20,
                derivative_points=4,
                prefix=prefix,
                MOD_numbers=[1, 2],
                randomize=false,
            )

            @test run.prefix == prefix
            @test size(run.MOD, 1) == 3
            @test run.model_numbers == [1, 2]
            @test run.channels === nothing
            @test length(calls) == 2
            @test all(call -> call[:load_results] == false, calls)
            @test all(call -> call[:channels] === nothing, calls)
            @test sort(basename.(call[:out_fn] for call in calls)) == [
                "structure_selection_01",
                "structure_selection_03",
            ]
            @test isfile(joinpath(prefix, "TAU_ALL__1_0"))
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "structure_selection_select reads cached outputs and defaults to all models and channels" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            MOD = make_MOD(1, 2)
            write(joinpath(prefix, "TAU_ALL__1_0"), "10\n20\n")
            write_structure_st(
                joinpath(prefix, "structure_selection_01_ST"),
                [0.1 0.4; 0.9 0.8],
                1,
            )
            write_structure_st(
                joinpath(prefix, "structure_selection_03_ST"),
                [0.3 0.5; 0.2 0.1],
                1,
            )

            run = StructureSelectionRun(
                prefix,
                MOD,
                2,
                2,
                [10, 20],
                nothing,
                [1, 2],
                4,
                "",
                "structure_selection",
            )

            selected_channel = structure_selection_select(run; channels=[2])
            @test selected_channel.best_model == [3]
            @test selected_channel.best_delays == [20]
            @test selected_channel.best_score == 0.1

            selected_all = structure_selection_select(run)
            @test selected_all.best_model == [3]
            @test selected_all.best_delays == [10]
            @test selected_all.best_score == 0.25

            selected_model = structure_selection_select(run; models=[1])
            @test selected_model.best_model == [1]
            @test length(selected_model.trials) == 1
            @test_throws ErrorException structure_selection_select(
                run;
                models=[1],
                MOD_numbers=[1],
            )
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "structure_selection_select maps requested channel names from compute run" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            MOD = make_MOD(1, 2)
            write(joinpath(prefix, "TAU_ALL__1_0"), "10\n20\n")
            write_structure_st(
                joinpath(prefix, "structure_selection_01_ST"),
                [0.1 0.4; 0.9 0.8],
                1,
            )
            write_structure_st(
                joinpath(prefix, "structure_selection_03_ST"),
                [0.3 0.5; 0.2 0.1],
                1,
            )

            run = StructureSelectionRun(
                prefix,
                MOD,
                2,
                2,
                [10, 20],
                [10, 20],
                [1, 2],
                4,
                "",
                "structure_selection",
            )

            selected = structure_selection_select(run; channel=[20])
            @test selected.best_model == [3]
            @test selected.best_delays == [20]
            @test selected.best_score == 0.1

            selected_all = structure_selection_select(run)
            @test selected_all.best_model == [3]
            @test selected_all.best_delays == [10]
            @test selected_all.best_score == 0.25

            per_channel = structure_selection_select(run; model_scope=:per_channel)
            @test [result.channel for result in per_channel.results] == [10, 20]
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "structure-selection plot data preserves delay-pair and time winners" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            MOD = make_MOD(2, 2)
            write(joinpath(prefix, "TAU_ALL__2_1"), "10 20\n10 30\n20 30\n")
            write(
                joinpath(prefix, "TAU_ALL__2_0"),
                "10 20\n10 30\n20 10\n20 30\n30 10\n30 20\n",
            )

            symmetric_errors = Array{Float64}(undef, 1, 3, 2)
            symmetric_errors[1, :, :] = [0.2 5.0; 5.0 5.0; 5.0 0.3]
            asymmetric_errors = Array{Float64}(undef, 1, 6, 2)
            asymmetric_errors[1, :, :] = [
                2.0 2.0
                1.5 1.5
                0.1 4.0
                2.0 2.0
                3.0 3.0
                4.0 4.0
            ]
            write_structure_st_series(
                joinpath(prefix, "structure_selection_01_02_ST"),
                symmetric_errors,
                2;
                T=[100, 200],
            )
            write_structure_st_series(
                joinpath(prefix, "structure_selection_01_04_ST"),
                asymmetric_errors,
                2;
                T=[100, 200],
            )

            run = StructureSelectionRun(
                prefix,
                MOD,
                2,
                2,
                [10, 20, 30],
                [1],
                [1, 3],
                4,
                "",
                "structure_selection",
            )
            module_ref = DelayDifferentialAnalysis.StructureSelection

            all_data = module_ref._structure_selection_plot_data(run; mode=:all)
            @test all_data.tau1 == [10, 20, 30]
            @test all_data.tau2 == [10, 20, 30]
            @test all(all_data.model_numbers[idx, idx] == 0 for idx in 1:3)
            @test all_data.model_numbers[2, 3] == 1
            @test all_data.model_numbers[3, 2] == 3

            time_data = module_ref._structure_selection_plot_data(run; mode=:time)
            @test time_data.T == [100.0, 200.0]
            @test time_data.model_numbers == [3, 1]
            @test time_data.tau1 == [20, 20]
            @test time_data.tau2 == [10, 30]

            one_model = module_ref._structure_selection_plot_data(
                run;
                mode=:time,
                models=[1],
            )
            @test one_model.model_numbers == [1, 1]
            @test one_model.tau1 == [10, 20]
            @test one_model.tau2 == [20, 30]
            @test_throws ErrorException module_ref._structure_selection_plot_data(run; mode=:other)
            @test_throws ErrorException module_ref._structure_selection_plot_data(
                run;
                models=[2],
            )

            write_structure_st_series(
                joinpath(prefix, "structure_selection_01_04_ST"),
                asymmetric_errors,
                2;
                T=[100, 201],
            )
            @test_throws ErrorException module_ref._structure_selection_plot_data(run; mode=:time)
        finally
            rm(out_dir; recursive=true, force=true)
        end
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
            calls_by_model = Dict(Tuple(call.model) => call for call in calls)
            @test dirname(calls_by_model[(1,)].tau_file) == result.artifacts_dir
            @test basename(result.artifacts_dir) |> startswith("structure_selection_")
            @test sort(basename.(calls[i].tau_file for i in eachindex(calls))) == [
                "TAU_ALL__1_0_run1",
                "TAU_ALL__2_1_run1",
            ]
            @test read(calls_by_model[(1,)].tau_file, String) == "10\n30\n50\n"
            @test read(calls_by_model[(4,)].tau_file, String) == "10 30\n10 50\n30 50\n"
            @test calls_by_model[(1,)].delays == [10]
            @test calls_by_model[(1,)].nr_tau == 1
            @test calls_by_model[(4,)].delays == [10, 30]
            @test calls_by_model[(4,)].nr_tau == 2
            @test calls_by_model[(1,)].kwargs[:flavors] == ["ST"]
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

    @testset "generated tau-file selection aggregates channel-stacked error rows" begin
        run_once = (; kwargs...) -> fake_errors(reshape([0.5, 0.1, 0.7, 0.4, 0.3, 0.9], 6, 1))

        out_dir = mktempdir()
        try
            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1, 2],
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
                prefix=prefix,
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
            @test sort(basename.(call.out_fn for call in calls)) == [
                "structure_selection_01",
                "structure_selection_04",
            ]
            @test isfile(joinpath(prefix, "TAU_ALL__1_0_run1"))
            @test isfile(joinpath(prefix, "TAU_ALL__2_1_run1"))
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "existing pool output is reused without rerunning" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            write(
                joinpath(prefix, "structure_selection_04_ST"),
                "0 0 1 2 3 0.5\n0 0 1 2 3 0.2\n0 0 1 2 3 0.9\n",
            )
            calls = Ref(0)
            run_once = (; kwargs...) -> begin
                calls[] += 1
                error("structure selection should reuse the existing output")
            end

            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                candidate_models=[[4]],
                delays=10:12,
                derivative_points=4,
                order=2,
                prefix=prefix,
            )

            @test calls[] == 0
            @test result.best_delays == [10, 12]
            @test result.best_score == 0.2
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "existing pool output does not rewrite tau file" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            tau_path = joinpath(prefix, "TAU_ALL__2_1")
            write(tau_path, "sentinel\n")
            write(
                joinpath(prefix, "structure_selection_04_ST"),
                "0 0 1 2 3 0.5\n0 0 1 2 3 0.2\n0 0 1 2 3 0.9\n",
            )
            run_once = (; kwargs...) -> error("structure selection should reuse the existing output")

            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                candidate_models=[[4]],
                delays=10:12,
                derivative_points=4,
                order=2,
                prefix=prefix,
            )

            @test read(tau_path, String) == "sentinel\n"
            @test result.best_delays == [10, 12]
            @test result.best_score == 0.2
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "identical tau file writes preserve existing file" begin
        out_dir = mktempdir()
        try
            tau_path = joinpath(out_dir, "TAU_ALL__2_0")
            tau_rows = [[10, 11], [10, 12], [11, 10]]
            DelayDifferentialAnalysis.StructureSelection._write_tau_file(tau_path, tau_rows)
            before = stat(tau_path).mtime
            sleep(1.1)
            DelayDifferentialAnalysis.StructureSelection._write_tau_file(tau_path, tau_rows)

            @test stat(tau_path).mtime == before
            @test read(tau_path, String) == "10 11\n10 12\n11 10\n"
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "partial existing pool output is skipped without overwriting" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            st_path = joinpath(prefix, "structure_selection_04_ST")
            write(st_path, "0 0 1 2 3 0.5\n")
            calls = Ref(0)
            run_once = (; kwargs...) -> begin
                calls[] += 1
                return fake_errors([0.5, 0.2, 0.9])
            end

            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=1:11,
                candidate_models=[[4], [5]],
                delays=10:12,
                derivative_points=4,
                order=2,
                prefix=prefix,
                randomize=false,
            )

            @test calls[] == 1
            @test result.best_model == [5]
            @test result.best_delays == [11]
            @test result.best_score == 0.2
            @test read(st_path, String) == "0 0 1 2 3 0.5\n"
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "all conflicted pool outputs give one clear error" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            st_path = joinpath(prefix, "structure_selection_04_ST")
            write(st_path, "0 0 1 2 3 0.5\n")
            calls = Ref(0)
            run_once = (; kwargs...) -> begin
                calls[] += 1
                return fake_errors([0.5, 0.2, 0.9])
            end

            err = try
                DelayDifferentialAnalysis.StructureSelection._structure_selection(
                    run_once;
                    file_path="data.ascii",
                    channels=1:11,
                    candidate_models=[[4]],
                    delays=10:12,
                    derivative_points=4,
                    order=2,
                    prefix=prefix,
                )
                nothing
            catch err
                err
            end

            @test err isa ErrorException
            @test occursin("No usable structure-selection candidates", sprint(showerror, err))
            @test calls[] == 0
            @test read(st_path, String) == "0 0 1 2 3 0.5\n"
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "existing info file without ST is skipped without overwriting" begin
        out_dir = mktempdir()
        try
            prefix = joinpath(out_dir, "RUN1")
            mkpath(prefix)
            info_path = joinpath(prefix, "structure_selection_04.info")
            write(info_path, "sentinel\n")
            calls = Ref(0)
            run_once = (; kwargs...) -> begin
                calls[] += 1
                return fake_errors([0.5, 0.2, 0.9])
            end

            result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                candidate_models=[[4], [5]],
                delays=10:12,
                derivative_points=4,
                order=2,
                prefix=prefix,
                randomize=false,
            )

            @test calls[] == 1
            @test result.best_model == [5]
            @test result.best_delays == [11]
            @test result.best_score == 0.2
            @test read(info_path, String) == "sentinel\n"
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "pool output replaces stale same-host lock" begin
        if Sys.iswindows()
            @test_skip false
        else
            out_dir = mktempdir()
            try
                output_base = joinpath(out_dir, "structure_selection_04")
                lock_path = "$(output_base).lock"
                mkdir(lock_path)
                write(
                    joinpath(lock_path, "owner"),
                    "pid=999999999\nhost=$(gethostname())\n",
                )
                calls = Ref(0)

                task = @async DelayDifferentialAnalysis.StructureSelection._run_or_reuse_pool_output(
                    output_base,
                    3,
                ) do
                    calls[] += 1
                    return fake_errors([0.5, 0.2, 0.9])
                end
                sleep(0.2)

                @test istaskdone(task)
                result = fetch(task)
                @test calls[] == 1
                @test DelayDifferentialAnalysis.StructureSelection._score_result(result, :minimum_error) == 0.2
                @test !ispath(lock_path)
            finally
                rm(out_dir; recursive=true, force=true)
            end
        end
    end

    @testset "ownerless lock stale threshold protects new locks" begin
        out_dir = mktempdir()
        try
            lock_path = joinpath(out_dir, "structure_selection_04.lock")
            mkdir(lock_path)
            mtime = stat(lock_path).mtime
            grace = DelayDifferentialAnalysis.StructureSelection._POOL_OWNERLESS_LOCK_GRACE_SECONDS

            @test !DelayDifferentialAnalysis.StructureSelection._ownerless_lock_is_stale(
                lock_path,
                mtime + grace,
            )
            @test DelayDifferentialAnalysis.StructureSelection._ownerless_lock_is_stale(
                lock_path,
                mtime + grace + 1,
            )
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "pool waiter reruns when lock disappears without output" begin
        out_dir = mktempdir()
        try
            output_base = joinpath(out_dir, "structure_selection_04")
            lock_path = "$(output_base).lock"
            mkdir(lock_path)
            calls = Ref(0)
            remover = @async begin
                sleep(0.1)
                rm(lock_path; recursive=true, force=true)
            end

            result = DelayDifferentialAnalysis.StructureSelection._run_or_reuse_pool_output(
                output_base,
                3,
            ) do
                calls[] += 1
                return fake_errors([0.5, 0.2, 0.9])
            end
            wait(remover)

            @test calls[] == 1
            @test DelayDifferentialAnalysis.StructureSelection._score_result(result, :minimum_error) == 0.2
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "pool candidate order can be randomized" begin
        out_dir = mktempdir()
        try
            calls = String[]
            run_once = (; out_fn, tau_file, kwargs...) -> begin
                push!(calls, basename(out_fn))
                return fake_errors(fill(1.0, length(readlines(tau_file))))
            end
            models = [[1], [2], [3], [4], [5]]
            expected = [
                "structure_selection_01",
                "structure_selection_02",
                "structure_selection_03",
                "structure_selection_04",
                "structure_selection_05",
            ]

            DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                candidate_models=models,
                delays=10:12,
                derivative_points=4,
                order=2,
                out_dir=out_dir,
                randomize=true,
                rng=Random.MersenneTwister(7),
            )

            @test sort(calls) == expected
            @test calls != expected
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "structure-selection output names use padded model term indices" begin
        calls = []

        run_once = (; out_fn, kwargs...) -> begin
            push!(calls, out_fn)
            return fake_result(1.0)
        end

        out_dir = mktempdir()
        try
            DelayDifferentialAnalysis.StructureSelection._structure_selection(
                run_once;
                file_path="data.ascii",
                channels=[1],
                candidate_models=[[1, 2, 10]],
                delays=[[7, 10]],
                derivative_points=4,
                order=3,
                out_dir=out_dir,
            )

            @test basename(only(calls)) == "structure_selection_01_02_10_d1"
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

            reported_call = only(filter(call -> basename(call.out_fn) == "structure_selection_01_03_06", calls))
            @test reported_call.nr_tau == 1
            @test reported_call.model == [1, 2, 3]
            @test basename(reported_call.tau_file) == "TAU_ALL__1_0"
            @test any(trial -> trial.model == [1, 3, 6], result.trials)
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

    @testset "omitted channels runs structure selection across all binary channels" begin
        calls = []

        run_once = (; channels, model, delays, kwargs...) -> begin
            push!(calls, (channels=channels, model=Int[model...], delays=Int[delays...]))
            return fake_result(1.0)
        end

        result = DelayDifferentialAnalysis.StructureSelection._structure_selection(
            run_once;
            file_path="data.ascii",
            candidate_models=[[1, 2, 6]],
            delays=[[7, 10]],
            derivative_points=4,
            order=3,
        )

        @test result.best_model == [1, 2, 6]
        @test length(calls) == 1
        @test calls[1].channels === nothing
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
            @test basename(calls[1].out_fn) == "structure_selection_ch1_01_02_06_d1"
            @test basename(calls[3].out_fn) == "structure_selection_ch2_01_02_06_d1"
            @test !haskey(calls[1].kwargs, :model_scope)
        finally
            rm(out_dir; recursive=true, force=true)
        end
    end

    @testset "tau_file_prefix is not a structure-selection keyword" begin
        @test_throws ErrorException DelayDifferentialAnalysis.StructureSelection._structure_selection(
            (; kwargs...) -> fake_result(1.0);
            file_path="data.ascii",
            channels=[1],
            candidate_models=[[1, 2, 6]],
            delays=[[7, 10]],
            derivative_points=4,
            order=3,
            tau_file_prefix="/tmp/old",
        )
        out_dir = mktempdir()
        try
            @test_throws ErrorException DelayDifferentialAnalysis.StructureSelection._structure_selection_compute(
                (; kwargs...) -> nothing;
                file_path="data.ascii",
                channels=[1],
                N_MOD=1,
                DDAorder=2,
                delays=7:10,
                derivative_points=4,
                prefix=out_dir,
                tau_file_prefix="/tmp/old",
            )
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
            @test basename(seen[1].out_fn) == "structure_selection_01_02_06_d1"
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

        print_structure_selection(io, MOD)
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
        @test_throws ErrorException print_structure_selection(zeros(Int, 1, 4))
        @test_throws MethodError print_structure_selection(MOD, 2)

        padded_io = IOBuffer()
        print_structure_selection(padded_io, repeat([1 0 0 0 0], 12, 1))
        table_lines = filter(line -> occursin('|', line), split(String(take!(padded_io)), '\n'))
        @test length(table_lines) == 13
        @test length(unique(findfirst('|', line) for line in table_lines)) == 1
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
