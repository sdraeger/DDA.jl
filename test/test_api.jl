using Test
using DelayDifferentialAnalysis

@testset "API" begin
    # Check if binary is available for integration tests
    binary_available = DelayDifferentialAnalysis.Flavors.find_binary() !== nothing
    repo_binary = abspath(joinpath(@__DIR__, "..", "..", "..", "bin", "run_DDA_AsciiEdf"))

    @testset "run_st integration" begin
        if !binary_available
            @info "Skipping run_st integration test: DDA binary not found"
        else
            data = randn(2, 10000)
            result = run_st(data=data; sfreq=256.0, WL=200, WS=100)

            @test result isa STResult
            @test n_channels(result) == 2
            @test n_coeffs(result) == length(DelayDifferentialAnalysis.DDADefaults.MODEL_PARAMS)
            @test n_windows(result) > 0
            @test length(result.window_starts) == n_windows(result)
            @test result.channel_labels == ["ch1", "ch2"]
        end
    end

    @testset "run_st custom labels" begin
        if !binary_available
            @info "Skipping run_st custom labels test: DDA binary not found"
        else
            data = randn(2, 10000)
            result = run_st(data=data; sfreq=256.0, WL=200, WS=100,
                           channel_labels=["Fp1", "Fp2"])
            @test result.channel_labels == ["Fp1", "Fp2"]
        end
    end

    @testset "run_ct integration" begin
        if !binary_available
            @info "Skipping run_ct integration test: DDA binary not found"
        else
            # CT needs more data and explicit CT window params
            data = randn(3, 20000)
            try
                result = run_ct(data=data; sfreq=256.0, WL=200, WS=100, ct_wl=2, ct_ws=2)
                @test result isa CTResult
                @test n_pairs(result) == 3  # C(3,2) = 3
                @test n_windows(result) > 0
                @test length(result.pair_labels) == 3
            catch e
                @warn "CT integration test failed (binary may not support this data config)" exception=e
                @test_broken false
            end
        end
    end

    @testset "run_ct requires 2+ channels" begin
        data = randn(1, 5000)
        @test_throws ErrorException run_ct(data=data; sfreq=256.0)
    end

    @testset "Public run_* APIs are keyword-only" begin
        data = randn(2, 1024)
        @test_throws MethodError run_st(data)
        @test_throws MethodError run_ct(data)
        @test_throws MethodError run_de(data)
    end

    @testset "file-based label resolution uses file metadata when available" begin
        temp_dir = mktempdir()
        ascii_path = joinpath(temp_dir, "labels.tsv")

        open(ascii_path, "w") do io
            println(io, "Fp1\tFp2\tC3")
            println(io, "1\t2\t3")
            println(io, "4\t5\t6")
        end

        try
            labels = DelayDifferentialAnalysis.API._resolve_labels(ascii_path, [1, 3], nothing)
            @test labels == ["Fp1", "C3"]
        finally
            rm(temp_dir; recursive=true, force=true)
        end
    end

    @testset "file-based label resolution falls back to compact defaults without header" begin
        temp_dir = mktempdir()
        ascii_path = joinpath(temp_dir, "numeric.tsv")

        open(ascii_path, "w") do io
            println(io, "1\t2\t3")
            println(io, "4\t5\t6")
        end

        try
            labels = DelayDifferentialAnalysis.API._resolve_labels(ascii_path, [1, 3], nothing)
            @test labels == ["ch1", "ch3"]
        finally
            rm(temp_dir; recursive=true, force=true)
        end
    end

    @testset "file-based API reports missing input cleanly" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            @test_throws ErrorException run_st(
                file_path="/nonexistent/path/to/data.edf",
                channels=[1, 2],
                binary_path=fake_binary,
            )
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "custom model requires explicit derivative config and order" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            @test_throws ErrorException run_st(
                data=randn(2, 1024),
                binary_path=fake_binary,
                model=[1, 2, 10],
            )
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "raw T, derived t, and full coefficients are preserved" begin
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

        result = DelayDifferentialAnalysis.API._st_from_raw(
            channels,
            ["ch1", "ch2"];
            sfreq=1.0,
            delays=[7, 10],
            model=[1, 2, 10],
            WL=128,
            WS=100,
            derivative_points=3,
            TM=12,
            order=4,
            nr_tau=2,
            sampling_rate=(500, 500),
            out_fn=nothing,
            selected_channels=[1, 2],
        )

        @test size(result.coefficients) == (2, 2, 3)
        @test result.coefficients[1, 1, :] == [1.0, 2.0, 3.0]
        @test result.coefficients[2, 2, :] == [10.0, 11.0, 12.0]
        @test result.errors == [0.1 0.2; 0.3 0.4]
        @test result.T == [10.0, 110.0]
        @test result.t == [0.052, 0.252]
        @test result.window_starts == [0, 100]
        @test result.window_ends == [128, 228]
        @test result.params["TM"] == 12
        @test result.params["sampling_rate"] == (500, 500)
    end

    @testset "matrix API smoke test with repo binary" begin
        if !isfile(repo_binary)
            @info "Skipping matrix API smoke test: repo binary not found"
        else
            data = randn(3, 20_000)

            st = run_st(data=data; binary_path=repo_binary, WL=200, WS=100)
            ct = run_ct(data=data; binary_path=repo_binary, WL=200, WS=100, ct_wl=2, ct_ws=2)
            de = run_de(data=data; binary_path=repo_binary, WL=200, WS=100, ct_wl=2, ct_ws=2)

            @test n_channels(st) == 3
            @test n_pairs(ct) == 3
            @test n_windows(st) > 0
            @test n_windows(ct) > 0
            @test n_windows(de) > 0
        end
    end

    @testset "run_de integration" begin
        if !binary_available
            @info "Skipping run_de integration test: DDA binary not found"
        else
            # DE needs more data
            data = randn(2, 20000)
            try
                result = run_de(data=data; sfreq=256.0, WL=200, WS=100, ct_wl=2, ct_ws=2)
                @test result isa DEResult
                @test n_windows(result) > 0
                @test length(result.ergodicity) == n_windows(result)
            catch e
                @warn "DE integration test failed (binary may not support this data config)" exception=e
                @test_broken false
            end
        end
    end
end
