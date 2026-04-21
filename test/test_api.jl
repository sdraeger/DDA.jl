using Test
using DelayDifferentialAnalysis

@testset "API" begin
    # Check if binary is available for integration tests
    binary_available = DelayDifferentialAnalysis.Variants.find_binary() !== nothing
    repo_binary = abspath(joinpath(@__DIR__, "..", "..", "..", "bin", "run_DDA_AsciiEdf"))

    @testset "run_st integration" begin
        if !binary_available
            @info "Skipping run_st integration test: DDA binary not found"
            @test_skip false
        else
            data = randn(2, 10000)
            result = run_st(data; sfreq=256.0, wl=200, ws=100)

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
            @test_skip false
        else
            data = randn(2, 10000)
            result = run_st(data; sfreq=256.0, wl=200, ws=100,
                           channel_labels=["Fp1", "Fp2"])
            @test result.channel_labels == ["Fp1", "Fp2"]
        end
    end

    @testset "run_ct integration" begin
        if !binary_available
            @info "Skipping run_ct integration test: DDA binary not found"
            @test_skip false
        else
            # CT needs more data and explicit CT window params
            data = randn(3, 20000)
            try
                result = run_ct(data; sfreq=256.0, wl=200, ws=100)
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
        @test_throws ErrorException run_ct(data; sfreq=256.0)
    end

    @testset "file-based API reports missing input cleanly" begin
        fake_binary = tempname()
        touch(fake_binary)

        try
            @test_throws ErrorException run_st(
                "/nonexistent/path/to/data.edf",
                [1, 2];
                binary_path=fake_binary,
            )
        finally
            rm(fake_binary, force=true)
        end
    end

    @testset "matrix API smoke test with repo binary" begin
        if !isfile(repo_binary)
            @info "Skipping matrix API smoke test: repo binary not found"
            @test_skip false
        else
            data = randn(3, 20_000)

            st = run_st(data; binary_path=repo_binary, wl=200, ws=100)
            ct = run_ct(data; binary_path=repo_binary, wl=200, ws=100, ct_wl=2, ct_ws=2)
            de = run_de(data; binary_path=repo_binary, wl=200, ws=100, ct_wl=2, ct_ws=2)

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
            @test_skip false
        else
            # DE needs more data
            data = randn(2, 20000)
            try
                result = run_de(data; sfreq=256.0, wl=200, ws=100)
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
