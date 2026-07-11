using Test
using DelayDifferentialAnalysis

@testset "Batch" begin
    @testset "GroupResult construction" begin
        coeffs = randn(3, 2, 5, 3)  # 3 subjects, 2 channels, 5 windows, 3 coeffs
        errs = randn(3, 2, 5)
        labels = ["subj1", "subj2", "subj3"]
        ch_labels = ["ch1", "ch2"]
        params = Dict{String,Any}("WL" => 200)

        g = GroupResult(coeffs, errs, labels, ch_labels, params, "ST")

        @test n_subjects(g) == 3
        @test n_channels(g) == 2
        @test n_windows(g) == 5
        @test n_coeffs(g) == 3
        @test g.variant == "ST"
    end

    @testset "mean_over_windows" begin
        coeffs = ones(2, 3, 4, 2)  # all ones
        errs = zeros(2, 3, 4)
        g = GroupResult(coeffs, errs, ["s1", "s2"], ["c1", "c2", "c3"], Dict{String,Any}(), "ST")

        m = mean_over_windows(g)
        @test size(m) == (2, 3, 2)
        @test all(m .== 1.0)
    end

    @testset "collect_results STResult" begin
        n_k = 3
        r1 = STResult(randn(2, 5, n_k), rand(2, 5), collect(1:5), collect(2:6), ["a", "b"], Dict{String,Any}())
        r2 = STResult(randn(2, 5, n_k), rand(2, 5), collect(1:5), collect(2:6), ["a", "b"], Dict{String,Any}())
        r3 = STResult(randn(2, 5, n_k), rand(2, 5), collect(1:5), collect(2:6), ["a", "b"], Dict{String,Any}())

        g = collect_results([r1, r2, r3]; labels=["s1", "s2", "s3"])

        @test n_subjects(g) == 3
        @test n_channels(g) == 2
        @test n_windows(g) == 5
        @test n_coeffs(g) == n_k
        @test g.variant == "ST"
        @test g.subject_labels == ["s1", "s2", "s3"]
    end

    @testset "collect_results truncates windows" begin
        n_k = 3
        r1 = STResult(randn(2, 10, n_k), rand(2, 10), collect(1:10), collect(2:11), ["a", "b"], Dict{String,Any}())
        r2 = STResult(randn(2, 5, n_k), rand(2, 5), collect(1:5), collect(2:6), ["a", "b"], Dict{String,Any}())

        g = collect_results([r1, r2])

        @test n_windows(g) == 5  # truncated to minimum
    end

    @testset "collect_results CTResult" begin
        n_k = 3
        r1 = CTResult(randn(3, 5, n_k), rand(3, 5), collect(1:5), collect(2:6), ["a-b", "a-c", "b-c"], Dict{String,Any}())
        r2 = CTResult(randn(3, 5, n_k), rand(3, 5), collect(1:5), collect(2:6), ["a-b", "a-c", "b-c"], Dict{String,Any}())

        g = collect_results([r1, r2])

        @test n_subjects(g) == 2
        @test n_channels(g) == 3  # pairs
        @test g.variant == "CT"
    end

    @testset "collect_results DEResult" begin
        r1 = DEResult(randn(8), collect(1:8), collect(2:9), Dict{String,Any}())
        r2 = DEResult(randn(6), collect(1:6), collect(2:7), Dict{String,Any}())

        g = collect_results([r1, r2])

        @test n_subjects(g) == 2
        @test n_windows(g) == 6  # truncated
        @test g.variant == "DE"
    end

    @testset "collect_results default labels" begin
        r1 = STResult(randn(1, 3, 2), rand(1, 3), [1, 2, 3], [2, 3, 4], ["ch1"], Dict{String,Any}())
        r2 = STResult(randn(1, 3, 2), rand(1, 3), [1, 2, 3], [2, 3, 4], ["ch1"], Dict{String,Any}())

        g = collect_results([r1, r2])
        @test g.subject_labels == ["subj_1", "subj_2"]
    end

    @testset "collect_results mismatched types" begin
        r1 = STResult(randn(2, 3, 3), rand(2, 3), [1, 2, 3], [2, 3, 4], ["a", "b"], Dict{String,Any}())
        r2 = DEResult(randn(3), [1, 2, 3], [2, 3, 4], Dict{String,Any}())

        @test_throws ErrorException collect_results([r1, r2])
    end

    @testset "collect_results empty" begin
        @test_throws ErrorException collect_results(Any[])
    end

    @testset "collect_results channel mismatch" begin
        r1 = STResult(randn(2, 3, 3), rand(2, 3), [1, 2, 3], [2, 3, 4], ["a", "b"], Dict{String,Any}())
        r2 = STResult(randn(3, 3, 3), rand(3, 3), [1, 2, 3], [2, 3, 4], ["a", "b", "c"], Dict{String,Any}())

        @test_throws ErrorException collect_results([r1, r2])
    end

    @testset "run_batch dispatches through the keyword-only API" begin
        directory = mktempdir()
        data_path = joinpath(directory, "data.ascii")
        binary_path = joinpath(directory, Sys.iswindows() ? "dda.exe" : "dda")
        write(data_path, "1 2\n3 4\n")
        touch(binary_path)

        try
            caught = try
                run_batch([data_path]; binary_path=binary_path, progress=false)
                nothing
            catch error
                error
            end
            @test caught !== nothing
            @test !(caught isa MethodError)
        finally
            rm(directory; recursive=true, force=true)
        end
    end
end
