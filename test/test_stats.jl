using Test
using DelayDifferentialAnalysis

@testset "Stats" begin
    # Helper to create mock STResults with controlled coefficients
    function make_st(coeffs_mean::Float64, n_ch::Int=2, n_win::Int=5, n_k::Int=3)
        coeffs = fill(coeffs_mean, n_ch, n_win, n_k) .+ 0.01 * randn(n_ch, n_win, n_k)
        errs = rand(n_ch, n_win)
        STResult(coeffs, errs, collect(1:n_win), collect(2:n_win+1),
                 ["ch$i" for i in 1:n_ch], Dict{String,Any}())
    end

    @testset "permutation_test basic" begin
        group_a = [make_st(1.0) for _ in 1:5]
        group_b = [make_st(0.0) for _ in 1:5]

        result = permutation_test(group_a, group_b; n_permutations=100, seed=42)

        @test size(result.observed_stat) == (2, 3)
        @test size(result.p_value) == (2, 3)
        @test size(result.null_distribution) == (100, 2, 3)
        @test result.n_permutations == 100
        @test result.tail == 0
        # With such a large difference, p-values should be small
        @test all(result.p_value .< 0.1)
    end

    @testset "permutation_test seed reproducibility" begin
        group_a = [make_st(1.0) for _ in 1:5]
        group_b = [make_st(0.0) for _ in 1:5]

        r1 = permutation_test(group_a, group_b; n_permutations=50, seed=123)
        r2 = permutation_test(group_a, group_b; n_permutations=50, seed=123)

        @test r1.null_distribution == r2.null_distribution
        @test r1.p_value == r2.p_value
    end

    @testset "permutation_test tail options" begin
        group_a = [make_st(1.0) for _ in 1:5]
        group_b = [make_st(0.0) for _ in 1:5]

        r_two = permutation_test(group_a, group_b; n_permutations=100, tail=0, seed=42)
        r_greater = permutation_test(group_a, group_b; n_permutations=100, tail=1, seed=42)
        r_less = permutation_test(group_a, group_b; n_permutations=100, tail=-1, seed=42)

        # group_a > group_b, so "greater" should have small p, "less" should have large p
        @test all(r_greater.p_value .< 0.1)
        @test all(r_less.p_value .> 0.5)
    end

    @testset "permutation_test identical groups" begin
        group = [make_st(1.0) for _ in 1:5]

        result = permutation_test(group, group; n_permutations=100, seed=42)

        # p-values should not be significant for identical groups
        @test all(result.p_value .> 0.05)
    end

    @testset "compute_effect_size basic" begin
        group_a = [make_st(2.0) for _ in 1:5]
        group_b = [make_st(0.0) for _ in 1:5]

        result = compute_effect_size(group_a, group_b)

        @test size(result.cohens_d) == (2, 3)
        @test size(result.mean_a) == (2, 3)
        @test size(result.mean_b) == (2, 3)
        # Large effect size expected
        @test all(abs.(result.cohens_d) .> 1.0)
        # means should be close to the input values
        @test all(abs.(result.mean_a .- 2.0) .< 0.1)
        @test all(abs.(result.mean_b) .< 0.1)
    end

    @testset "compute_effect_size identical groups" begin
        group = [make_st(1.0) for _ in 1:10]

        result = compute_effect_size(group, group)
        # Identical groups → near-zero effect size
        @test all(abs.(result.cohens_d) .< 0.5)
    end

    @testset "compare_windows ttest" begin
        # Create a result with distinct baseline and test windows
        n_ch, n_win, n_k = 2, 20, 3
        coeffs = zeros(n_ch, n_win, n_k)
        coeffs[:, 1:10, :] .= 0.0   # baseline = 0
        coeffs[:, 11:20, :] .= 2.0  # test = 2
        coeffs .+= 0.01 * randn(n_ch, n_win, n_k)
        errs = rand(n_ch, n_win)
        r = STResult(coeffs, errs, collect(1:n_win), collect(2:n_win+1),
                     ["a", "b"], Dict{String,Any}())

        comp = compare_windows(r, 1:10, 11:20; method="ttest")

        @test size(comp.p_value) == (2, 3)
        @test size(comp.statistic) == (2, 3)
        @test comp.method == "ttest"
        # Should detect the difference
        @test all(comp.p_value .< 0.05)
        # Baseline mean ≈ 0, test mean ≈ 2
        @test all(abs.(comp.baseline_mean) .< 0.5)
        @test all(abs.(comp.test_mean .- 2.0) .< 0.5)
    end

    @testset "compare_windows ranksum" begin
        n_ch, n_win, n_k = 2, 20, 3
        coeffs = zeros(n_ch, n_win, n_k)
        coeffs[:, 1:10, :] .= 0.0
        coeffs[:, 11:20, :] .= 3.0
        errs = rand(n_ch, n_win)
        r = STResult(coeffs, errs, collect(1:n_win), collect(2:n_win+1),
                     ["a", "b"], Dict{String,Any}())

        comp = compare_windows(r, 1:10, 11:20; method="ranksum")

        @test comp.method == "ranksum"
        @test all(comp.p_value .< 0.05)
    end

    @testset "compare_windows with CTResult" begin
        n_p, n_win, n_k = 3, 20, 3
        coeffs = zeros(n_p, n_win, n_k)
        coeffs[:, 1:10, :] .= 1.0
        coeffs[:, 11:20, :] .= 5.0
        errs = rand(n_p, n_win)
        r = CTResult(coeffs, errs, collect(1:n_win), collect(2:n_win+1),
                     ["a-b", "a-c", "b-c"], Dict{String,Any}())

        comp = compare_windows(r, 1:10, 11:20)
        @test size(comp.p_value) == (3, 3)
    end

    @testset "compare_windows invalid method" begin
        r = STResult(randn(1, 10, 1), rand(1, 10), collect(1:10), collect(2:11),
                     ["ch1"], Dict{String,Any}())
        @test_throws ErrorException compare_windows(r, 1:5, 6:10; method="invalid")
    end
end
