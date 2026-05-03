using Test
using DelayDifferentialAnalysis

@testset "Results" begin
    @testset "STResult construction and accessors" begin
        n_ch, n_win, n_k = 3, 5, 3
        coeffs = randn(n_ch, n_win, n_k)
        errs = rand(n_ch, n_win)
        ws = collect(1:100:500)
        we = ws .+ 99
        labels = ["ch1", "ch2", "ch3"]
        params = Dict{String,Any}("WL" => 100)

        r = STResult(coeffs, errs, ws, we, labels, params)

        @test n_channels(r) == 3
        @test n_windows(r) == 5
        @test n_coeffs(r) == 3
        @test r.channel_labels == labels
        @test r.params["WL"] == 100
        @test size(r.coefficients) == (3, 5, 3)
        @test size(r.errors) == (3, 5)
        @test r.T == Float64.(ws)
        @test r.t == Float64.(ws)
    end

    @testset "CTResult construction and accessors" begin
        n_p, n_win, n_k = 6, 10, 3
        coeffs = randn(n_p, n_win, n_k)
        errs = rand(n_p, n_win)
        ws = collect(1:50:500)
        we = ws .+ 49
        pair_labels = ["ch1-ch2", "ch1-ch3", "ch1-ch4", "ch2-ch3", "ch2-ch4", "ch3-ch4"]
        params = Dict{String,Any}()

        r = CTResult(coeffs, errs, ws, we, pair_labels, params)

        @test n_pairs(r) == 6
        @test n_windows(r) == 10
        @test n_coeffs(r) == 3
        @test r.pair_labels == pair_labels
        @test r.T == Float64.(ws)
        @test r.t == Float64.(ws)
    end

    @testset "DEResult construction and accessors" begin
        n_win = 20
        erg = randn(n_win)
        ws = collect(1:100:2000)
        we = ws .+ 99
        params = Dict{String,Any}()

        r = DEResult(erg, ws, we, params)

        @test n_windows(r) == 20
        @test r.ergodicity == erg
        @test length(r.window_starts) == 20
        @test r.T == Float64.(ws)
        @test r.t == Float64.(ws)
    end

    @testset "STResult single-channel edge case" begin
        coeffs = randn(1, 1, 3)
        errs = rand(1, 1)
        r = STResult(coeffs, errs, [1], [100], ["only"], Dict{String,Any}())

        @test n_channels(r) == 1
        @test n_windows(r) == 1
        @test n_coeffs(r) == 3
    end
end
