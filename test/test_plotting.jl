using Test
using DelayDifferentialAnalysis

@testset "Plotting" begin
    # Plotting tests are smoke tests that only run if Plots.jl is available
    plots_available = try
        DelayDifferentialAnalysis.Plotting._ensure_plots()
        true
    catch
        false
    end

    if !plots_available
        @info "Skipping plotting tests: Plots.jl not available"
    else
        @testset "plot_coefficients" begin
            r = STResult(
                randn(2, 10, 3), rand(2, 10),
                collect(1:10), collect(2:11),
                ["ch1", "ch2"], Dict{String,Any}()
            )
            p = plot_coefficients(r)
            @test p !== nothing
        end

        @testset "plot_coefficients with time" begin
            r = STResult(
                randn(2, 10, 3), rand(2, 10),
                collect(0:100:900), collect(99:100:999),
                ["ch1", "ch2"], Dict{String,Any}()
            )
            p = plot_coefficients(r; use_time=true, sfreq=256.0)
            @test p !== nothing
        end

        @testset "plot_errors" begin
            r = STResult(
                randn(2, 10, 3), rand(2, 10),
                collect(1:10), collect(2:11),
                ["ch1", "ch2"], Dict{String,Any}()
            )
            p = plot_errors(r)
            @test p !== nothing
        end

        @testset "plot_heatmap" begin
            r = STResult(
                randn(3, 15, 3), rand(3, 15),
                collect(1:15), collect(2:16),
                ["ch1", "ch2", "ch3"], Dict{String,Any}()
            )
            p = plot_heatmap(r; coeff_index=1)
            @test p !== nothing
        end

        @testset "plot_ergodicity" begin
            r = DEResult(
                randn(10), collect(1:10), collect(2:11), Dict{String,Any}()
            )
            p = plot_ergodicity(r)
            @test p !== nothing
        end

        @testset "plot_model" begin
            p = plot_model([1, 2, 10])
            @test p !== nothing
        end
    end
end
