function native_parity_fixture()
    samples = Matrix{Float64}(undef, 96, 3)
    for row in axes(samples, 1)
        t = row - 1
        x = sin(t * 0.11) + 0.03cos(t * 0.37)
        y = 0.42sin(t * 0.11 + 0.4) + 0.18cos(t * 0.07)
        z = 0.25x + 0.35y + 0.08sin(t * 0.19 + 0.2)
        samples[row, :] = [x, y, z]
    end
    return samples
end

function run_native_parity_fixture(device)
    return run_dda_matrix(
        native_parity_fixture();
        device=device,
        channels=[1, 2, 3],
        flavors=["ST", "CT", "CD", "DE", "SY"],
        window_length=32,
        window_step=16,
        delays=[1, 2],
        model_terms=[1, 2, 4],
        derivative_points=3,
        order=3,
        nr_tau=2,
        ct_channel_pairs=[(1, 2), (2, 3)],
        cd_channel_pairs=[(1, 2), (2, 1), (3, 2)],
        ct_window_length=2,
        ct_window_step=1,
    )
end

const NATIVE_EXPECTED = Dict(
    "ST" => [
        1.973444341356 1.841545857354 2.003185043142 1.885291325996
        2.025111106589 1.847634039673 1.986256200640 1.837801225124
        1.902336013334 1.777259029120 2.161582204063 1.899445940956
    ],
    "CT" => [
        1.990182542284 1.845707597562 1.991823686485 1.860031840421
        1.950545902080 1.821266554039 2.073120721445 1.865684403341
    ],
    "CD" => [
        0.005953420797 0.005251563906 0.011315953838 0.010937784290
        0.001965274791 0.012051424377 0.001590852203 0.015272098482
        0.002282108860 0.018393432760 0.003182694468 0.010282281326
    ],
    "DE" => [
        0.484785019689 0.054670849622 0.280102825445 0.016186927079
        0.606556023271 0.081611340763 0.612910698557 0.100614186601
    ],
    "SY" => [0.001195397813 -0.001673205084 0.001008255784 0.030790748346],
)

@testset "Native Julia DDA" begin
    result = run_native_parity_fixture("cpu")

    @test result.window_markers == [40.0, 56.0, 72.0, 88.0]
    @test [flavor.id for flavor in result.flavors] == ["ST", "CT", "CD", "DE", "SY"]
    @test result["ST"].row_labels == ["Ch 1", "Ch 2", "Ch 3"]
    @test result["CT"].row_labels == ["Ch 1-Ch 2", "Ch 2-Ch 3"]
    @test result["CD"].row_labels == ["Ch 1 <- Ch 2", "Ch 2 <- Ch 1", "Ch 3 <- Ch 2"]
    @test result["DE"].row_labels == ["Ch 1-Ch 2", "Ch 2-Ch 3"]
    @test result["SY"].row_labels == ["Ch 1 <-> Ch 2"]
    for (flavor, expected) in NATIVE_EXPECTED
        @test result[flavor].matrix ≈ expected rtol=1e-9 atol=1e-9
    end

    @test_throws ErrorException run_native_parity_fixture("metal")

    cuda = try
        Base.require(Main, :CUDA)
    catch
        nothing
    end
    if cuda !== nothing && cuda.functional()
        gpu_result = run_native_parity_fixture("cuda:0")
        for flavor in keys(NATIVE_EXPECTED)
            @test gpu_result[flavor].matrix ≈ result[flavor].matrix rtol=1e-8 atol=1e-10
        end
    else
        @info "Skipping native CUDA parity test: CUDA.jl or a functional GPU is unavailable"
    end
end
