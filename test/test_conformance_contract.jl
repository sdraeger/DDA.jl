using Test
using JSON3
using DelayDifferentialAnalysis

@testset "Cross-Language Conformance Contract" begin
    dda_defaults = DelayDifferentialAnalysis.DDADefaults

    contract_path = normpath(joinpath(@__DIR__, "..", "conformance", "dda_conformance_contract.json"))
    contract = JSON3.read(read(contract_path, String))

    defaults = contract.defaults
    @test dda_defaults.WINDOW_LENGTH == Int(defaults.window_length)
    @test dda_defaults.WINDOW_STEP == Int(defaults.window_step)
    @test dda_defaults.DERIVATIVE_POINTS == Int(defaults.derivative_points)
    @test dda_defaults.POLYNOMIAL_ORDER == Int(defaults.polynomial_order)
    @test dda_defaults.NUM_TAU == Int(defaults.num_tau)
    @test collect(dda_defaults.MODEL_PARAMS) == Int.(defaults.model_terms)
    @test collect(dda_defaults.DELAYS) == Int.(defaults.delays)

    @test collect(VARIANT_ORDER) == String.(contract.variant_order)
    @test [v.abbreviation for v in active_variants()] == String.(contract.active_variants)

    for case in contract.select_mask_cases
        got = generate_select_mask(String.(case.variants))
        @test got == Int.(case.mask)
    end

    for abbrev in String.(contract.ct_window_required_for)
        variant = get_variant_by_abbrev(abbrev)
        @test variant !== nothing
        @test "-WL_CT" in variant.required_params
        @test "-WS_CT" in variant.required_params
    end
end
