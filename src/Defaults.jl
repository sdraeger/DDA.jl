"""Default parameter values for DDA analysis."""
module DDADefaults
    const MODEL_PARAMS = [1, 2, 10]
    const DERIVATIVE_POINTS = 3
    const POLYNOMIAL_ORDER = 4
    const NUM_TAU = 2
    const WL = nothing
    const WS = nothing
    const WINDOW_LENGTH = WL
    const WINDOW_STEP = WS
    const DELAYS = (7, 10)
    const SAMPLING_RATE = nothing
end
