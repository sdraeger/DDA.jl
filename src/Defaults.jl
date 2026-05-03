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

"""CLI flag constants for the DDA binary."""
module DDAFlags
    const DATA_FILE = "-DATA_FN"
    const OUTPUT_FILE = "-OUT_FN"
    const CHANNEL_LIST = "-CH_list"
    const SELECT_MASK = "-SELECT"
    const MODEL = "-MODEL"
    const DELAY_VALUES = "-TAU"
    const DERIVATIVE_POINTS = "-dm"
    const POLYNOMIAL_ORDER = "-order"
    const NUM_TAU = "-nr_tau"
    const WL = "-WL"
    const WS = "-WS"
    const WINDOW_LENGTH = WL
    const WINDOW_STEP = WS
    const CT_WINDOW_LENGTH = "-WL_CT"
    const CT_WINDOW_STEP = "-WS_CT"
    const TIME_BOUNDS = "-StartEnd"
    const SAMPLING_RATE = "-SR"
    const EDF_FLAG = "-EDF"
    const ASCII_FLAG = "-ASCII"
end
