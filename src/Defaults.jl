"""Default parameter values for DDA analysis."""
module DDADefaults
    const MODEL_PARAMS = [1, 2, 10]
    const MODEL_DIMENSION = 4
    const POLYNOMIAL_ORDER = 4
    const NUM_TAU = 2
    const WINDOW_LENGTH = 200
    const WINDOW_STEP = 100
    const DELAYS = (7, 10)
    const SAMPLING_RATE = (500, 1000)
end

"""CLI flag constants for the DDA binary."""
module DDAFlags
    const DATA_FILE = "-DATA_FN"
    const OUTPUT_FILE = "-OUT_FN"
    const CHANNEL_LIST = "-CH_list"
    const SELECT_MASK = "-SELECT"
    const MODEL = "-MODEL"
    const DELAY_VALUES = "-TAU"
    const MODEL_DIMENSION = "-dm"
    const POLYNOMIAL_ORDER = "-order"
    const NUM_TAU = "-nr_tau"
    const WINDOW_LENGTH = "-WL"
    const WINDOW_STEP = "-WS"
    const CT_WINDOW_LENGTH = "-WL_CT"
    const CT_WINDOW_STEP = "-WS_CT"
    const TIME_BOUNDS = "-StartEnd"
    const SAMPLING_RATE = "-SR"
    const EDF_FLAG = "-EDF"
    const ASCII_FLAG = "-ASCII"
end
