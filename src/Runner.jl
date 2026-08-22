#=
DDA Binary Runner Module

Provides functions to execute DDA analysis by running the
run_DDA_AsciiEdf binary and parsing results.
=#

module Runner

import UUIDs
import Dates
using ..Flavors: SELECT_MASK_SIZE, DEFAULT_DELAYS, FileType, VariantMetadata,
    get_variant_by_abbrev, generate_select_mask, parse_select_mask,
    requires_ct_params, Individual, Pairs, DirectedPairs,
    require_binary, REQUIRES_SHELL_WRAPPER, SHELL_COMMAND,
    EDF, ASCII, get_flag, file_type_from_extension
using ..DDADefaults
using ..ModelEncoding: model_matrix_to_encoding

include("Runner/Types.jl")
include("Runner/Request.jl")
include("Runner/Labels.jl")
include("Runner/Results.jl")
include("Runner/Command.jl")
include("Runner/Parsing.jl")
include("Runner/Execution.jl")

end # module Runner
