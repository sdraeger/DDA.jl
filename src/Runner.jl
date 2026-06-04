#=
DDA Binary Runner Module

Provides functions to execute DDA analysis by running the
run_DDA_AsciiEdf binary and parsing results.
=#

module Runner

import UUIDs
import Dates
using ..Variants
using ..DDADefaults
using ..ModelEncoding: model_matrix_to_encoding

export DDARequest, DDAResult, VariantResultData
export DDARunner, run_DDA
export StructuredTimepoint, StructuredChannelData
export run_analysis_structured, parse_output_file_structured

include("Runner/Types.jl")
include("Runner/Request.jl")
include("Runner/Labels.jl")
include("Runner/Results.jl")
include("Runner/Command.jl")
include("Runner/Parsing.jl")
include("Runner/Execution.jl")

end # module Runner
