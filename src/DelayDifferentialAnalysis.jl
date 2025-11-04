module DelayDifferentialAnalysis

using Dates
using UUIDs

include("error.jl")
include("types.jl")
include("parser.jl")
include("runner.jl")

export DDARunner, run_dda, binary_path

export DDARequest, DDAResult, VariantResult
export Bounds, AlgorithmSelection
export WindowParameters, DelayParameters

export DDAError, BinaryNotFoundError, FileNotFoundError, UnsupportedFileTypeError
export ExecutionFailedError, ParseError, InvalidParameterError

end
