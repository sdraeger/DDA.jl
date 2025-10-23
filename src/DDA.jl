module DDA

using Dates
using UUIDs

include("error.jl")
include("types.jl")
include("parser.jl")
include("runner.jl")

export DDARunner, run, binary_path

export DDARequest, DDAResult, VariantResult
export TimeRange, PreprocessingOptions, AlgorithmSelection
export WindowParameters, ScaleParameters

export DDAError, BinaryNotFoundError, FileNotFoundError, UnsupportedFileTypeError
export ExecutionFailedError, ParseError, InvalidParameterError

end
