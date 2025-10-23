abstract type DDAError <: Exception end

struct BinaryNotFoundError <: DDAError
    path::String
end

struct FileNotFoundError <: DDAError
    path::String
end

struct UnsupportedFileTypeError <: DDAError
    file_type::String
end

struct ExecutionFailedError <: DDAError
    message::String
end

struct ParseError <: DDAError
    message::String
end

struct InvalidParameterError <: DDAError
    message::String
end

Base.showerror(io::IO, e::BinaryNotFoundError) = print(io, "DDA binary not found at: $(e.path)")
Base.showerror(io::IO, e::FileNotFoundError) = print(io, "Input file not found: $(e.path)")
Base.showerror(io::IO, e::UnsupportedFileTypeError) = print(io, "Unsupported file type: $(e.file_type)")
Base.showerror(io::IO, e::ExecutionFailedError) = print(io, "DDA execution failed: $(e.message)")
Base.showerror(io::IO, e::ParseError) = print(io, "Failed to parse DDA output: $(e.message)")
Base.showerror(io::IO, e::InvalidParameterError) = print(io, "Invalid parameter: $(e.message)")
