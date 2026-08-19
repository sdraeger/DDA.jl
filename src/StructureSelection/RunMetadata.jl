const _RUN_METADATA_FILENAME = "structure_selection.toml"
const _RUN_METADATA_VERSION = 1

function _write_structure_selection_run(run::StructureSelectionRun)::String
    metadata = Dict{String,Any}(
        "format_version" => _RUN_METADATA_VERSION,
        "MOD" => [collect(row) for row in eachrow(run.MOD)],
        "DDAorder" => run.DDAorder,
        "nr_delays" => run.nr_delays,
        "delays" => run.delays,
        "model_numbers" => run.model_numbers,
        "derivative_points" => run.derivative_points,
        "tau_file_suffix" => run.tau_file_suffix,
        "trial_prefix" => run.trial_prefix,
    )
    run.channels === nothing || (metadata["channels"] = run.channels)

    buffer = IOBuffer()
    TOML.print(buffer, metadata; sorted=true)
    contents = String(take!(buffer))
    path = joinpath(run.prefix, _RUN_METADATA_FILENAME)
    isfile(path) && read(path, String) == contents && return path

    temporary, io = mktemp(run.prefix)
    try
        write(io, contents)
        close(io)
        mv(temporary, path; force=true)
    catch
        isopen(io) && close(io)
        rm(temporary; force=true)
        rethrow()
    end
    return path
end

function _structure_selection_read(prefix)::StructureSelectionRun
    output_root = expanduser(String(prefix))
    path = joinpath(output_root, _RUN_METADATA_FILENAME)
    isfile(path) || error(
        "No structure-selection metadata found at $path; run `structure_selection_compute` first",
    )
    metadata = TOML.parsefile(path)
    get(metadata, "format_version", nothing) == _RUN_METADATA_VERSION ||
        error("Unsupported structure-selection metadata format in $path")

    MOD = _read_MOD_rows(metadata["MOD"], path)
    run = StructureSelectionRun(
        output_root,
        MOD,
        Int(metadata["DDAorder"]),
        Int(metadata["nr_delays"]),
        Int[metadata["delays"]...],
        haskey(metadata, "channels") ? Int[metadata["channels"]...] : nothing,
        _resolve_MOD_numbers(MOD, metadata["model_numbers"]),
        Int(metadata["derivative_points"]),
        String(metadata["tau_file_suffix"]),
        String(metadata["trial_prefix"]),
    )
    _LAST_STRUCTURE_SELECTION_RUN[] = run
    return run
end

function _read_MOD_rows(rows, path::AbstractString)::Matrix{Int}
    rows isa AbstractVector && !isempty(rows) || error("Invalid MOD data in $path")
    column_count = length(first(rows))
    column_count > 0 || error("Invalid MOD data in $path")
    MOD = Matrix{Int}(undef, length(rows), column_count)
    for (row_index, row) in enumerate(rows)
        length(row) == column_count || error("Inconsistent MOD row lengths in $path")
        MOD[row_index, :] = Int[row...]
    end
    _validate_binary_MOD(MOD)
    return MOD
end
