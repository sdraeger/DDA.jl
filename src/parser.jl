function parse_dda_output(content::String)::Matrix{Float64}
    lines = split(content, '\n')
    matrix = Vector{Vector{Float64}}()

    for line in lines
        line_stripped = strip(line)

        if isempty(line_stripped) || startswith(line_stripped, '#')
            continue
        end

        values = Float64[]
        for token in split(line_stripped)
            try
                val = parse(Float64, token)
                if isfinite(val)
                    push!(values, val)
                end
            catch
                continue
            end
        end

        if !isempty(values)
            push!(matrix, values)
        end
    end

    if isempty(matrix)
        throw(ParseError("No valid data found in DDA output"))
    end

    matrix_array = reduce(hcat, matrix)'

    @info "Loaded DDA output shape: $(size(matrix_array, 1)) rows × $(size(matrix_array, 2)) columns"

    if size(matrix_array, 1) > 0 && size(matrix_array, 2) >= 10
        @debug "First row sample (first 10 values): $(matrix_array[1, 1:10])"
    end

    if size(matrix_array, 2) > 2
        after_skip = matrix_array[:, 3:end]

        @debug "After skipping first 2 columns: $(size(after_skip, 1)) rows × $(size(after_skip, 2)) columns"

        if size(after_skip, 1) > 0 && size(after_skip, 2) >= 10
            @debug "After skip, first row (first 10 values): $(after_skip[1, 1:10])"
        end

        col_indices = 1:4:size(after_skip, 2)
        extracted = after_skip[:, col_indices]

        if size(extracted, 1) > 0 && size(extracted, 2) >= 5
            @debug "First extracted row sample (first 5 values): $(extracted[1, 1:min(5, size(extracted, 2))])"
        end

        if isempty(extracted) || size(extracted, 2) == 0
            throw(ParseError("No data after column extraction"))
        end

        @info "Extracted matrix shape: $(size(extracted, 1)) rows × $(size(extracted, 2)) columns (time windows × delays/scales)"

        transposed = extracted'

        if isempty(transposed) || size(transposed, 2) == 0
            throw(ParseError("Transpose resulted in empty data"))
        end

        @info "Transposed to: $(size(transposed, 1)) channels × $(size(transposed, 2)) timepoints"

        return transposed
    else
        return reshape(vec(matrix_array), 1, :)
    end
end
