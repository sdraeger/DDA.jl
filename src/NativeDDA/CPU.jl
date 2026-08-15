function _solve_cpu(problem::RegressionProblem)::SolvedBlock
    isempty(problem.design) && return _nan_block(size(problem.design, 2))
    factorization = try
        svd(problem.design; full=false)
    catch
        return _nan_block(size(problem.design, 2))
    end
    sigma_max = isempty(factorization.S) ? 0.0 : maximum(factorization.S)
    tolerance = max(size(problem.design)...) * eps(Float64) * max(sigma_max, 1.0)
    projected = factorization.U' * problem.fit_target
    scaled = zeros(Float64, length(projected))
    keep = factorization.S .> tolerance
    scaled[keep] = projected[keep] ./ factorization.S[keep]
    coefficients = factorization.V * scaled
    residual = problem.residual_target - problem.design * coefficients
    rmse = sqrt(sum(abs2, residual) / size(problem.design, 1))
    return SolvedBlock(coefficients, rmse)
end

_solve_cpu(problems::Vector{RegressionProblem}) = SolvedBlock[_solve_cpu(problem) for problem in problems]
