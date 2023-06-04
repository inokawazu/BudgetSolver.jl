using JuMP, Ipopt

raw"""
Budgeting scenario where the fixed costs and (excepted) income are known in advance.

Parameters:
`ntimesteps`
`starting_budget`
`income ≡ I[1..ntimesteps]`
`fixed_costs ≡ I[1..ntimesteps]`
"""
struct FixedICC{T <: Real}
    ntimesteps::Int
    starting_budget::T
    income::Vector{T}
    fixed_costs::Vector{T}

    function FixedICC{T}(starting_budget::T, income::Vector{T}, fixed_costs::Vector{T}) where T <: Real
        ntimesteps = length(income)
        ntimesteps > 1 || throw(ArgumentError("ntimesteps must be greater than 1"))
        
        length(income) == length(fixed_costs) || throw(ArgumentError("length(income) must be equal to length(fixed_costs)"))

        all(>=(0), income) || throw(ArgumentError("income must be greater than or equal to 0"))
        all(>=(0), fixed_costs) || throw(ArgumentError("fixed_costs must be greater than or equal to 0"))

        new{T}(ntimesteps, starting_budget, income, fixed_costs)
    end
end

function FixedICC(starting_budget::B, income::Vector{U}, fixed_costs::Vector{V}) where {B <: Real, U <: Real, V <: Real}
    FixedICC{promote_type(B, U, V)}(starting_budget, income, fixed_costs)
end

ntimesteps(ficc::FixedICC) = ficc.ntimesteps
starting_budget(ficc::FixedICC) = ficc.starting_budget
income(ficc::FixedICC) = ficc.income
fixed_costs(ficc::FixedICC) = ficc.fixed_costs

struct BudgetSolverResult{T <: Real}
    costs::Vector{T}
    wealth::Vector{T}
end

function ipopt_solve(
        bs, ideal_costs::AbstractVector; verbosity::Int=0, inflation_rate::Real=0.0, interest_rate::Real=0.0)
    m = Model()
    set_optimizer(m, Ipopt.Optimizer)
    set_attribute(m, "print_level", verbosity)

    local t, B, C
    @variable(m, 0 <= t[0:ntimesteps(bs)])
    @variable(m, 0 <= B[t in 0:ntimesteps(bs)])
    @variable(m, 0 <= C[t in 0:ntimesteps(bs)])

    @constraint(m, B[0] == starting_budget(bs))
    @constraint(m, C[0] == 0)

    I = income(bs)
    Cf = fixed_costs(bs)
    rho = inflation_rate
    rha = interest_rate

    @constraint(m, [t in 1:ntimesteps(bs)], B[t] == (1 - rho)*(1 + rha)*B[t-1] + I[t] - C[t])
    @constraint(m, [t in 1:ntimesteps(bs)], C[t] >= Cf[t])
    
    local Max
    @objective(m, Max, -sum((ideal_costs[t] - C[t])^2 for t in 1:ntimesteps(bs)))

    optimize!(m)

    return BudgetSolverResult([value(C[i]) for i in 1:ntimesteps(bs)],
                       [value(B[i]) for i in 1:ntimesteps(bs)])
end
