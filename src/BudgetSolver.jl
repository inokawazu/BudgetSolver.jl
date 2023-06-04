module BudgetSolver

include("fixedicc.jl")

struct PlannedIncomeScenario
    income
    starting_budget
    budget_interest
end

times(pis::PlannedIncomeScenario) = eachindex(pis.income)
starting_time(pis::PlannedIncomeScenario) = firstindex(pis.income)
starting_cost(::PlannedIncomeScenario) = 0

budget_interest((; budget_interest)::PlannedIncomeScenario) = budget_interest
income((; income)::PlannedIncomeScenario) = income
starting_budget((; starting_budget)::PlannedIncomeScenario) = starting_budget

using JuMP

function budget_model(scenario; vb = 0.01)
    t = times(scenario)

    m = Model()
    # set_optimizer(m, Ipopt.Optimizer)
    # set_attribute(m, "print_level", verbosity)

    local B, C
    @variable(m, 0 <= B[t=t])
    @variable(m, 0 <= C[t=t])

    I = income(scenario)

    t0 = starting_time(scenario)
    @constraint(m, B[t0] == starting_budget(scenario))
    @constraint(m, C[t0] == starting_cost(scenario))
    
    
    rho = budget_interest(scenario)
    for (ti, tip1) in zip(t, Iterators.drop(t, 1))
        @constraint(m, B[tip1] == (1+rho)*B[ti] + I[ti] - C[ti])
    end

    @NLobjective(m, MAX_SENSE, sum(log(1 + c) for c in C) + vb * sum(log(1 + b) for b in B))

    return m
end

end # module BudgetSolver
