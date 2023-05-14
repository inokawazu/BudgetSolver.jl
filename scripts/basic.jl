using BudgetSolver

# Parameters
nmonths = 12*40 # 40 years
starting_budget = 100.0
lifestyle_interest = 0.01/12
overall_interest = 0.02/12

income = [starting_budget for i in 1:nmonths]
fixed_costs = income ./ 3

ideal_costs = [90.0 * (1 + lifestyle_interest) ^ i for i in 1:nmonths]

scenario = BudgetSolver.FixedICC(starting_budget, income, fixed_costs)

# println(scenario)

solution = BudgetSolver.ipopt_solve(scenario, ideal_costs, inflation_rate = overall_interest)
