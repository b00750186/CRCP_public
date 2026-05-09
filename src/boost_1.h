#include "global_variables.h"

//subproblem_return boost_shortest11(instance* inst);

//subproblem_return tboost_shortest11(instance* inst);

///*subproblem_return*/ void boost_shortest11(instance* inst, vector<int> road);

void precompute_min_f2c_via_f(const instance* inst);

subproblem_return_single boost_shortest(int customer, instance* inst, IloNumArray dualCapacity, IloArray<NumMatrix> dualArcs, IloNumArray dualCoverage);

void boost_shortest11(instance* inst, vector<int>& road, const std::vector<char>& inRoad, const std::vector<char>& avoid);

//void boost_shortest11(instance* inst, vector<int> &road, vector<int> &next_to_avoid);
