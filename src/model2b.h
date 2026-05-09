#include "global_variables.h"
#include "boost_1.h"
#include "global_functions.h"

double dummy(instance* inst);

void master2b(instance* inst);

LP_return solveCG_LP(instance* inst, vector<int>& road, vector<int>& next_to_avoid);

void master2c(instance* inst);
