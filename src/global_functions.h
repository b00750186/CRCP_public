#include "global_variables.h"

void readData(instance* inst);

void generate_latex_model2(instance* inst, IloNumArray solution_p);

void generate_latex_model1(instance* inst, vector<vector<double>> solution);

void write_log(instance* inst, int method, double objective, double bestlp, double solution_time, int status);

void write_log_comp(instance* inst, int method, double objective, double bestlp, double rootlb, double solution_time, int status, long long nodes);

void write_log1(instance* inst, int method, double objective, double solution_time, int numberofcycles, int numberoflambda, int status);

void write_log21(instance* inst, int method, double objective, double bestMIP, double LB, double solution_time, int numberofcycles, int numberoflambda, int status);

int arc_coef(int i, int j, instance* inst);

int findNearestNode(instance* inst, int currentNode, const vector<bool>& visited);









