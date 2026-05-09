#include "model2b.h"

int calculateModifiedHammingDistance(const std::vector<int>& vec1, const std::vector<int>& vec2) {
	int distance = 0;
	size_t minLength = std::min(vec1.size(), vec2.size());

	// Calculate distance for the common length part
	for (size_t i = 0; i < minLength; ++i) {
		if (vec1[i] != vec2[i]) {
			++distance;
		}
	}

	// Add the difference in lengths to the distance
	distance += std::max((int)vec1.size(), (int)vec2.size()) - (int)minLength;

	return distance;
}

vector<vector<int>> generatetours(vector<int> road_to_take, vector<int> next_to_avoid, int size, int number, instance* inst) {
	vector<vector<int>> result;

	int size_return;
	if (size > inst->dimension_f)
		size_return = inst->dimension_f;
	else
		size_return = size;


	std::vector<int> v;

	for (int i = 1; i < inst->dimension_f; i++) {
		if (find(road_to_take.begin(), road_to_take.end(), i) == road_to_take.end())
			v.push_back(i);
	}


	// Seed the random number generator
	std::srand(static_cast<unsigned int>(std::time(0)));

	for (int i = 0; i < number; i++) {


		// Shuffle the vector
		while (true) {
			std::random_shuffle(v.begin(), v.end());
			if (find(next_to_avoid.begin(), next_to_avoid.end(), v[0]) == next_to_avoid.end())
				break;
		}

		vector<int> v0 = road_to_take;
		v0.insert(v0.end(), v.begin(), v.begin() + max(int(size_return - road_to_take.size()), 0));

		double capacity = 0.0;
		for (auto vind : v0) capacity += inst->capacity[vind];
		if (capacity > inst->totaldemand) {
			result.push_back(v0);
		}
	}

	return result;
}



vector<path> generateInitialPaths11_withconstr_0(instance* inst, vector<int> road_to_take, vector<int> next_to_avoid);
vector<path> generateInitialPaths11_withconstr_mip(instance* inst, vector<int> road_to_take, vector<int> next_to_avoid);

void print_path11(int i, instance* inst) {

	vector<int> visitedfacilities = inst->paths[i].getFacilityPath();
	for (int j = 0; j < visitedfacilities.size(); j++) {

		std:: cout << visitedfacilities[j] << " ";

	}
	std::cout << endl;
}

static ofstream time_profile;
static map<std::vector<int>, double> tours_solved;
static int tour_repetition_counter = 0;

//solves with fixed tour
double tour_solver1(instance* inst, vector<int> tour,bool print, double cutoff = -1.0) {

	IloEnv env = IloEnv();
	IloModel model1 = IloModel(env);

	IloNumVar::Type inttype = IloNumVar::Int;
	
	//IloNumVar::Type inttype = IloNumVar::Float;


	//Objective
	IloObjective  Lateness = IloAdd(model1, IloMinimize(env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(env, inst->dimension_c, 0, 1, inttype);
		/*
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());


		}
		*/
	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(env, inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			/*
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());

			}*/
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		//string name = "s." + to_string(k);
		//s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			//string name = "arc." + to_string(i) + "." + to_string(j);

			arcs[arc_coef(i, j, inst)] = IloNumVar(env, 0.0, 0.0, inttype);

			/*
			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(env, 0.0, 0.0, inttype);
			}
			else
				arcs[arc_coef(i, j, inst)] = IloNumVar(env, 0.0, 1.0, inttype);
			*/
			//arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}


	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {

		IloExpr v(env);
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += inst->c_f[i][ii] * f[k][i][ii];

			}

		for (IloInt i = 0; i < inst->dimension_f; i++) {

			v += z[i][k] * inst->c_fc[i][k];

		}

		model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}

			if (i == 0)
				model1.add(v == -1 + z[0][k]);
			else
				model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++)
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				model1.add(f[k][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(env);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				//v += arcs[ii][i];
				v += arcs[arc_coef(ii, i, inst)];

			}
			model1.add(z[i][k] <= v);
			v.end();

		}

	for (int k = 0; k < inst->dimension_c; k++) {
		IloExpr v(env);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			v += z[i][k];
		}
		model1.add(v == 1.0);
	}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		model1.add(v <= inst->capacity[i]);
		v.end();
	}

	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(env);
		IloExpr Xdeltaplus(env);
		for (int ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		model1.add(Xdeltaminus <= 1);
		model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}
	int prev = 0;
	for (auto v : tour) if (v != 0) {
		model1.add(arcs[arc_coef(prev, v, inst)] == 1);
		arcs[arc_coef(prev, v, inst)].setUb(1.0);
		prev = v;
	}
	model1.add(arcs[arc_coef(prev, 0, inst)] == 1);
	arcs[arc_coef(prev, 0, inst)].setUb(1.0);

	IloCplex pathSolver(model1);
	//string filename = inst->input_file + "solution";
	//pathSolver.readMIPStarts(filename.c_str());
	//pathSolver.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	pathSolver.setParam(IloCplex::Param::TimeLimit, 600);
	if (cutoff > 0) {
		pathSolver.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, cutoff);
	}
	//pathSolver.exportModel("finalMIP1.lp");
	bool feasible = pathSolver.solve();
	double solution;
	if (feasible) {
		solution = pathSolver.getObjValue();

		if (print) {
			cout << "ARCS Values: " << endl;
			for (int i = 0; i < inst->dimension_f; i++)
				for (int j = 0; j < inst->dimension_f; j++) {
					if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {
						inst->active_arcs.push_back(pairs(i, j));
						cout << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << ", ";
					}
				}

			vector<vector<double>> solutions = vector<vector<double>>(inst->dimension_f);
			cout << "Assignment Values: " << endl;
			for (int i = 0; i < inst->dimension_f; i++) {
				solutions[i] = vector<double>(inst->dimension_c);
				for (int k = 0; k < inst->dimension_c; k++) {
					if(pathSolver.getValue(z[i][k]) > RC_EPS)
						cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << ", ";
					solutions[i][k] = pathSolver.getValue(z[i][k]);
				}
			}

			generate_latex_model1(inst, solutions);		

		}
	}
	else
		solution = IloInfinity;



	pathSolver.clearModel();
	env.end();

	return solution;
}

double tour_forest_solver1(vector<int>* best_tour,instance* inst, bool print, double cutoff = -1.0) {

	IloEnv env = IloEnv();
	IloModel model1 = IloModel(env);

	IloNumVar::Type inttype = IloNumVar::Int;

	//IloNumVar::Type inttype = IloNumVar::Float;


	//Objective
	IloObjective  Lateness = IloAdd(model1, IloMinimize(env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(env, inst->dimension_c, 0, 1, inttype);
		/*
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());
		}
		*/
	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(env, inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			/*
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());

			}*/
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		//string name = "s." + to_string(k);
		//s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			//string name = "arc." + to_string(i) + "." + to_string(j);

			arcs[arc_coef(i, j, inst)] = IloNumVar(env, 0.0, 0.0, inttype);

			/*
			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(env, 0.0, 0.0, inttype);
			}
			else
				arcs[arc_coef(i, j, inst)] = IloNumVar(env, 0.0, 1.0, inttype);
			*/
			//arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}


	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {

		IloExpr v(env);
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += inst->c_f[i][ii] * f[k][i][ii];

			}

		for (IloInt i = 0; i < inst->dimension_f; i++) {

			v += z[i][k] * inst->c_fc[i][k];

		}

		model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}

			if (i == 0)
				model1.add(v == -1 + z[0][k]);
			else
				model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++)
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				model1.add(f[k][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(env);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				//v += arcs[ii][i];
				v += arcs[arc_coef(ii, i, inst)];

			}
			model1.add(z[i][k] <= v);
			v.end();

		}

	for (int k = 0; k < inst->dimension_c; k++) {
		IloExpr v(env);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			v += z[i][k];
		}
		model1.add(v == 1.0);
	}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		model1.add(v <= inst->capacity[i]);
		v.end();
	}

	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(env);
		IloExpr Xdeltaplus(env);
		for (int ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		model1.add(Xdeltaminus <= 1);
		model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}

	/*
	for (int i = 0; i < inst->dimension_f; ++i) {
		for (int j = 0; j < inst->dimension_f; ++j) {
			if (arcs_allowed[i][j] > RC_EPS) {
				arcs[arc_coef(i, j, inst)].setUb(1.0);
			}
		}
	}
	*/

	set<vector<int>> uniquepaths;
	for (auto p : inst->paths) {
		uniquepaths.insert(p.getFacilityPath());
	}

	for (auto path : uniquepaths) {
		int prev = 0;
		if(path.size()>1)
			for (int i = 1; i < path.size(); i++) {
				arcs[arc_coef(path[prev], path[i], inst)].setUb(1.0);
				prev = i;
			}
		arcs[arc_coef(path[prev], 0, inst)].setUb(1.0);
	}

	IloCplex pathSolver(model1);
	//string filename = inst->input_file + "solution";
	//pathSolver.readMIPStarts(filename.c_str());
	//pathSolver.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	pathSolver.setParam(IloCplex::Param::TimeLimit, 600);
	if (cutoff > 0) {
		pathSolver.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, cutoff);
	}
	//pathSolver.exportModel("finalMIP1.lp");
	bool feasible = pathSolver.solve();
	double solution;
	if (feasible) {
		solution = pathSolver.getObjValue();
		if (solution < inst->bestMIP) {
			inst->bestMIP = solution;
		}
		best_tour->clear();
		int previ = 0;
		best_tour->push_back(previ);
		set<int> nextcandidate;
		for (int i = 1; i < inst->dimension_f; i++)
			nextcandidate.insert(i);
		while (best_tour->size() < inst->dimension_f) {
			int next;
			double max = -INFINITY;
			for (auto i : nextcandidate) {
				if (pathSolver.getValue(arcs[arc_coef(previ, i, inst)])> max) {
					max = pathSolver.getValue(arcs[arc_coef(previ,i,inst)]);
					next = i;
				}
			}
			best_tour->push_back(next);
			previ = next;
			nextcandidate.erase(next);
		}



		if (print) {
			cout << "ARCS Values: " << endl;
			for (int i = 0; i < inst->dimension_f; i++)
				for (int j = 0; j < inst->dimension_f; j++) {
					if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {
						inst->active_arcs.push_back(pairs(i, j));
						cout << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << ", ";
					}
				}

			vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
			cout << "Assignment Values: " << endl;
			for (int i = 0; i < inst->dimension_f; i++) {
				solution[i] = vector<double>(inst->dimension_c);
				for (int k = 0; k < inst->dimension_c; k++) {
					if (pathSolver.getValue(z[i][k]) > RC_EPS)
						cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << ", ";
					solution[i][k] = pathSolver.getValue(z[i][k]);
				}
			}

			generate_latex_model1(inst, solution);

		}
	}
	else
		solution = IloInfinity;



	pathSolver.clearModel();
	env.end();
	inst->branchtree << "MIP solution with non-zero LP arcs allowed: " << solution << endl;
	return solution;
}

double local_search(instance* inst, vector<int> tour) {
	//first operation swappping with facility within radius not in tour
	if (tours_solved.find(tour) != tours_solved.end()) {
		tour_repetition_counter++;
		return tours_solved[tour];
	}
	vector<int> besttour = tour;
	double radius = 10.0;
	set<int> free;

	for (int i = 0; i < inst->dimension_f; i++) {
		free.insert(i);
	}

	for (int i : tour) {
		free.erase(i);
	}
	
	
	double currentCost = tour_solver1(inst,tour,false);
	tours_solved[tour] = currentCost;
	inst->branchtree << "LOCAL start: " << currentCost<<"|" <<endl;


	//DELETING NODES ? how many to delete
	if (false) {
		for (int counter = 0; counter < 5; counter++) {
			vector<int> tempv = tour;
			int todelete = -1;
			for (int i = 0; i < tour.size(); i++) {
				tempv = tour;
				tempv.erase(std::remove(tempv.begin(), tempv.end(), i), tempv.end());
				double newCost;
				if (tours_solved.find(tour) == tours_solved.end()) {
					newCost = tour_solver1(inst, tempv, false);
					tours_solved[tempv] = newCost;
				}
				else {
					tour_repetition_counter++;
					newCost = tours_solved[tour];
				}
				if (newCost < currentCost) {
					todelete = i;
					currentCost = newCost;
					besttour = tempv;
				}
			}
			tour.erase(std::remove(tour.begin(), tour.end(), todelete), tour.end());
			if (todelete > -1)
				free.insert(todelete);
		}

		inst->branchtree << "LOCAL deleting I: " << currentCost << "|" << endl;
	}
	//INSERTING

	if (false) {
		set<int> tempfree = free;
		bool have_time = true;
		for (int repeat = 0; repeat < 10 && have_time; repeat++) {
			int bestins=0;
			for(auto ins: tempfree) {

				for (int i = 1; i < tour.size(); i++) {
					vector<int> tempv = tour;
					//tempv = tour;
					tempv.insert(tempv.begin() + i, ins);
					double newCost;
					if (tours_solved.find(tempv) == tours_solved.end()) {
						newCost = tour_solver1(inst, tempv, false);
						tours_solved[tempv] = newCost;
					}
					else {
						tour_repetition_counter++;
						newCost = tours_solved[tempv];
					}
					if (newCost < currentCost) {
						currentCost = newCost;
						besttour = tempv;
						bestins = ins;
					}
				}

				tour = besttour;
				tempfree.erase(bestins);
				if ((double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC > inst->timelimit) {
					inst->branchtree << "Timelimit in local search." << endl;
					have_time = false;
					break;
				}

			}			
		}
		inst->branchtree << "LOCAL insert: " << currentCost << "|" << endl;
	}
	//Update free
	for (int i : tour) {
		free.erase(i);
	}

	if (false) {
		//EXTERNAL SWAPPING
		for (int i = 0; i < tour.size(); i++) {
			int lastentry = 0;
			for (int j : free) {

				if (inst->c_f[i][j] < radius) {
					std::swap(tour[i], j);

					// Calculate the new cost after the swap
					double newCost;
					if (tours_solved.find(tour) == tours_solved.end()) {
						newCost = tour_solver1(inst, tour, false);
						tours_solved[tour] = newCost;
					}
					else {
						tour_repetition_counter++;
						newCost = tours_solved[tour];
					}

					// If the new cost is lower, keep the swap; otherwise, revert the swap
					if (newCost < currentCost) {
						currentCost = newCost;
						besttour = tour;
						lastentry = j;
					}
					else {
						// Revert the swap
						std::swap(tour[i], j);
					}
				}
			}
			free.erase(lastentry);
		}

		inst->branchtree << "LOCAL external swap: " << currentCost << "|" << endl;
	}
	
	//INTERNAL INVERSING

	if (false) {
		for (int i = 0; i < tour.size() - 1; i++) {
			for (int j = i + 1; j < tour.size() && j < i + 4; j++) {

				std::reverse(tour.begin() + i, tour.begin() + j + 1);
				double newCost;
				if (tours_solved.find(tour) == tours_solved.end()) {
					newCost = tour_solver1(inst, tour, false, currentCost);
					tours_solved[tour] = newCost;
				}
				else {
					tour_repetition_counter++;
					newCost = tours_solved[tour];
				}

				// If the new cost is lower, keep the swap; otherwise, revert the swap
				if (newCost < currentCost - RC_EPS) {
					currentCost = newCost;
					besttour = tour;
				}
				else {
					//Revert
					std::reverse(tour.begin() + i, tour.begin() + j + 1);

				}

			}
			inst->branchtree << "LOCAL total invert intermed i " << i << " : " << currentCost << "|" << endl;
		}
		inst->branchtree << "LOCAL total invert: " << currentCost << "|" << endl;
	}
	bool improved = true;
		
	//INTERNAL SWAPPING
	improved = true;
	bool have_time = true;
	besttour = tour;
	for (int counter = 0; counter < 7 && improved && have_time; counter++) {
		improved = false;
		vector<int> currenttour = besttour;
		for (int i = 0; i < currenttour.size(); i++) {
			if ((double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC > inst->timelimit) {
				inst->branchtree << "Timelimit in local search." << endl;
				have_time = false;
				break;
			}
			for (int j = i + 1; j < currenttour.size(); j++) {
				// Swap elements i and j in the tour
				std::swap(currenttour[i], currenttour[j]);

				// Calculate the new cost after the swap
				double newCost;
				if (tours_solved.find(currenttour) == tours_solved.end()) {
					newCost = tour_solver1(inst, currenttour, false, currentCost);
					tours_solved[currenttour] = newCost;
				}
				else {
					tour_repetition_counter++;
					newCost = tours_solved[currenttour];
				}

				
				if (newCost < currentCost - RC_EPS) {
					currentCost = newCost;
					besttour = currenttour;
					improved = true;
				}else
					
					std::swap(currenttour[i], currenttour[j]);

			}
			
		}
		inst->branchtree << "LOCAL internal swap pass " << counter << " : " << currentCost << "|" << endl;
	}
	inst->branchtree << "LOCAL internal swap: " << currentCost << "|" << endl;

	if (false) {
		for (int counter = 0; counter < 1; counter++) {
			vector<int> tempv = tour;
			int todelete = -1;
			for (int i = 0; i < tour.size(); i++) {
				tempv = tour;
				tempv.erase(std::remove(tempv.begin(), tempv.end(), i), tempv.end());
				double newCost;
				if (tours_solved.find(tour) == tours_solved.end()) {
					newCost = tour_solver1(inst, tempv, false);
					tours_solved[tempv] = newCost;
				}
				else {
					tour_repetition_counter++;
					newCost = tours_solved[tempv];
				}
				if (newCost < currentCost) {
					todelete = i;
					currentCost = newCost;
					besttour = tempv;
				}
			}
			tour.erase(std::remove(tour.begin(), tour.end(), todelete), tour.end());
		}

		inst->branchtree << "LOCAL deleting II: " << currentCost << "|";
	}
	inst->branchtree << "LOCAL Best tour: "; 
	for (int i : besttour) inst->branchtree << i << ","; 
	inst->branchtree << endl;

	return currentCost;
}

vector<int> calculate_sorted_c(instance* inst) {

	vector<int> index_sorted = vector<int>(inst->dimension_c);
	
	for (int i = 0; i < inst->dimension_c; i++)
		iota(index_sorted.begin(), index_sorted.end(), 0);

	for (int j = 0; j < inst->dimension_c; j++) {
		sort(index_sorted.begin(), index_sorted.end(), [&](int a, int b) { return inst->u[a] < inst->u[b]; });
	}

	return index_sorted;
}

int createLPCGsub(instance* inst, vector<int> road, vector<int> next_to_avoid) {
	inst->paths.clear();
	inst->path_counter = 0;
	vector<path> initialPaths;


	//Objective
	inst->Lateness = IloAdd(inst->model1, IloMinimize(inst->env));

	//Constraints
	//Constraints: Arc
	inst->zerosArcs = IloNumArray(inst->env);
	for (int i = 0; i < inst->dimension_f; i++) {
		inst->zerosArcs.add(0.0);
	}

	inst->arcconstr = NumConstrMatrix2(inst->env, inst->dimension_f);
	for (IloInt i = 0; i < inst->dimension_f; i++) {
		inst->arcconstr[i] = IloAdd(inst->model1, IloRangeArray(inst->env, inst->zerosArcs, IloInfinity));
	}
	

	//Variables: arc
	inst->arcs = NumVarMatrix(inst->env, inst->dimension_f);
	for (IloInt i = 0; i < inst->dimension_f; i++)
		inst->arcs[i] = IloNumVarArray(inst->env, inst->dimension_f, 0.0, 1.0, ILOFLOAT);


	for (IloInt i = 0; i < inst->dimension_f; i++)
		for (IloInt j = 0; j < inst->dimension_f; j++) {

			//string name = "arc." + to_string(i) + "." + to_string(j);
			inst->arcs[i][j] = IloNumVar(inst->env, 0.0, 1.0, ILOFLOAT);

			//inst->arcs[i][j].setName(name.c_str());

			for (IloInt k = 0; k < inst->dimension_c; k++) {
				inst->arcconstr[i][j].setLinearCoef(inst->arcs[i][j], inst->dimension_c);
			}
			//time to traverse an arc
			//cout << "Arc " << i << "  " << j << " : " << inst->c_f[i][j] << endl;
		}

	int prev = 0;
	for (int ind = 1; ind < road.size();ind++) {
		inst->arcs[prev][road[ind]].setLB(1.0);
		prev = road[ind];
	}


	//prohibited
	for (auto jj : next_to_avoid) {
		inst->arcs[road.back()][jj].setUB(0.0);
	}


	// 
	//Constraints: Arcs II

	inst->arcconstr2 = NumConstrMatrix3(inst->env, 2);
	inst->arcconstr2[0] = NumConstrMatrix2(inst->env, inst->dimension_f);
	inst->arcconstr2[1] = NumConstrMatrix2(inst->env, inst->dimension_f);
	for (IloInt i = 0; i < inst->dimension_f; i++) {

		inst->arcconstr2[0][i] = IloAdd(inst->model1, IloRangeArray(inst->env, inst->zerosArcs, IloInfinity));
		inst->arcconstr2[1][i] = IloAdd(inst->model1, IloRangeArray(inst->env, inst->zerosArcs, IloInfinity));

		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			if (i != 0) {
				inst->arcconstr2[0][i][ii] = IloAdd(inst->model1, inst->arcs[i][ii] <= 0.0);
			}
			else {
				inst->arcconstr2[0][i][ii] = IloAdd(inst->model1, inst->arcs[i][ii] <= 1.0);
			}

			if (ii != 0) {
				inst->arcconstr2[1][i][ii] = IloAdd(inst->model1, inst->arcs[i][ii] <= 0.0);

			}
			else {
				inst->arcconstr2[1][i][ii] = IloAdd(inst->model1, inst->arcs[i][ii] <= 1.0);
			}

		}

	}

	//Constraints: Coverage
	inst->Coverage = IloAdd(inst->model1, IloRangeArray(inst->env, inst->dimension_c, 1.0, IloInfinity));

	//Constraints: Capacity
	IloNumArray facility_capacities(inst->env);
	//inst->capacity_constraint_initialized = vector<bool>(inst->dimension_f, false);
	for (int i = 0; i < inst->dimension_f; i++) {
		facility_capacities.add(inst->capacity[i]);
	}

	inst->Capacity = IloRangeArray(inst->env, inst->dimension_f);

	//Variables
	//Variables: lambda path
	initialPaths = generateInitialPaths11_withconstr_0(inst, road,next_to_avoid);

	if (initialPaths.size() < 1) {
		initialPaths = generateInitialPaths11_withconstr_mip(inst, road, next_to_avoid);
	}

	if (initialPaths.size() < 1) {
		return -1;
	}

	inst->la = IloNumVarArray(inst->env, initialPaths.size(), 0, IloInfinity, ILOFLOAT);
	
	for (int j = 0; j < inst->dimension_f; j++) {

		IloExpr v(inst->env);

		for (int i = 0; i < inst->dimension_f; i++) {
			v -= inst->capacity[j] * inst->arcs[i][j];
		}

		inst->Capacity[j] = IloAdd(inst->model1, IloRange(v <= 0.0));
		v.end();
	}

			
	

	for (int j = 0; j < initialPaths.size(); j++) {

		inst->paths.push_back(initialPaths[j]);
		inst->Lateness.setLinearCoef(inst->la[j], initialPaths[j].getLateness());

		//string name = "la." + to_string(initialPaths[j].getCustomer()) + ".id" + to_string(inst->path_counter);
		inst->path_counter++;
		//inst->la[j].setName(name.c_str());
		inst->Coverage[initialPaths[j].getCustomer()].setLinearCoef(inst->la[j], 1);

		inst->Capacity[initialPaths[j].getWorkingFacility()].setLinearCoef(inst->la[j], inst->demand[initialPaths[j].getCustomer()]);
		
		vector<int> path_ = initialPaths[j].getFacilityPath();
		int pos = 1;
		if (path_.size() > 1) {
			int prev = path_[0];
			for (int i = 1; i < path_.size(); i++) {
				inst->arcconstr[prev][path_[i]].setLinearCoef(inst->la[j], -1);
				pos++;
				prev = path_[i];
			}

		}

		for (int ii = 0; ii < inst->dimension_f; ii++) {
			if (initialPaths[j].getWorkingFacility() != 0) {
				inst->arcconstr2[0][initialPaths[j].getWorkingFacility()][ii].setLinearCoef(inst->la[j], -1.0);
				inst->arcconstr2[1][ii][initialPaths[j].getWorkingFacility()].setLinearCoef(inst->la[j], -1.0);
			}
		}

	}
	initialPaths.clear();

	//X constraints

	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(inst->env);
		IloExpr Xdeltaplus(inst->env);
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			Xdeltaminus += inst->arcs[ii][i];
			Xdeltaplus += inst->arcs[i][ii];
		}

		inst->model1.add(Xdeltaminus <= 1);
		inst->model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			inst->model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}

	
	return 0;
}

LP_return solveCG_LP(instance* inst, vector<int> &road, vector<int> &next_to_avoid) {
	
	clock_t time_start_cycle = clock();
	
	double boost_time = 0.0; 
	double CGLP_time = 0.0;
	double add_cols_time = 0.0;
	int CG_count = 0;

	std::vector<char> inRoad(inst->dimension_f, 0), avoid(inst->dimension_f, 0);
	for (int v : road) if (0 <= v && v < inst->dimension_f) inRoad[v] = 1;
	for (int v : next_to_avoid) if (0 <= v && v < inst->dimension_f) avoid[v] = 1;

	for (;;) {				
		CG_count++;

		clock_t CGLP_start = clock();
		inst->CGsolver.setParam(IloCplex::Param::TimeLimit, max(50.0, inst->timelimit - (double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC));
		
		
		inst->CGsolver.setParam(IloCplex::Param::RandomSeed,100);
		inst->CGsolver.setParam(IloCplex::Param::Threads,1);
		

		inst->CGsolver.solve();

		CGLP_time += (double)(clock() - CGLP_start) / (double)CLOCKS_PER_SEC;

		clock_t boost_start = clock();

		boost_shortest11(inst, road, inRoad, avoid);
		
		boost_time += (double)(clock() - boost_start) / (double)CLOCKS_PER_SEC;
		
		if ((double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC > inst->timelimit) {
			inst->branchtree << "Timelimit. Maybe there are paths with negative reduced cost left." << endl;
			break;
		}

		if (inst->sr_return_has_negative) {
			clock_t add_cols_start = clock();
			
			for (auto p : inst->sr_return_paths) {
				
				//generating column
				IloNumVar tempp = IloAdd(inst->model1, IloNumVar(inst->Lateness(p.getLateness())));
				tempp.setBounds(0.0, IloInfinity);
				//string name = "la." + to_string(p.getCustomer()) + ".id" + to_string(inst->path_counter);
				inst->path_counter++;
				//tempp.setName(name.c_str());
				inst->Coverage[p.getCustomer()].setLinearCoef(tempp, 1);				

					inst->Capacity[p.getWorkingFacility()].setLinearCoef(tempp, inst->demand[p.getCustomer()]);						

				inst->la.add(tempp);
				inst->paths.push_back(p);
				vector<int> path_ = p.getFacilityPath();
				if (path_.size() > 1) {
					int prev = path_[0];
					for (int i = 1; i < path_.size(); i++) {
						inst->arcconstr[prev][path_[i]].setLinearCoef(tempp, -1);
						prev = path_[i];
					}

					for (int ii = 0; ii < inst->dimension_f; ii++) {
						if (p.getWorkingFacility() != 0) {
							inst->arcconstr2[0][p.getWorkingFacility()][ii].setLinearCoef(tempp, -1.0);
							inst->arcconstr2[1][ii][p.getWorkingFacility()].setLinearCoef(tempp, -1.0);
						}
					}				

				}
			}
			add_cols_time += (double)(clock() - add_cols_start) / (double)CLOCKS_PER_SEC;
			
		}
		else
		{
			break;
		}		
	}
	double result = INFINITY;
	
	inst->arcsLPCGsol.clear();
	if (inst->CGsolver.getStatus() == IloAlgorithm::Optimal || inst->CGsolver.getStatus() == IloAlgorithm::Feasible) {
		result = inst->CGsolver.getObjValue();
		try {
			
			for (int i = 0; i < inst->dimension_f; i++) {
				
				inst->arcsLPCGsol.push_back(inst->CGsolver.getValue(inst->arcs[road.back()][i]));
				//cout << inst->arcsLPsol[i] << endl;
			}
		}
		catch (const std::runtime_error& e) {
			inst->CGsolver.exportModel("error_model.lp");
			// Handle the exception
			std::cout << "Exception caught: " << e.what() << std::endl;			
		}	
		
		
	}

	bool status;
	if (inst->CGsolver.getStatus() == IloAlgorithm::Optimal) {
		status = true;
	}
	else
		status = false;

	inst->boost_time_global += boost_time;
	inst->CG_LP_time_global += CGLP_time;
	inst->pricing_number_of_cycles += CG_count;
	inst->branchtree << " Number of lambdas: " << inst->la.getSize() << " ";
	inst->branchtree << " Boost: " << boost_time << " ";
	inst->branchtree << " Solving restricted LP: " << CGLP_time << " ";
	inst->branchtree << " Time to add cols to restricted LP: " << add_cols_time << " ";
	inst->branchtree << " Number of CG cycles: " << CG_count << endl;
	return LP_return(result, status, {});
}

void master2c(instance* inst) {

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 
	inst->env = IloEnv();
	inst->model1 = IloModel(inst->env);
	inst->sorted = calculate_sorted_c(inst);

	//log files
	string filename = inst->input_file + "tree";
	inst->branchtree.open(filename, std::ios_base::app);
	inst->branchtree << endl;
	inst->branchtree << "Solving instance " << filename << " Flags: ";
	if (inst->use_LP_basedMIP)
		inst->branchtree << "use lambda in {0,1}, ";
	else if (inst->use_MIP_basedMIP)
		inst->branchtree << "use x in A_p";
	if (inst->use_LS)
		inst->branchtree << "use Local Search use " << inst->numSol_limit << "solutions";
	inst->branchtree << endl;

	inst->boost_time_global = 0.0;
	inst->CG_LP_time_global = 0.0;
	inst->pricing_number_of_cycles = 0;
	inst->createCGLPtime = 0.0;
	inst->MIPsolvetimeatnode = 0.0;
	inst->total_LS_time = 0.0;

	vector<path> initialPaths;

	inst->CGsolver = IloCplex(inst->model1);

	inst->startglobal = clock();
	clock_t time_end_cycle;
	int ncycles = 0;

	inst->bestMIP = INFINITY;

	inst->CGsolver.setParam(IloCplex::Param::TimeLimit, 3600);

	vector<int> road = { 0 };	

	queue<Lnode> L;
	
	class nodeLess
	{
	public:
		bool operator()(const Lnode a1, const Lnode a2) const {
			return a1.LP > a2.LP;
		}
	};

	vector<Lnode> base;

	if (inst->tree.size() > 1)
		for (const auto& pair : inst->tree) {
			if (pair.second.currentpath.size() > 2 && pair.second.LP <= inst->bestMIPfromcompact) {
				road = pair.second.currentpath;
				double curLP = pair.second.LP;
				base.push_back(Lnode(road, curLP, 1, 1, true, {}, {}));
			}
		}
	else {		
		base.push_back(Lnode(road, 0, 1, 1, true, {}, {}));
	}

	std::priority_queue<Lnode, vector<Lnode>, nodeLess> pL(base.begin(), base.end());

	inst->arcsLPCGsol = vector<double>(inst->dimension_f, 0.0);

	int depth;	

	Lnode lastnode;
	double solution_time_cycle;

	int LPexiststat = 0;
	double LB = 0.0;
	int number_of_nodes = 0;

	while (!pL.empty()) {
		number_of_nodes++;

		LP_return LPSol_CG = LP_return(-1.0, false, {});

		Lnode current_node = pL.top();

		lastnode = current_node;
		pL.pop();
		
		if (current_node.LP > inst->bestMIP) {
			continue; 
		}

		LB = std::min(LB, current_node.LP);

		vector<int> road_local = current_node.road;
		depth = road_local.size();
		inst->branchtree << "Current node: ";
		for (auto v : road_local) inst->branchtree << " " << v;
		inst->branchtree << "\t";
		inst->branchtree << "Parent LP: " << current_node.LP << " best MIP: " << inst->bestMIP << "\t";
		inst->branchtree << "Tree size: " << pL.size() << "Time in cycle: " << (double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC << "\t";
		inst->branchtree << "Prohibit:";
		for (auto v : current_node.prohibited_next) inst->branchtree << " " << v;
		inst->branchtree << "\t";

		clock_t CG_time_start = clock();
		clock_t createCGLP_start = clock();

		inst->CGsolver.clearModel();
		inst->env.end();
		inst->env = IloEnv();
		inst->model1 = IloModel(inst->env);

		clock_t createCG_start = clock();

		LPexiststat = createLPCGsub(inst, road_local, current_node.prohibited_next);
		inst->CGsolver = IloCplex(inst->model1);
		if (LPexiststat != 0) {
			continue;
		}

		inst->createCGLPtime += (double)(clock() - createCGLP_start) / (double)CLOCKS_PER_SEC;

		LPSol_CG = solveCG_LP(inst, road_local, current_node.prohibited_next);

#ifndef impose_integrality

		write_log1(inst, 21, LPSol_CG.LP, 0, ncycles, inst->la.getSize(), (int)inst->CGsolver.getStatus());
		return;

#endif 
		inst->branchtree << "LPCG time : " << (double)(clock() - CG_time_start) / (double)CLOCKS_PER_SEC << " Status: " << LPSol_CG.status << " LPCG obj: " << LPSol_CG.LP << endl;

		double currentLB = INFINITY;
		if (LPSol_CG.status) {
			currentLB = LPSol_CG.LP;
		}

		if (LPSol_CG.status && currentLB > inst->bestMIP)
			if ((double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC > inst->timelimit) {
				break;
			}
			else {
				continue; 
			}

		/*
		bool arcsareint = true;
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt j = 0; j < inst->dimension_f; j++)
				if (i != j) {
					if (inst->CGsolver.getValue(inst->arcs[i][j]) > RC_EPS &&
						inst->CGsolver.getValue(inst->arcs[i][j]) < 1 - RC_EPS) {
						arcsareint = false;
					}
				}
		*/

		time_end_cycle = clock();
		solution_time_cycle = (double)(time_end_cycle - inst->startglobal) / (double)CLOCKS_PER_SEC;
		if (solution_time_cycle > inst->timelimit) {
			break;
		}

		//find the greatest x value
		double maxarc = 1.0;
		int maxnext = -1;

		for (int j = 1; j < inst->dimension_f; j++) {
			if (find(road_local.begin(), road_local.end(), j) == road_local.end()) {
				if (inst->arcsLPCGsol[j] > 1 - RC_EPS)
				{
					maxarc = 1;
					maxnext = j;
					break;
				}				
				double const frac = ::fabs(0.5 - inst->arcsLPCGsol[j]);
				if (maxarc > frac)
				{
					maxarc = frac;
					maxnext = j;
				}
			}
		}

		if (current_node.requiredCG) {			

			if (maxnext >= 0) {
				vector<int> road_local_copy = road_local;
				road_local_copy.push_back(maxnext);
				if (road_local_copy.size() < inst->dimension_f) {
					pL.push(Lnode(road_local_copy, currentLB, 0, 0, true, {}, {}));
					pL.push(Lnode(road_local, currentLB, 0, 0, false, {}, { maxnext }));
				}
			}
		}
		else {

			if (maxnext >= 0) {
				vector<int> road_local_copy = road_local;
				road_local_copy.push_back(maxnext);
				vector<int> prohibit_copy = current_node.prohibited_next;
				prohibit_copy.push_back(maxnext);
				if (road_local_copy.size() < inst->dimension_f) {
					pL.push(Lnode(road_local_copy, currentLB, 0, 0, true, {}, {}));
				}

				if (road_local.size() + prohibit_copy.size() < inst->dimension_f) {
					pL.push(Lnode(road_local, currentLB, 0, 0, false, {}, prohibit_copy));
				}
			}
		}

		inst->branchtree << "Best next is: " << maxnext << " arc value: " << maxarc << endl;

		//SOLVING MIP

		clock_t MIPsolve_start;

		//SOLVE ALLOWING ONLY ARCS IN A_p
		if (inst->use_LP_basedMIP && number_of_nodes % inst->MIPthr == 0) {	
			MIPsolve_start = clock();
			vector<int> best_tour;
			tour_forest_solver1(&best_tour, inst, false);			
			inst->MIPsolvetimeatnode += (double)(clock() - MIPsolve_start) / (double)CLOCKS_PER_SEC;
			inst->branchtree << "Integral arcs MIP_CG time : " << (double)(clock() - MIPsolve_start) / (double)CLOCKS_PER_SEC << endl;
			
			clock_t LS_start = clock();
			if (inst->use_LS && (depth > inst->LS_depth_limit || depth == 1) && number_of_nodes % inst->LSthr == 0) {
				double localsearchsolution = local_search(inst, best_tour);

				if (localsearchsolution < inst->bestMIP) {
					inst->bestMIP = localsearchsolution;
				}
				
				inst->branchtree << " MIP_CG obj LOCAL SEARCH : " << localsearchsolution << endl;
				
			}
			inst->total_LS_time += (double)(clock() - LS_start) / (double)CLOCKS_PER_SEC;
			inst->branchtree << "LS time : " << (double)(clock() - LS_start) / (double)CLOCKS_PER_SEC << endl;

		}
		//END NON-ZERO ARCS IN A_p

		//SOLVING MIP WITH LAMBDA IN {0,1}
		if (inst->use_MIP_basedMIP && number_of_nodes % inst->MIPthr == 0) {

			MIPsolve_start = clock();

			IloConversion typeconvcycle = IloConversion(inst->env, inst->la, IloNumVar::Int);

			inst->model1.add(typeconvcycle);

			for (IloInt i = 0; i < inst->dimension_f; i++)
				for (IloInt j = 0; j < inst->dimension_f; j++)
					if (i != j) {
						inst->model1.add(IloConversion(inst->env, inst->arcs[i][j], ILOINT));
					}

			inst->CGsolver.setParam(IloCplex::Param::TimeLimit, max(50.0, inst->timelimit - (double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC));

			inst->CGsolver.solve();

			inst->branchtree << "Integral arcs MIP_CG time : " << (double)(clock() - MIPsolve_start) / (double)CLOCKS_PER_SEC << "\t";
			
			if (inst->CGsolver.getStatus() == IloAlgorithm::Optimal || inst->CGsolver.getStatus() == IloAlgorithm::Feasible) {
				if (inst->CGsolver.getObjValue() < inst->bestMIP)
					inst->bestMIP = inst->CGsolver.getObjValue();
			}
			
			inst->branchtree << " MIP_CG status : " << inst->CGsolver.getStatus() << endl;

			inst->MIPsolvetimeatnode += (double)(clock() - MIPsolve_start) / (double)CLOCKS_PER_SEC;


			//LOCAL SEARCH
			if (inst->use_LS && (depth > inst->LS_depth_limit || depth == 1) && number_of_nodes % inst->LSthr == 0) {
				int numSolutions = inst->CGsolver.getSolnPoolNsolns();
				int counter_sol = 0;
				vector<vector<int>> prevtours;
				for (int nSol = 0; nSol < numSolutions && counter_sol < inst->numSol_limit; nSol++) {
					inst->branchtree << "Solution: " << nSol << " Obj: " << inst->CGsolver.getObjValue(nSol) << "\t";

					//construct starttour
					vector<int> starttour;
					int previ = 0;
					starttour.push_back(previ);
					set<int> nextcandidate;
					for (int i = 1; i < inst->dimension_f; i++)
						nextcandidate.insert(i);
					double accruedcapacity = 0.0;
					accruedcapacity += inst->capacity[0];
					while (starttour.size() < inst->dimension_f) {
						int next;
						double max = -INFINITY;
						for (auto i : nextcandidate) {
							if (inst->CGsolver.getValue(inst->arcs[previ][i], nSol) > max) {
								max = inst->CGsolver.getValue(inst->arcs[previ][i], nSol);
								next = i;
							}
						}
						starttour.push_back(next);
						accruedcapacity += inst->capacity[next];
						previ = next;
						nextcandidate.erase(next);
					}

					//check if diverse
					bool isdiverse = true;
					for (auto prevt : prevtours) {
						if (calculateModifiedHammingDistance(prevt, starttour) < 10) {
							isdiverse = false;
							break;
						}
					}

					if (isdiverse) {
						counter_sol++;
						prevtours.push_back(starttour);

						double localsearchsolution = local_search(inst, starttour);

						if (localsearchsolution < inst->bestMIP) {
							inst->bestMIP = localsearchsolution;
							inst->branchtree << " MIP_CG obj LOCAL SEARCH : " << localsearchsolution << endl;
						}
						else {
							inst->branchtree << " MIP_CG obj LOCAL SEARCH : " << localsearchsolution << endl;
						}
					}
					else {
						inst->branchtree << "Not diverse" << endl;
					}
				}
			}

		}
		//END solving MIP lambda in {0,1}

	}

	clock_t finalMIP_start = clock();

	inst->CGsolver.clearModel();
	inst->env.end();
	inst->env = IloEnv();
	inst->model1 = IloModel(inst->env);
		
	if (pL.empty()) {
		// Full enumeration: optimal solution found; global LB = bestMIP
		LB = inst->bestMIP;
	}
	else {
		LB = std::numeric_limits<double>::infinity();
		while (!pL.empty()) {
			Lnode cand = pL.top();
			pL.pop();
			if (cand.LP < LB) {
				LB = cand.LP;
			}
		}
	}

	inst->path_counter = 0;

	double solution_time = (double)(clock() - inst->startglobal) / (double)CLOCKS_PER_SEC;
	inst->branchtree << "Best LP: " << LB << endl;
	write_log21(inst, 22, 0, inst->bestMIP, LB, solution_time, ncycles, number_of_nodes, 0);

	inst->branchtree << "boost time: " << inst->boost_time_global << endl;
	inst->branchtree << "solving restricted LP time: " << inst->CG_LP_time_global << endl;
	inst->branchtree << "Number of pricing cycles: " << inst->pricing_number_of_cycles << endl;
	inst->branchtree << "create CGLP at nodes: " << inst->createCGLPtime << endl;
	inst->branchtree << "MIP solve at nodes if int arcs: " << inst->MIPsolvetimeatnode << endl;
	inst->branchtree << "Local Search time: " << inst->total_LS_time << endl;

	//IloNumArray solution = IloNumArray(inst->env);
	//pathSolver.getValues(inst->la, solution);
	/*
	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			if (inst->CGsolver.getValue(inst->arcs[i][j]) > 1 - RC_EPS) {


				inst->active_arcs.push_back(pairs(i, j));

			}
		}

	vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {
		solution[i] = vector<double>(inst->dimension_c,0.0);
	}


	for (IloInt i = 0; i < inst->la.getSize(); i++)
		if (inst->CGsolver.getValue(inst->la[i]) != 0.0)
		{
			solution[inst->paths[i].getWorkingFacility()][inst->paths[i].getCustomer()] = 1.0;
		}


	generate_latex_model1(inst, solution);
	*/
	inst->branchtree << "Tours repeted in local search: " << tour_repetition_counter;
	inst->branchtree << "Solved paths: " << tours_solved.size();
}



double dummy(instance* inst) {
	vector<int> tour;
	vector<bool> visited(inst->dimension_f, false);

	int currentNode = 0;
	tour.push_back(currentNode);
	visited[currentNode] = true;

	for (int i = 1; i < inst->dimension_f; i++) {
		int nearestNode = findNearestNode(inst, currentNode, visited);
		if (nearestNode == -1) break; // If no nearest node is found (should not happen in a complete graph)
		tour.push_back(nearestNode);
		visited[nearestNode] = true;
		currentNode = nearestNode;
	}
	
	return tour_solver1(inst,tour,false);
}

vector<path> generateInitialPaths11_withconstr_0(instance* inst, vector<int>road_to_take, vector<int> next_to_avoid) {

	double lateness = 0.0;
	vector<path> paths = vector<path>();
	vector<int> nearestF = vector<int>(inst->dimension_c);
	vector<set<int>> serve = vector<set<int>>(inst->dimension_f);
	set<int> opened_facilities;
	vector<double> depleted_capacity = inst->capacity;
	set<int> vert_in_road;
	set<int> customers_to_serve;
	double mincapacity = INFINITY;
	int mincap;
	vector<int> all = vector<int>(inst->dimension_f);
	iota(all.begin(), all.end(), 0);
	//first use obligatory facilities (current path)
	for (auto v : road_to_take) {
		vert_in_road.insert(v);
	}
	vert_in_road.erase(0);

	for (int k = 0; k < inst->dimension_c; k++) {

		customers_to_serve.insert(k);

	}
	//we need to use all facilities in current path
	while (!vert_in_road.empty()) {
		mincapacity = INFINITY; //select facility with min capacity
		for (auto v : vert_in_road) {
			if (inst->capacity[v] < mincapacity) {
				mincapacity = inst->capacity[v];
				mincap = v;
			}
		}
		vert_in_road.erase(mincap);
		bool found = false; //assign a customer if possible
		for (auto k : customers_to_serve) {
			if (inst->demand[k] <= mincapacity) {
				found = true;
				depleted_capacity[mincap] -= inst->demand[k];
				nearestF[k] = mincap;
				serve[mincap].insert(k);
				opened_facilities.insert(mincap);
				customers_to_serve.erase(k);
				break;
			}
		}
		if (!found)
			return paths;
	}

	//liason facilities
	bool success = true;// at least one liason is open
	if (next_to_avoid.size() > 0) {
		success = false;
		std::vector<int> liason;
				
		std::set<int>temp;
		temp.insert(road_to_take.cbegin(), road_to_take.cend());
		temp.insert(next_to_avoid.cbegin(), next_to_avoid.cend());

		for (int i = 0; i < inst->dimension_f;i++) {
			if (temp.find(i) == temp.end()) {
				liason.push_back(i);
			}
		}
		 
		for (auto k : customers_to_serve) {

			//nearest facility
			int nearest = -1;
			double incumbent = IloInfinity;
			for (auto i:liason)
			{
				if (inst->c_fc[i][k] < incumbent && depleted_capacity[i] >= inst->demand[k]) {
					nearest = i;
					incumbent = inst->c_fc[i][k];
				}				
			}
			if (nearest == -1) {
				continue;
			}
			else {
				depleted_capacity[nearest] -= inst->demand[k];
				nearestF[k] = nearest;
				serve[nearest].insert(k);
				opened_facilities.insert(nearest);
				customers_to_serve.erase(k);
				success = true;
			}
			if (success) break;
		}
	}
	vector<int> path_;
	if (success) {
		//allocate remaining customers using all available facilities
		for (auto k : customers_to_serve) {

			//nearest facility
			int nearest = -1;
			double incumbent = IloInfinity;
			for (int i = 0; i < inst->dimension_f; i++)
			{
				if (inst->c_fc[i][k] < incumbent && depleted_capacity[i] >= inst->demand[k]) {
					nearest = i;
					incumbent = inst->c_fc[i][k];
				}
				//cout << "From facility " << i << " to customer " << k << ": " << inst->c_fc[i][k] << endl;
			}
			if (nearest == -1) {

				cout << "Instance is infeasible. Not enough capacity.";
				return {};

			}
			depleted_capacity[nearest] -= inst->demand[k];
			nearestF[k] = nearest;
			serve[nearest].insert(k);
			opened_facilities.insert(nearest);
		}

		vector<int> nearest_facility = vector<int>(inst->dimension_c);

		//root node
		
		path_.push_back(0);
		int prev = 0;
		if (!serve[0].empty()) {

			for (auto j : serve[0]) {

				path p = path(j, path_, inst);
				lateness += p.getLateness();
				paths.push_back(p);

			}

			opened_facilities.erase(0);
		}

		for (auto ind : road_to_take) if (ind != 0) {

			path_.push_back(ind);
			for (auto j : serve[ind]) {

				path p = path(j, path_, inst);
				lateness += p.getLateness();
				paths.push_back(p);
				//cout << "path to customer " << j;
				//for (int i = 0; i < path_.size(); i++)
				//	cout << "  " << path_[i];
				//cout << endl;
			}
			opened_facilities.erase(ind);
			prev = ind;
		}
		bool is_next_after_road = true;
		while (!opened_facilities.empty()) {

			//next nearest facility

			double minnext = INFINITY;
			int next;

			if (is_next_after_road) {
				is_next_after_road = false;
				for (auto i : opened_facilities)
					if (find(next_to_avoid.begin(), next_to_avoid.end(), i) == next_to_avoid.end() && inst->c_f[prev][i] < minnext) {
						next = i;
						minnext = inst->c_f[prev][i];
					}
					else {
						if (opened_facilities.size() == 1) {
							break;
						}
					}
			}
			else
			{
				for (auto i : opened_facilities)
					if (inst->c_f[prev][i] < minnext) {
						next = i;
						minnext = inst->c_f[prev][i];
					}
			}

			path_.push_back(next);

			for (auto j : serve[next]) {

				path p = path(j, path_, inst);
				paths.push_back(p);
				lateness += p.getLateness();
				//cout << "path to customer " << j;
				//for (int i = 0; i < path_.size(); i++)
				//	cout << "  " << path_[i];
				//cout << endl;
			}

			opened_facilities.erase(next);
			prev = next;
		}

		if (lateness < inst->bestMIP)
			inst->bestMIP = lateness;
	}

	//II set
	vector<set<int>> servedbyfacility(inst->dimension_f);
	vector<double> disttofacility = vector<double>(inst->dimension_f);
	vector<int> servingfacility = vector<int>(inst->dimension_c);
	set<int> roadv, notinroad;
	vector<double> resid_capacity = inst->capacity;
	vector<double> dist_in_road = vector<double>(road_to_take.size());

	int prev = 0;
	dist_in_road[0] = 0.0;
	path_.clear();
	path_.push_back(0);
	for (int ind = 1; ind < road_to_take.size(); ind++) {
		int v = road_to_take[ind];
		dist_in_road[ind] = dist_in_road[ind - 1] + inst->c_f[prev][v];
		path_.push_back(v);
		disttofacility[v] = dist_in_road[ind];
		prev = v;
	}

	lateness = 0.0;
	for (int i = 0; i < inst->dimension_f; i++)
		notinroad.insert(i);

	for (auto v : road_to_take) {
		roadv.insert(v);
		notinroad.erase(v);
	}

	prev = road_to_take[road_to_take.size() - 1];
	bool is_next_after_road = true;
	while (!notinroad.empty()) {
		//next nearest facility
		double minnext = INFINITY;
		int next = -1;

		if (is_next_after_road) {
			is_next_after_road = false;
			for (auto i : notinroad)
				if (find(next_to_avoid.begin(), next_to_avoid.end(), i) == next_to_avoid.end() && inst->c_f[prev][i] < minnext) {
					next = i;
					minnext = inst->c_f[prev][i];
				}
		}
		else {
			for (auto i : notinroad)
				if (inst->c_f[prev][i] < minnext) {
					next = i;
					minnext = inst->c_f[prev][i];
				}
		}
		disttofacility[next] = disttofacility[prev] + minnext;
		path_.push_back(next);
		notinroad.erase(next);
		prev = next;
	}


	for (int ind = 0; ind < inst->dimension_c; ind++) {
		int k = inst->sorted[ind];
		bool attachedtofacility = false;
		double mindist = INFINITY;
		for (int i = 0; i < inst->dimension_f; i++) {
			if (disttofacility[i] + inst->c_fc[i][k] < mindist && resid_capacity[i] >= inst->demand[k]) {
				attachedtofacility = true;
				mindist = disttofacility[i] + inst->c_fc[i][k];
				servingfacility[k] = i;
			}
		}
		if (attachedtofacility) {
			resid_capacity[servingfacility[k]] -= inst->demand[k];
			servedbyfacility[servingfacility[k]].insert(k);
		}
	}

	vector<int> cleanpath;
	for (auto v : path_) {
		if (servedbyfacility[v].size() == 0 && v != 0) {}
		else {
			cleanpath.push_back(v);
		}
	}

	for (int k = 0; k < inst->dimension_c; k++) {
		vector<int> p_(find(cleanpath.begin(), cleanpath.end(), servingfacility[k]) - cleanpath.begin() + 1);
		copy(cleanpath.begin(), find(cleanpath.begin(), cleanpath.end(), servingfacility[k]) + 1, p_.begin());
		//path temppath = path(k, p_, inst);
		paths.push_back(path(k, p_, inst));
		lateness += max(disttofacility[servingfacility[k]] + inst->c_fc[servingfacility[k]][k] - inst->u[k], 0.0);
		//cout << lateness<<"  " <<temppath.getLateness()<<endl;
	}

	if (lateness < inst->bestMIP) {
		inst->bestMIP = lateness;
		inst->bestheur = paths;
	}

	return paths;
}

vector<path> generateInitialPaths11_withconstr_mip(instance* inst, vector<int>road_to_take, vector<int> next_to_avoid) {

	//construct tour greedy nearest
	vector<int> tour;
	set<int> candidates;
	for (int i = 0; i < inst->dimension_f; i++)
		candidates.insert(i);

	for (auto i : road_to_take) {
		tour.push_back(i);
		candidates.erase(i);
	}

	int prev = tour.back();
	double mindist = INFINITY;
	int next = -1;
	for (int i:candidates) {
		if (inst->c_f[prev][i] < mindist && find(next_to_avoid.begin(),next_to_avoid.end(), i) == next_to_avoid.end()) {
			mindist = inst->c_f[prev][i];
			next = i;
		}
	}
	if (next > 0) {
		tour.push_back(next);
		candidates.erase(next);
	}
	prev = tour.back();
	while (tour.size() < inst->dimension_f) {
		mindist = INFINITY;
		next = -1;
		for (auto i : candidates) {
			if (inst->c_f[prev][i] < mindist) {
				mindist = inst->c_f[prev][i];
				next = i;
			}
		}
		if (next > 0) {
			tour.push_back(next);
			candidates.erase(next);
			prev = next;
		}
	}

	IloEnv env = IloEnv();
	IloModel model1 = IloModel(env);

	IloNumVar::Type inttype = IloNumVar::Int;

	//IloNumVar::Type inttype = IloNumVar::Float;


	//Objective
	IloObjective  Lateness = IloAdd(model1, IloMinimize(env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(env, inst->dimension_c, 0, 1, inttype);
		/*
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());


		}
		*/
	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(env, inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			/*
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());

			}*/
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		//string name = "s." + to_string(k);
		//s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			//string name = "arc." + to_string(i) + "." + to_string(j);

			arcs[arc_coef(i, j, inst)] = IloNumVar(env, 0.0, 0.0, inttype);

			//arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}


	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {

		IloExpr v(env);
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += inst->c_f[i][ii] * f[k][i][ii];

			}

		for (IloInt i = 0; i < inst->dimension_f; i++) {

			v += z[i][k] * inst->c_fc[i][k];

		}

		model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}

			if (i == 0)
				model1.add(v == -1 + z[0][k]);
			else
				model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++)
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				model1.add(f[k][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(env);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				//v += arcs[ii][i];
				v += arcs[arc_coef(ii, i, inst)];

			}
			model1.add(z[i][k] <= v);
			v.end();

		}

	for (int k = 0; k < inst->dimension_c; k++) {
		IloExpr v(env);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			v += z[i][k];
		}
		model1.add(v == 1.0);
	}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		model1.add(v <= inst->capacity[i]);
		v.end();
	}

	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(env);
		IloExpr Xdeltaplus(env);
		for (int ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		model1.add(Xdeltaminus <= 1);
		model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}
	prev = 0;
	for (auto v : tour) if (v != 0) {
		arcs[arc_coef(prev, v, inst)].setUB(1.0);
		model1.add(arcs[arc_coef(prev, v, inst)] == 1);
		prev = v;
	}
	arcs[arc_coef(tour.back(), 0, inst)].setUB(1.0);

	//condition 3 (all open facilities have customers)
	for (auto i : road_to_take) if(i!=0){
		IloExpr zsumi(env);
		for (IloInt ii = 0; ii < inst->dimension_c; ii++) {
			zsumi += z[i][ii];
		}

		model1.add(zsumi >= 1);
		zsumi.end();
	}

	IloCplex pathSolver(model1);
	//string filename = inst->input_file + "solution";
	//pathSolver.readMIPStarts(filename.c_str());
	//pathSolver.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	pathSolver.setParam(IloCplex::Param::TimeLimit, 600);
	
	//pathSolver.exportModel("finalMIP1.lp");
	bool feasible = pathSolver.solve();
	vector<path> solution;
	if (feasible) {
		double obj = pathSolver.getObjValue();
		if (obj < inst->bestMIP)
			inst->bestMIP = obj;

		//generate paths

		vector<int> path_;
		int prev = 0;
		//zero paths
		path_.push_back(0);
		for (int k = 0; k < inst->dimension_c; k++) {
			if (pathSolver.getValue(z[0][k]) > 0.999) {
				path p = path(k, path_, inst);
				solution.push_back(p);
			}
		}

		for (auto i : tour) if (i != 0) {
			if (pathSolver.getValue(arcs[arc_coef(prev, i, inst)]) > 0.999) {
				path_.push_back(i);
				bool found = false;
				for (int k = 0; k < inst->dimension_c; k++) {
					if (pathSolver.getValue(z[i][k])>0.999) {
						path p = path(k, path_, inst);
						solution.push_back(p);
						found = true;
					}
				}
				prev = i;
				if (!found) path_.pop_back();
			}
		}

		if (false) {
			cout << "ARCS Values: " << endl;
			for (int i = 0; i < inst->dimension_f; i++)
				for (int j = 0; j < inst->dimension_f; j++) {
					if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {
						inst->active_arcs.push_back(pairs(i, j));
						cout << i << " " << j << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << ", ";
					}
				}

			vector<vector<double>> solutions = vector<vector<double>>(inst->dimension_f);
			cout << "Assignment Values: " << endl;
			for (int i = 0; i < inst->dimension_f; i++) {
				solutions[i] = vector<double>(inst->dimension_c);
				for (int k = 0; k < inst->dimension_c; k++) {
					if (pathSolver.getValue(z[i][k]) > RC_EPS)
						cout << "z" << i << " " << k << " = " << pathSolver.getValue(z[i][k]) << ", ";
					solutions[i][k] = pathSolver.getValue(z[i][k]);
				}
			}

			generate_latex_model1(inst, solutions);

		}
	}
	else
	{
		solution = {};
	}

	pathSolver.clearModel();
	env.end();

	return solution;
}

