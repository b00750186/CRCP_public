 #include "model2.h"

vector<path> generateInitialPaths(instance * inst);
vector<path> generateInitialPaths_fixed(instance* inst);

void print_path(int i, instance * inst) {



	vector<int> visitedfacilities = inst->paths[i].getFacilityPath();
	for (int j = 0; j < visitedfacilities.size(); j++) {
	
		cout << visitedfacilities[j] <<" ";
	
	}
	cout << endl;
}
void master2(instance* inst) {

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 
	inttype = IloNumVar::Float;
	inst->model1 =  IloModel(inst->env);
	
	//Objective
	IloObjective  Lateness = IloAdd(inst->model1, IloMinimize(inst->env));
	
	//Constraints
	//Constraints: Arc
	IloNumArray zerosArcs(inst->env);
	for (int i = 0; i < inst->dimension_f; i++) {
		zerosArcs.add(0.0);
	}

	NumConstrMatrix3 arcconstr(inst->env, inst->dimension_c);	
	for (IloInt k = 0; k < inst->dimension_c; k++) {
		arcconstr[k] = NumConstrMatrix2(inst->env, inst->dimension_f);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			
			arcconstr[k][i] = IloAdd(inst->model1, IloRangeArray(inst->env, zerosArcs, IloInfinity));

		}
		
	}
	//Constraints: Coverage
	IloNumArray ones(inst->env);
	for (int j = 0; j < inst->dimension_c; j++) {
		ones.add(1.0);
	}
	IloRangeArray  Coverage = IloAdd(inst->model1, IloRangeArray(inst->env, ones, IloInfinity));
	

	//Constraints: Capacity
	IloNumArray facility_capacities(inst->env);
	vector<bool> capacity_constraint_initialized = vector<bool>(inst->dimension_f, false);
	for (int i = 0; i < inst->dimension_f; i++) {
		facility_capacities.add(inst->capacity[i]);
	}
	//IloRangeArray  Capacity = IloAdd(inst->model1, IloRangeArray(inst->env, 0, facility_capacities));
	//IloNumVarArray dummy = IloNumVarArray(inst->env, inst->dimension_f, 0, 1, ILOFLOAT);
	IloRangeArray  Capacity = IloRangeArray(inst->env,inst->dimension_f);
	/*
	for (IloInt j = 0; j < inst->dimension_f; ++j) {

		//Capacity[j] =  IloAdd(inst->model1,IloRange(dummy[j] <= inst->capacity[j]));
		Capacity[j] = IloAdd(inst->model1, IloRange());

	}
	*/
	//Variables
	//Variables: lambda path
	
	int path_counter = 0;
	vector<path> initialPaths = generateInitialPaths(inst);
	
	
	IloNumVarArray la = IloNumVarArray(inst->env, initialPaths.size(), 0, IloInfinity, ILOFLOAT);
	
	for (int j = 0; j < initialPaths.size(); j++) {
	
		inst->paths.push_back(initialPaths[j]);
		Lateness.setLinearCoef(la[j], initialPaths[j].getLateness());

		string name = "la." + to_string(initialPaths[j].getCustomer())+".id"+to_string(path_counter);
		path_counter++;
		la[j].setName(name.c_str());
		Coverage[initialPaths[j].getCustomer()].setLinearCoef(la[j],1);

		if (!capacity_constraint_initialized[initialPaths[j].getWorkingFacility()]) {
			
			Capacity[initialPaths[j].getWorkingFacility()] = IloAdd(inst->model1, IloRange(inst->demand[initialPaths[j].getCustomer()] * la[j] <= inst->capacity[initialPaths[j].getWorkingFacility()]));
			capacity_constraint_initialized[initialPaths[j].getWorkingFacility()] = true;
		}
		else {
			Capacity[initialPaths[j].getWorkingFacility()].setLinearCoef(la[j], inst->demand[initialPaths[j].getCustomer()]);
		}
		vector<int> path_ = initialPaths[j].getFacilityPath();
		if (path_.size() > 1) {
			int prev = path_[0];
			for (int i = 1; i < path_.size(); i++) {
				arcconstr[initialPaths[j].getCustomer()][prev][path_[i]].setLinearCoef(la[j], -1);
				prev = path_[i];
			}
		}
		
	}
	initialPaths.clear();
	//Variables: arc
	NumVarMatrix arcs(inst->env, inst->dimension_f);
	for (IloInt i = 0; i < inst->dimension_f; i++)
		arcs[i] = IloNumVarArray(inst->env, inst->dimension_f, 0.0, 1.0, ILOFLOAT);


	for (IloInt i = 0; i < inst->dimension_f; i++)
	for (IloInt j = 0; j < inst->dimension_f; j++) {
		
		string name = "arc." + to_string(i) + "." + to_string(j);
		arcs[i][j] = IloNumVar(inst->env,0.0,1.0,ILOFLOAT);
		arcs[i][j].setName(name.c_str());

		for (IloInt k = 0; k < inst->dimension_c; k++) {
			arcconstr[k][i][j].setLinearCoef(arcs[i][j],1.0);
		}
		//time to traverse an arc
		cout << "Arc " << i << "  " << j << " : " << inst->c_f[i][j] << endl;
	}

	//X constraints
	
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(inst->env);
		IloExpr Xdeltaplus(inst->env);
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			Xdeltaminus += arcs[ii][i];
			Xdeltaplus += arcs[i][ii];
		}

		inst->model1.add(Xdeltaminus <= 1);
		inst->model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			inst->model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}


	//fix some arcs
	//inst->model1.add(arcs[0][10] == 1);
	//inst->model1.add(arcs[10][2] == 1);
	/*
	//SEC MTZ
	IloNumVarArray t = IloNumVarArray(inst->env, inst->dimension_f, 0, inst->dimension_f, ILOFLOAT);

	for (int i = 0; i < inst->dimension_f; i++) {
		string name = "t." + to_string(i);
		t[i].setName(name.c_str());
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int ii = 0; ii < inst->dimension_f; ii++) {

			inst->model1.add(t[ii] >= t[i] - inst->dimension_f + inst->dimension_f * arcs[i][ii]);

		}
	*/
	IloCplex pathSolver(inst->model1);
	pathSolver.exportModel("model1_path_toy.lp");
	
	clock_t time_start = clock();
	int ncycles = 0;
	for (;;) {
		
		ncycles++;

		pathSolver.solve();

		//Printing solution and reduced costs
		//cout << "Solution:";
		//cout << "Objective: " << pathSolver.getObjValue() << endl;
		//for (IloInt i = 0; i < la.getSize(); i++) {
		//	std::cout << "  Lambda: " << la[i].getName() << " = " << pathSolver.getValue(la[i]) << " \t Reduced cost: "<<pathSolver.getReducedCost(la[i])<<std::endl;
		//	cout << "Facilities: "; 
		//	print_path(i, inst);
		//}


		//for (IloInt i = 0; i < inst->dimension_f; i++)
		//	for (IloInt j = 0; j < inst->dimension_f; j++) {
		//		
		//		std::cout << "Arc: " << arcs[i][j].getName() << " = " << pathSolver.getValue(arcs[i][j]) << endl;
		//		
		//}
		

		IloNumArray dualCoverage(inst->env, inst->dimension_c);
		cout << "Dual Coverage:"<<endl;
		for (int j = 0; j < inst->dimension_c; j++) {
			dualCoverage[j] = pathSolver.getDual(Coverage[j]);
			//cout << dualCoverage[j] <<endl;
		
		}
		IloNumArray dualCapacity(inst->env, inst->dimension_f);
		cout << "Dual Capacity:" << endl;
		for (int i = 0; i < inst->dimension_f; i++) {
			if (capacity_constraint_initialized[i])
				dualCapacity[i] = pathSolver.getDual(Capacity[i]);
			else
				dualCapacity[i] = 0.0;
			//cout << dualCapacity[i] << endl;
		}
		IloArray<NumMatrix> dualArcs(inst->env, inst->dimension_c);
		
		for (int k = 0; k < inst->dimension_c; k++) {
		
			dualArcs[k] = NumMatrix(inst->env, inst->dimension_f);
			for (int i = 0; i < inst->dimension_f; i++)
				dualArcs[k][i] = IloNumArray(inst->env, inst->dimension_f);

		}
		
		//cout << "Dual Arcs:" << endl;
		for (int k = 0; k < inst->dimension_c; k++) {
			//cout << "dual of arcs for customer " << k << ":" << endl;
			for (int i = 0; i < inst->dimension_f; i++)
				for (int ii = 0; ii < inst->dimension_f; ii++) {

					dualArcs[k][i][ii] = pathSolver.getDual(arcconstr[k][i][ii]);
					/*if(dualArcs[k][i][ii] < -RC_EPS || dualArcs[k][i][ii] > RC_EPS)
						cout << i << "  " << ii << "  :" << dualArcs[k][i][ii] << endl;*/
				}
		}
		
		subproblem_return_single sr;
		bool has_negative = false;
		for (int k = 0; k < inst->dimension_c; k++) {
			
			//cout << "start searching subproblem path: " << k <<endl;
			sr = boost_shortest(k, inst, dualCapacity, dualArcs, dualCoverage);
			//cout << "finished searching" << endl;
			
			if (sr.reduced_cost > -RC_EPS) {

				continue;

			}
			else {
				has_negative = true;
				//generating column
				IloNumVar tempp = IloAdd(inst->model1, IloNumVar(Lateness(sr.p.getLateness())));
				tempp.setBounds(0.0, IloInfinity);
				string name = "la." + to_string(k) + ".id" + to_string(path_counter);
				path_counter++;
				tempp.setName(name.c_str());
				Coverage[k].setLinearCoef(tempp, 1);


				if (!capacity_constraint_initialized[sr.p.getWorkingFacility()]) {

					Capacity[sr.p.getWorkingFacility()] = IloAdd(inst->model1, IloRange(inst->demand[sr.p.getCustomer()] * tempp <= inst->capacity[sr.p.getWorkingFacility()]));
					capacity_constraint_initialized[sr.p.getWorkingFacility()] = true;
				}
				else {

					Capacity[sr.p.getWorkingFacility()].setLinearCoef(tempp, inst->demand[sr.p.getCustomer()]);

				}

				la.add(tempp);
				inst->paths.push_back(sr.p);
				vector<int> path_ = sr.p.getFacilityPath();
				if (path_.size() > 1) {
					
					int prev = path_[0];
					for (int i = 1; i < path_.size(); i++) {

						arcconstr[k][prev][path_[i]].setLinearCoef(tempp, -1);
						prev = path_[i];
					}
				}
				//pathSolver.exportModel("update.lp");
				//cin.get();
			}
		}

		if (!has_negative)
			break;
		clock_t time_end_cycle = clock();
		double solution_time_cycle = (double)(time_end_cycle - time_start) / (double)CLOCKS_PER_SEC;
		if (solution_time_cycle > inst->timelimit)
			break;
	
	}
	
	inst->model1.add(IloConversion(inst->env, la, inttype));
	for (int i = 0; i < la.getSize(); i++)
		la[i].setBounds(0,1);

	for (IloInt i = 0; i < inst->dimension_f; i++)
		for (IloInt j = 0; j < inst->dimension_f; j++) {

			inst->model1.add(IloConversion(inst->env,arcs[i][j], inttype));
		}
	
	pathSolver.solve();
	
	clock_t time_end = clock();
	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	pathSolver.exportModel("finalMIP2.lp");
	std::cout << "Solution status: " << pathSolver.getStatus() << std::endl;
	std::cout << "Solution:  "<< pathSolver.getObjValue() << " " << std::endl;
	std::cout << std::endl;
	for (int i = 0; i < la.getSize(); i++) {
		//std::cout << la[i].getName() << " = " << pathSolver.getValue(la[i]) <<std::endl;
		if (pathSolver.getValue(la[i]) > 0.0) {
			std::cout << la[i].getName() << " = " << pathSolver.getValue(la[i]) << " Customer :" << inst->paths[i].getCustomer() << " Lateness: " << inst->paths[i].getLateness() << std::endl;
			print_path(i, inst);
		}
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			if (pathSolver.getValue(arcs[i][j]) > 1 - RC_EPS) {
			
				
				inst->active_arcs.push_back(pairs(i, j));
			
			}
			if(pathSolver.getValue(arcs[i][j])>0.0)
				std::cout << "Arc: " << arcs[i][j].getName() << " = " << pathSolver.getValue(arcs[i][j]) << endl;

		}
	write_log1(inst, 2, pathSolver.getObjValue(), solution_time,ncycles,(int)la.getSize(),(int)pathSolver.getStatus());
	IloNumArray solution = IloNumArray(inst->env);
	pathSolver.getValues(la,solution);
	generate_latex_model2(inst,solution);
}


vector<path> generateInitialPaths(instance* inst) {

	vector<path> paths;
	vector<double> depleted_capacity = inst->capacity;
	
	for (int j = 0; j < inst->dimension_c; j++) {
		
		//nearest facility
		int nearest = -1;
		double incumbent = IloInfinity;
		for (int i = 0; i < inst->dimension_f; i++) {
			if (inst->c_fc[i][j] < incumbent && depleted_capacity[i] >= inst->demand[j]) {
				nearest = i;
				incumbent = inst->c_fc[i][j];
			}
			cout << "From facility " << i << " to customer " << j << ": " << inst->c_fc[i][j] << endl;
		}
		if (nearest == -1) {
		
			cout << "Instance is infeasible. Not enough capacity.";
		
		}
		depleted_capacity[nearest] -= inst->demand[j];

		//greedy search from the root to the nearest facility
		vector<int> path_;
		vector<bool> visited = vector<bool>(inst->dimension_f, false);
		path_.push_back(0);
		int next = 0;
		int curr = 0;
		visited[0] = true;
		
		while (curr != nearest) {
			double minnext = INFINITY;
			for (int i = 1; i < inst->dimension_f; i++) {
				if (!visited[i]) {
					if (inst->c_f[curr][i] < minnext) {
						minnext = inst->c_f[curr][i];
						next = i;
					}
				}
			}
			
			curr = next;
			visited[curr] = true;
			path_.push_back(curr);
		}
		cout << "path to customer " << j;
		for (int i = 0; i < path_.size(); i++)
			cout << "  " << path_[i];
		cout << endl;
		path p = path(j, path_,inst);
		paths.push_back(p);
	}
	return paths;
}


