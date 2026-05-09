#include "model1.h"

// Generic callback that implements most infeasible branching.

class BranchCallback : public IloCplex::Callback::Function {
	IloNumVarArray x;
	int calls;
	int branches;
	instance* inst;
	int lastvisited;
public:
	BranchCallback(IloNumVarArray _x,instance* inst) : x(_x), calls(0), branches(0), inst(inst) {
	}

	void invoke(IloCplex::Callback::Context const& context) ILO_OVERRIDE {
		++calls;

		IloInt depth = context.getLongInfo(IloCplex::Callback::Context::Info::NodeDepth);
		/*if (depth > 1000) {
			context.pruneCurrentNode();
			return;
		}*/

		IloCplex::CplexStatus status = context.getRelaxationStatus(0);
		double obj = context.getRelaxationObjective();
		
		if (status != IloCplex::Optimal &&
			status != IloCplex::OptimalInfeas) {
			return;
		}

		double bestMIP = context.getIncumbentObjective();
		cout << "Best MIP : " << bestMIP << endl;
		if (obj > bestMIP) {
			context.pruneCurrentNode();
			return;
		}

		IloNumArray v(context.getEnv());
		//IloNumArray local_lb(context.getEnv());
		int currentnodeID = context.getLongInfo(IloCplex::Callback::Context::Info::NodeUID);
		if (inst->compact_node[currentnodeID].empty()) {

			lastvisited = 0;

			cout << "Empty branch. Node ID: " << currentnodeID<< endl;
			
			//context.getLocalLB(x, local_lb);
			//for (IloInt i = 0; i < x.getSize(); ++i) {
				//if (local_lb[i] > RC_EPS) {
					//cout << x[i].getName() << " >= " << local_lb[i] << "  ";
				//}
			//}
			
			cout << endl;
			//cout << context.getLongInfo(IloCplex::Callback::Context::Info::Feasible) << endl;
			//cout << context.getLongInfo(IloCplex::Callback::Context::Info::BestBound) << endl;
			//cout << context.getLongInfo(IloCplex::Callback::Context::Info::BestSolution) << endl;
			//cout << context.getLongInfo(IloCplex::Callback::Context::Info::CandidateSource) << endl;
			//cout << context.getLongInfo(IloCplex::Callback::Context::Info::NodeDepth) << endl;
			//cout << context.getLongInfo(IloCplex::Callback::Context::Info::NodesLeft) << endl;

			//context.getRelaxationPoint(x, v);

			
			//for (IloInt i = 0; i < x.getSize(); ++i) {
				//if (x[i].getType() != IloNumVar::Float && v[i] > 1 - RC_EPS) {
					//cout << x[i].getName() << " = " << v[i] << endl;
				//}
			//}
			
		}else {
			
			lastvisited = inst->compact_node[currentnodeID].back();
		
		}
		
			context.getRelaxationPoint(x, v);
			IloInt maxVar = -1;
			IloNum maxFrac = 1.0;

			//context.getLocalLB(x, local_lb);
			
			//for (IloInt i = 0; i < x.getSize(); ++i) {
				//if (local_lb[i] > RC_EPS) {
					//cout << x[i].getName() << " >= " << local_lb[i] << "  ";
				//}
			//}
			
			//cout << endl;
			
			for (IloInt i = inst->dimension_f * lastvisited; i < inst->dimension_f * (lastvisited + 1); ++i) {

				
				//for (IloInt i = 0; i < x.getSize(); ++i) {
					//if (x[i].getType() != IloNumVar::Float && v[i] > 1 - RC_EPS) {
						//cout << x[i].getName() << " = " << v[i] << endl;
					//}
				//}
				
				if (x[i].getType() != IloNumVar::Float && find(inst->compact_node[currentnodeID].begin(), inst->compact_node[currentnodeID].end(),i - inst->dimension_f * lastvisited)== inst->compact_node[currentnodeID].end()) {
					double const frac = ::fabs(0.5 - v[i]);
					//double max = 0.0;
					if (frac < maxFrac) {
						maxFrac = v[i];
						maxVar = i;
					}
				}
			}

			if (maxVar <= 0) {
				//context.pruneCurrentNode();
				return;
			}

			IloNum minFrac = 0.001;
			if (maxFrac > minFrac) {
				CPXLONG upChild, downChild;
				double const up = ::ceil(v[maxVar]);
				double const down = ::floor(v[maxVar]);
				IloNumVar branchVar = x[maxVar];				

				// Create UP branch (branchVar >= up)
				upChild = context.makeBranch(branchVar, up, IloCplex::BranchUp, obj);
				++branches;

				// Create DOWN branch (branchVar <= down)
				downChild = context.makeBranch(branchVar, down, IloCplex::BranchDown, obj);
				++branches;

				inst->compact_node[upChild] = inst->compact_node[context.getLongInfo(IloCplex::Callback::Context::Info::NodeUID)];
				inst->compact_node[upChild].push_back(maxVar - inst->dimension_f * lastvisited);

				node_store temp = node_store(inst->compact_node[upChild], context.getRelaxationObjective());
				inst->tree[upChild] = temp;

				inst->compact_node[downChild] = inst->compact_node[context.getLongInfo(IloCplex::Callback::Context::Info::NodeUID)];
				temp = node_store(inst->compact_node[downChild], context.getRelaxationObjective());
				inst->tree[downChild] = temp;

				temp = node_store({},-1.0);
				inst->tree[context.getLongInfo(IloCplex::Callback::Context::Info::NodeUID)] = temp;

				cout << "Current node ID: " << context.getLongInfo(IloCplex::Callback::Context::Info::NodeUID) << endl;
				cout << "LP: " << context.getRelaxationObjective() << endl;
				cout << "Depth: " << depth << endl;
				size_t sizeInBytes = sizeof(inst->compact_node);
				std::cout << "Size of the map in memory: " << sizeInBytes << " bytes" << endl;
				/*for (auto v0 : inst->compact_node[currentnodeID]) {
					cout << v0 << " ";
				}*/
				cout << endl;
				//cout << "Down child: " << downChild << endl;
				cout << "Up child: " << upChild << endl;
				cout << "Variable: " << branchVar.getName() << endl;
				//cout << "Dimension: " << inst->dimension_f << endl;
				(void)downChild;
				(void)upChild;
			}
			v.end();
		
	}

	int getCalls() const { return calls; }
	int getBranches() const { return branches; }
};

void buildPathDFS(instance* inst,const std::vector<std::vector<double>>& arcs, int current, std::vector<int>& path, std::vector<bool>& visited) {
	visited[current] = true;
	path.push_back(current);

	for (int next = 0; next < inst->dimension_f; next++) {
		if (arcs[current][next] > 0.5 && !visited[next]) {
			buildPathDFS(inst,arcs, next, path, visited);
		}
	}
}

// Generic callback Dummy
/*
class BranchCallback : public IloCplex::Callback::Function {
	
public:
	BranchCallback(IloNumVarArray _x, instance* inst){
	}

	void invoke(IloCplex::Callback::Context const& context) ILO_OVERRIDE {

		//do nothing

	}
	int getCalls() const { return 0; }
	int getBranches() const { return 0; }
};
*/

void saveSolutionWithArcCoef(instance* inst, IloCplex& pathSolver, IloNumVarArray& arcs, IloArray<IloNumVarArray>& z) {

	std::string solution_file = inst->input_file + ".sol";
	std::ofstream outfile(solution_file);

	if (!outfile.is_open()) {
		std::cerr << "Error: Cannot open solution file: " << solution_file << std::endl;
		return;
	}

	outfile << std::fixed << std::setprecision(6);

	// Write objective value
	outfile << "OBJECTIVE " << pathSolver.getObjValue() << std::endl;

	// Write facility-facility arcs
	outfile << "ARCS" << std::endl;
	for (int i = 0; i < inst->dimension_f; i++) {
		for (int j = 0; j < inst->dimension_f; j++) {
			double value = pathSolver.getValue(arcs[arc_coef(i, j, inst)]);
			if (value > RC_EPS) {
				outfile << i << " " << j << " " << value << std::endl;
			}
		}
	}

	// Write facility-customer assignments
	outfile << "ASSIGNMENTS" << std::endl;
	for (int i = 0; i < inst->dimension_f; i++) {
		for (int k = 0; k < inst->dimension_c; k++) {
			double value = pathSolver.getValue(z[i][k]);
			if (value > RC_EPS) {
				outfile << i << " " << k << " " << value << std::endl;
			}
		}
	}

	outfile << "END" << std::endl;
	outfile.close();

	std::cout << "Solution saved to: " << solution_file << std::endl;
}

void master11(instance* inst) {
	
	inst->model1 = IloModel(inst->env);

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 
	
	//Objective
	IloObjective  Lateness = IloAdd(inst->model1, IloMinimize(inst->env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(inst->env,inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {
	
		z[i] = IloNumVarArray(inst->env,inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {
		
			z[i][k] = IloNumVar(inst->env, 0, 1,inttype);
			
			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());
			
			
		}
	
	}
	
	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(inst->env,inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(inst->env,inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(inst->env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			for (int ii = 0; ii < inst->dimension_f; ii++) {
			
				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());
			
			}
		}

	}
	
	//Variable s
	IloNumVarArray s = IloNumVarArray(inst->env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {
	
		string name = "s." + to_string(k);
		s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k],1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(inst->env, inst->dimension_f*inst->dimension_f, 0, 1, inttype);
	
	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			string name = "arc." + to_string(i) + "." + to_string(j);
			
			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 0.0, inttype);
			}else
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 1.0, inttype);
			
			arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}


	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {
	
		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++) 
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			
				v += inst->c_f[i][ii] * f[k][i][ii];
			
			}
		
		for (IloInt i = 0; i < inst->dimension_f; i++) {
		
			v += z[i][k] * inst->c_fc[i][k];
		
		}
	
		inst->model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}
	
	for (int k = 0; k < inst->dimension_c; k++) 
		for (int i = 0; i < inst->dimension_f; i++) {
		
			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}
			
			if (i == 0)
				inst->model1.add(v == -1 + z[0][k]);
			else
				inst->model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++) 
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				inst->model1.add(f[k][i][ii] <= arcs[arc_coef(i,ii,inst)]);
			}

		for (int k = 0; k < inst->dimension_c; k++)
			for (int i = 0; i < inst->dimension_f; i++) {

				IloExpr v(inst->env);
				for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

					//v += arcs[ii][i];
					v += arcs[arc_coef(ii,i,inst)];

				}
				inst->model1.add(z[i][k] <= v);
				v.end();

			}

		for (int k = 0; k < inst->dimension_c; k++) {
			IloExpr v(inst->env);
			for (IloInt i = 0; i < inst->dimension_f; i++) {
				v += z[i][k];
			}
			inst->model1.add(v == 1.0);
		}


		//Constraints:Capacity
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (int k = 0; k < inst->dimension_c; k++) {

				v += inst->demand[k] * z[i][k];

			}

			for (int j = 0; j < inst->dimension_f; j++) {
				v -= inst->capacity[i] * arcs[arc_coef(j, i, inst)];
			}
			inst->model1.add(v <= 0.0);
			v.end();
		}
		
		//X constraints
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr Xdeltaminus(inst->env);
			IloExpr Xdeltaplus(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
				//Xdeltaminus += arcs[ii][i];
				Xdeltaminus += arcs[arc_coef(ii,i,inst)];
				//Xdeltaplus += arcs[i][ii];
				Xdeltaplus += arcs[arc_coef(i,ii,inst)];
			}

			inst->model1.add(Xdeltaminus <= 1);
			inst->model1.add(Xdeltaminus == Xdeltaplus);
			if (i == 0) {

				inst->model1.add(Xdeltaminus == 1);

			}
			Xdeltaminus.end();
			Xdeltaplus.end();

		}
		
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
		//string filename = inst->input_file + "solution";
		//pathSolver.readMIPStarts(filename.c_str());
		pathSolver.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 14 * 1024);
		pathSolver.setParam(IloCplex::Param::TimeLimit, inst->timelimit);

		//MIPstart
		IloNumArray arc_values(inst->env, inst->dimension_f* inst->dimension_f);
		for (int i = 0; i < inst->dimension_f * inst->dimension_f; ++i) {
			arc_values[i] = IloInfinity;  // Initialize with IloInfinity (invalid value for CPLEX to ignore)
		}
		//IloNumVarArray arc_vars(inst->env);
		vector<int> tour;
		vector<bool> visited(inst->dimension_f, false);

		int currentNode = 0;
		tour.push_back(currentNode);
		visited[currentNode] = true;

		for (int i = 1; i < inst->dimension_f; i++) {
			int nearestNode = findNearestNode(inst, currentNode, visited);
			if (nearestNode == -1) break; // If no nearest node is found (should not happen in a complete graph)
			tour.push_back(nearestNode);
			arc_values[arc_coef(currentNode, nearestNode, inst)] = 1.0;
			visited[nearestNode] = true;
			currentNode = nearestNode;
		}
		arc_values[arc_coef(currentNode, 0, inst)] = 1.0;

		pathSolver.addMIPStart(arcs, arc_values, IloCplex::MIPStartAuto, "withfixedpath");
		//pathSolver.exportModel("finalMIP1.lp");
		
		//pathSolver.writeSolution(filename.c_str());
		//pathSolver.exportModel("finalMIP1.lp");
		//std::cout << "Solution status: " << pathSolver.getStatus() << std::endl;
		//std::cout << "Solution:  " << pathSolver.getObjValue() << " " << std::endl;

		//IloEnv env;

		// Read in the model
		//IloModel model(env);
		//IloCplex cplex(model);
		//IloObjective obj(env);
		//IloNumVarArray x(env);
		//IloRangeArray cons(env);
		//cplex.importModel(model, argv[a], obj, x, cons);

#ifdef impose_integrality
		
		// Register the callback function.
		BranchCallback cb(arcs,inst);
		pathSolver.use(&cb, IloCplex::Callback::Context::Id::Branching);

		//pathSolver.setParam(IloCplex::Param::MIP::Limits::Nodes, 1000);
		pathSolver.setParam(IloCplex::Param::MIP::Interval, -1);
		pathSolver.setParam(IloCplex::Param::Threads, 1);

		pathSolver.setParam(IloCplex::Param::Preprocessing::Reduce, CPX_PREREDUCE_PRIMALONLY);
		pathSolver.setParam(IloCplex::Param::Preprocessing::Linear, 0);
		
#endif impose_integrality
		//Search tree
		inst->compact_node[0] = {0};
		pathSolver.setParam(IloCplex::Param::MIP::Display, 2);
		pathSolver.setParam(IloCplex::Param::RandomSeed, 100);
		// Solve the model and report some statistics
		clock_t time_start = clock();
		bool feasible = pathSolver.solve();
		clock_t time_end = clock();

		cout << "Open nodes: " << endl;
		for (const auto& pair : inst->tree) {		
			if (pair.second.currentpath.size() > 0) {
				cout << pair.second.LP<< " : ";
				for (auto v : pair.second.currentpath) {
					cout << v << "  ";
				}
				cout << endl;
			}
		}

		cout << "Model solved, status = " << pathSolver.getStatus()
			<< endl;
		if (feasible)
			cout << "Objective value: " << pathSolver.getObjValue() << endl;
		//cout << "Callback was invoked " << cb.getCalls() << " times and created "
		//	<< cb.getBranches() << " branches" << endl;
		

		for (int i = 0; i < inst->dimension_f; i++)
			for (int j = 0; j < inst->dimension_f; j++) {

				if (pathSolver.getValue(arcs[arc_coef(i,j,inst)]) > RC_EPS) {

					inst->active_arcs.push_back(pairs(i, j));
					std::cout << "Arc: " << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << endl;
				}

			}

		vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {
			solution[i] = vector<double>(inst->dimension_c);
			for (int k = 0; k < inst->dimension_c; k++) {
				if(pathSolver.getValue(z[i][k])!=0)
					cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << endl;
				solution[i][k] = pathSolver.getValue(z[i][k]);
			}
		}

		for (int k = 0; k < inst->dimension_c; k++) {
			cout << "Customer " << k << " served by: ";
			for (int i = 0; i < inst->dimension_f; i++) {
				if (pathSolver.getValue(z[i][k]) != 0)
					cout << i << " ";
			}
			cout << endl;
		}
		
		/*for (int i = 0; i < inst->dimension_f; i++) {
			cout << t[i].getName() << " = " << pathSolver.getValue(t[i]) << endl;
		}*/		

		for (int k = 0; k < inst->dimension_c; k++) {

			for (int i = 0; i < inst->dimension_f; i++) {

				for (int ii = 0; ii < inst->dimension_f; ii++) {

					if(pathSolver.getValue(f[k][i][ii]) != 0)
					cout << f[k][i][ii].getName() << " = " << pathSolver.getValue(f[k][i][ii]) << endl;
		
				}
			}

		}

		for (int k = 0; k < inst->dimension_c; k++) {
			cout << s[k].getName() << " = " << pathSolver.getValue(s[k]) << endl;
		}
		
		double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

		write_log(inst, 11,pathSolver.getObjValue(), pathSolver.getBestObjValue(),solution_time,(int)pathSolver.getStatus());

		//generate_latex_model1(inst,solution,1);
		//inst->bestMIPfromcompact = pathSolver.getObjValue();
		inst->env.end();
}

void master1(instance* inst) {

	//log files
	string filename = inst->input_file + "tree";
	inst->branchtree.open(filename,std::ios_base::app);

	//log files
	inst->model1 = IloModel(inst->env);

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 

	//Objective
	IloObjective  Lateness = IloAdd(inst->model1, IloMinimize(inst->env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(inst->env, inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(inst->env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());


		}

	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(inst->env, inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(inst->env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());

			}
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(inst->env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		string name = "s." + to_string(k);
		s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(inst->env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			string name = "arc." + to_string(i) + "." + to_string(j);

			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 0.0, inttype);
			}
			else
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 1.0, inttype);

			arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}


	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {

		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += inst->c_f[i][ii] * f[k][i][ii];

			}

		for (IloInt i = 0; i < inst->dimension_f; i++) {

			v += z[i][k] * inst->c_fc[i][k];

		}

		inst->model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}

			if (i == 0)
				inst->model1.add(v == -1 + z[0][k]);
			else
				inst->model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++)
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				inst->model1.add(f[k][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				//v += arcs[ii][i];
				v += arcs[arc_coef(ii, i, inst)];

			}
			inst->model1.add(z[i][k] <= v);
			v.end();

		}

	for (int k = 0; k < inst->dimension_c; k++) {
		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			v += z[i][k];
		}
		inst->model1.add(v == 1.0);
	}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(inst->env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		inst->model1.add(v <= inst->capacity[i]);
		v.end();
	}



	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(inst->env);
		IloExpr Xdeltaplus(inst->env);
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		inst->model1.add(Xdeltaminus <= 1);
		inst->model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			inst->model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}

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
	pathSolver.setParam(IloCplex::Param::TimeLimit, inst->timelimit);

	pathSolver.setParam(IloCplex::Param::MIP::Display, 3);
	pathSolver.setParam(IloCplex::Param::RandomSeed, 100);
	// Solve the model and report some statistics
	clock_t time_start = clock();
	bool feasible = pathSolver.solve();
	clock_t time_end = clock();

	cout << "Open nodes: " << endl;
	for (const auto& pair : inst->tree) {
		if (pair.second.currentpath.size() > 0) {
			cout << pair.second.LP << " : ";
			for (auto v : pair.second.currentpath) {
				cout << v << "  ";
			}
			cout << endl;
		}
	}

	cout << "Model solved, status = " << pathSolver.getStatus()
		<< endl;
	inst->bestMIP = INFINITY;
	if (feasible) {
		cout << "Objective value: " << pathSolver.getObjValue() << endl;
		inst->bestMIP = pathSolver.getObjValue();
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {

				inst->active_arcs.push_back(pairs(i, j));
				std::cout << "Arc: " << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << endl;
			}

		}


	vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {
		solution[i] = vector<double>(inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << endl;
			solution[i][k] = pathSolver.getValue(z[i][k]);
		}
	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << "Customer " << k << " served by: ";
		for (int i = 0; i < inst->dimension_f; i++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << i << " ";
		}
		cout << endl;
	}

	/*for (int i = 0; i < inst->dimension_f; i++) {
		cout << t[i].getName() << " = " << pathSolver.getValue(t[i]) << endl;
	}*/


	for (int k = 0; k < inst->dimension_c; k++) {

		for (int i = 0; i < inst->dimension_f; i++) {

			for (int ii = 0; ii < inst->dimension_f; ii++) {

				if (pathSolver.getValue(f[k][i][ii]) != 0)
					cout << f[k][i][ii].getName() << " = " << pathSolver.getValue(f[k][i][ii]) << endl;

			}
		}

	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << s[k].getName() << " = " << pathSolver.getValue(s[k]) << endl;
	}

	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	write_log(inst, 1, inst->bestMIP, pathSolver.getBestObjValue(), solution_time, (int)pathSolver.getStatus());

	//generate_latex_model1(inst, solution, 1);
	inst->env.end();
}

void master13(instance* inst) {

	//log files
	string filename = inst->input_file + "tree";
	inst->branchtree.open(filename, std::ios_base::app);

	//log files
	inst->model1 = IloModel(inst->env);

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 

	//Objective
	IloObjective  Lateness = IloAdd(inst->model1, IloMinimize(inst->env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(inst->env, inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(inst->env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());


		}

	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(inst->env, inst->dimension_f);
	for (int m = 0; m < inst->dimension_f; m++) {

		f[m] = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[m][i] = IloNumVarArray(inst->env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(m) + "." + to_string(i) + "." + to_string(ii);
				f[m][i][ii].setName(name.c_str());

			}
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(inst->env, inst->dimension_c, 0.0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		string name = "s." + to_string(k);
		s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variable y
	IloNumVarArray y = IloNumVarArray(inst->env, inst->dimension_f, 0.0, 1.0, ILOINT);
	for (int i = 0; i < inst->dimension_f; i++) {

		string name = "y." + to_string(i);
		y[i].setName(name.c_str());
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(inst->env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			string name = "arc." + to_string(i) + "." + to_string(j);

			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 0.0, inttype);
			}
			else
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 1.0, inttype);

			arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}

	//big M
	double bigM = 0.0;
	for (IloInt i = 0; i < inst->dimension_f; i++)
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			bigM += inst->c_f[i][ii];
	}
	
	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {
		for (IloInt m = 1; m < inst->dimension_f; m++) {

			IloExpr v(inst->env);
			for (IloInt i = 0; i < inst->dimension_f; i++)
				for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

					v += inst->c_f[i][ii] * f[m][i][ii];

				}
			
			v += (z[m][k]-1)*(bigM + inst->c_fc[m][k]);

			v += z[m][k] * inst->c_fc[m][k];

			inst->model1.add(v <= inst->u[k] + s[k]);
			v.end();
		}

		inst->model1.add(z[0][k] * inst->c_fc[0][k] <= inst->u[k] + s[k]);
	}


	//z_ik <= y_i
	for (int k = 0; k < inst->dimension_c; k++) {
		for (IloInt i = 1; i < inst->dimension_f; i++) {
			inst->model1.add(z[i][k] <= y[i]);
		}
	}

	for (int m = 1; m < inst->dimension_f; m++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[m][ii][i];
				v -= f[m][i][ii];
			}

			if (i == 0)
				inst->model1.add(v == -y[m]);
			else if(i==m)
				inst->model1.add(v == y[m]);
			else 
				inst->model1.add(v == 0.0);
			v.end();
		}

	for (int m = 1; m < inst->dimension_f; m++)
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				inst->model1.add(f[m][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}
	inst->model1.add(y[0] == 1.0);

	for (int k = 0; k < inst->dimension_c; k++) {
		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			v += z[i][k];
		}
		inst->model1.add(v == 1.0);
	}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(inst->env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		inst->model1.add(v <= inst->capacity[i] * y[i]);
		v.end();
	}



	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(inst->env);
		IloExpr Xdeltaplus(inst->env);
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		inst->model1.add(Xdeltaminus <= 1);
		inst->model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			inst->model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}	

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
	pathSolver.setParam(IloCplex::Param::TimeLimit, inst->timelimit);

	pathSolver.setParam(IloCplex::Param::MIP::Display, 3);
	pathSolver.setParam(IloCplex::Param::RandomSeed, 100);
	pathSolver.setParam(IloCplex::Param::Threads, 1);
	//MIPstart
	IloNumArray arc_values(inst->env, inst->dimension_f* inst->dimension_f);
	for (int i = 0; i < inst->dimension_f * inst->dimension_f; ++i) {
		arc_values[i] = IloInfinity;  // Initialize with IloInfinity (invalid value for CPLEX to ignore)
	}
	//IloNumVarArray arc_vars(inst->env);
	vector<int> tour;
	vector<bool> visited(inst->dimension_f, false);

	int currentNode = 0;
	tour.push_back(currentNode);
	visited[currentNode] = true;

	for (int i = 1; i < inst->dimension_f; i++) {
		int nearestNode = findNearestNode(inst, currentNode, visited);
		if (nearestNode == -1) break; // If no nearest node is found (should not happen in a complete graph)
		tour.push_back(nearestNode);
		arc_values[arc_coef(currentNode, nearestNode, inst)] = 1.0;
		visited[nearestNode] = true;
		currentNode = nearestNode;
	}
	arc_values[arc_coef(currentNode, 0, inst)] = 1.0;

	pathSolver.addMIPStart(arcs, arc_values, IloCplex::MIPStartAuto, "withfixedpath");

	// Solve the model and report some statistics
	clock_t time_start = clock();
	bool feasible = pathSolver.solve();
	clock_t time_end = clock();

	cout << "Open nodes: " << endl;
	for (const auto& pair : inst->tree) {
		if (pair.second.currentpath.size() > 0) {
			cout << pair.second.LP << " : ";
			for (auto v : pair.second.currentpath) {
				cout << v << "  ";
			}
			cout << endl;
		}
	}

	cout << "Model solved, status = " << pathSolver.getStatus()
		<< endl;
	inst->bestMIP = INFINITY;
	if (feasible) {
		cout << "Objective value: " << pathSolver.getObjValue() << endl;
		inst->bestMIP = pathSolver.getObjValue();
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {

				inst->active_arcs.push_back(pairs(i, j));
				std::cout << "Arc: " << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << endl;
			}

		}


	vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {
		solution[i] = vector<double>(inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << endl;
			solution[i][k] = pathSolver.getValue(z[i][k]);
		}
	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << "Customer " << k << " served by: ";
		for (int i = 0; i < inst->dimension_f; i++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << i << " ";
		}
		cout << endl;
	}

	/*for (int i = 0; i < inst->dimension_f; i++) {
		cout << t[i].getName() << " = " << pathSolver.getValue(t[i]) << endl;
	}*/


	for (int m = 1; m < inst->dimension_f; m++) {

		for (int i = 0; i < inst->dimension_f; i++) {

			for (int ii = 0; ii < inst->dimension_f; ii++) {

				if (pathSolver.getValue(f[m][i][ii]) != 0)
					cout << f[m][i][ii].getName() << " = " << pathSolver.getValue(f[m][i][ii]) << endl;

			}
		}

	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << s[k].getName() << " = " << pathSolver.getValue(s[k]) << endl;
	}
	pathSolver.exportModel("Model13.lp");
	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	write_log(inst, 13, inst->bestMIP, pathSolver.getBestObjValue(), solution_time, (int)pathSolver.getStatus());

	//generate_latex_model1(inst, solution, 1);
	inst->env.end();
}

void master14(instance* inst) {

	//log files
	string filename = inst->input_file + "tree";
	inst->branchtree.open(filename, std::ios_base::app);

	//log files
	inst->model1 = IloModel(inst->env);

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 

	//Objective
	IloObjective  Lateness = IloAdd(inst->model1, IloMinimize(inst->env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(inst->env, inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(inst->env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());

		}

	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(inst->env, inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(inst->env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());

			}
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(inst->env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		string name = "s." + to_string(k);
		s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(inst->env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			string name = "arc." + to_string(i) + "." + to_string(j);

			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 0.0, inttype);
			}
			else
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 1.0, inttype);

			arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}

	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {

		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += inst->c_f[i][ii] * f[k][i][ii];

			}

		for (IloInt i = 0; i < inst->dimension_f; i++) {

			v += z[i][k] * inst->c_fc[i][k];

		}

		inst->model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}

			if (i == 0)
				inst->model1.add(v == -1 + z[0][k]);
			else
				inst->model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++)
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				inst->model1.add(f[k][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}

	//for (int k = 0; k < inst->dimension_c; k++)
	//	for (int i = 0; i < inst->dimension_f; i++) {

	//		IloExpr v(inst->env);
	//		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

	//			//v += arcs[ii][i];
	//			v += arcs[arc_coef(ii, i, inst)];

	//		}
	//		inst->model1.add(z[i][k] <= v);
	//		v.end();

	//	}

	//for (int k = 0; k < inst->dimension_c; k++) {
	//	IloExpr v(inst->env);
	//	for (IloInt i = 0; i < inst->dimension_f; i++) {
	//		v += z[i][k];
	//	}
	//	inst->model1.add(v == 1.0);
	//}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(inst->env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		
		for (int j = 0; j < inst->dimension_f; j++) {
			v -= inst->capacity[i] * arcs[arc_coef(j, i, inst)];
		}
		inst->model1.add(v <= 0.0);
		v.end();
	}



	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(inst->env);
		IloExpr Xdeltaplus(inst->env);
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		inst->model1.add(Xdeltaminus <= 1);
		inst->model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			inst->model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}

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
	pathSolver.setParam(IloCplex::Param::TimeLimit, inst->timelimit);

	pathSolver.setParam(IloCplex::Param::MIP::Display, 3);
	pathSolver.setParam(IloCplex::Param::RandomSeed, 100);
	pathSolver.setParam(IloCplex::Param::Threads, 1);
	//MIPstart
	IloNumArray arc_values(inst->env, inst->dimension_f* inst->dimension_f);
	for (int i = 0; i < inst->dimension_f * inst->dimension_f; ++i) {
		arc_values[i] = IloInfinity;  // Initialize with IloInfinity (invalid value for CPLEX to ignore)
	}
	//IloNumVarArray arc_vars(inst->env);
	vector<int> tour;
	vector<bool> visited(inst->dimension_f, false);

	int currentNode = 0;
	tour.push_back(currentNode);
	visited[currentNode] = true;

	for (int i = 1; i < inst->dimension_f; i++) {
		int nearestNode = findNearestNode(inst, currentNode, visited);
		if (nearestNode == -1) break; // If no nearest node is found (should not happen in a complete graph)
		tour.push_back(nearestNode);
		arc_values[arc_coef(currentNode, nearestNode, inst)] = 1.0;
		visited[nearestNode] = true;
		currentNode = nearestNode;
	}
	arc_values[arc_coef(currentNode, 0, inst)] = 1.0;

	pathSolver.addMIPStart(arcs, arc_values, IloCplex::MIPStartAuto, "withfixedpath");


	//root LP calculation
	//pathSolver.setParam(IloCplex::Param::MIP::Limits::Nodes, 0);   // process root only
	pathSolver.setParam(IloCplex::Param::TimeLimit, inst->timelimit);
	
	//bool ok_root = pathSolver.solve();

	double rootLB = IloInfinity;
	//if (ok_root) {
		// root relaxation bound (with root cuts etc.)
		//rootLB = pathSolver.getBestObjValue();		
	//}

	//pathSolver.setParam(IloCplex::Param::MIP::Limits::Nodes,IloIntMax);


	// Solve the model
	clock_t time_start = clock();
	bool feasible = pathSolver.solve();
	clock_t time_end = clock();

	cout << "Open nodes: " << endl;
	for (const auto& pair : inst->tree) {
		if (pair.second.currentpath.size() > 0) {
			cout << pair.second.LP << " : ";
			for (auto v : pair.second.currentpath) {
				cout << v << "  ";
			}
			cout << endl;
		}
	}

	cout << "Model solved, status = " << pathSolver.getStatus()
		<< endl;
	inst->bestMIP = INFINITY;
	if (feasible) {
		cout << "Objective value: " << pathSolver.getObjValue() << endl;
		inst->bestMIP = pathSolver.getObjValue();
		saveSolutionWithArcCoef(inst, pathSolver, arcs, z);
	}
	/*
	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {

				inst->active_arcs.push_back(pairs(i, j));
				std::cout << "Arc: " << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << endl;
			}

		}
	*/

	/*
	vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {
		solution[i] = vector<double>(inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << endl;
			solution[i][k] = pathSolver.getValue(z[i][k]);
		}
	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << "Customer " << k << " served by: ";
		for (int i = 0; i < inst->dimension_f; i++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << i << " ";
		}
		cout << endl;
	}

	//for (int i = 0; i < inst->dimension_f; i++) {
		//cout << t[i].getName() << " = " << pathSolver.getValue(t[i]) << endl;
	//}


	for (int k = 0; k < inst->dimension_c; k++) {

		for (int i = 0; i < inst->dimension_f; i++) {

			for (int ii = 0; ii < inst->dimension_f; ii++) {

				if (pathSolver.getValue(f[k][i][ii]) != 0)
					cout << f[k][i][ii].getName() << " = " << pathSolver.getValue(f[k][i][ii]) << endl;

			}
		}

	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << s[k].getName() << " = " << pathSolver.getValue(s[k]) << endl;
	}
	*/
	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	double nnodes = pathSolver.getNnodes();

	//write_log(inst, 14, inst->bestMIP, pathSolver.getBestObjValue(), solution_time, (int)pathSolver.getStatus());

	write_log_comp(inst, 14, inst->bestMIP, pathSolver.getBestObjValue(), rootLB, solution_time, (int)pathSolver.getStatus(), nnodes);

	//generate_latex_model1(inst, solution, 1);
	inst->env.end();
}

void master1_start(instance* inst) {

	//log files
	string filename = inst->input_file + "tree";
	inst->branchtree.open(filename, std::ios_base::app);

	//log files
	inst->model1 = IloModel(inst->env);

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 

	//Objective
	IloObjective  Lateness = IloAdd(inst->model1, IloMinimize(inst->env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(inst->env, inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(inst->env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());


		}

	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(inst->env, inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(inst->env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());

			}
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(inst->env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		string name = "s." + to_string(k);
		s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(inst->env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			string name = "arc." + to_string(i) + "." + to_string(j);

			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 0.0, inttype);
			}
			else
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 1.0, inttype);

			arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}


	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {

		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += inst->c_f[i][ii] * f[k][i][ii];

			}

		for (IloInt i = 0; i < inst->dimension_f; i++) {

			v += z[i][k] * inst->c_fc[i][k];

		}

		inst->model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}

			if (i == 0)
				inst->model1.add(v == -1 + z[0][k]);
			else
				inst->model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++)
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				inst->model1.add(f[k][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				//v += arcs[ii][i];
				v += arcs[arc_coef(ii, i, inst)];

			}
			inst->model1.add(z[i][k] <= v);
			v.end();

		}

	for (int k = 0; k < inst->dimension_c; k++) {
		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			v += z[i][k];
		}
		inst->model1.add(v == 1.0);
	}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(inst->env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		inst->model1.add(v <= inst->capacity[i]);
		v.end();
	}



	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(inst->env);
		IloExpr Xdeltaplus(inst->env);
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		inst->model1.add(Xdeltaminus <= 1);
		inst->model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			inst->model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}

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
	pathSolver.setParam(IloCplex::Param::TimeLimit, inst->timelimit);

	pathSolver.setParam(IloCplex::Param::MIP::Display, 3);
	pathSolver.setParam(IloCplex::Param::RandomSeed, 100);
	pathSolver.setParam(IloCplex::Param::Threads, 1);

	//MIPstart
	IloNumArray arc_values(inst->env, inst->dimension_f * inst->dimension_f);
	for (int i = 0; i < inst->dimension_f * inst->dimension_f; ++i) {
		arc_values[i] = IloInfinity;  // Initialize with IloInfinity (invalid value for CPLEX to ignore)
	}
	//IloNumVarArray arc_vars(inst->env);
	vector<int> tour;
	vector<bool> visited(inst->dimension_f, false);

	int currentNode = 0;
	tour.push_back(currentNode);
	visited[currentNode] = true;

	for (int i = 1; i < inst->dimension_f; i++) {
		int nearestNode = findNearestNode(inst, currentNode, visited);
		if (nearestNode == -1) break; // If no nearest node is found (should not happen in a complete graph)
		tour.push_back(nearestNode);
		arc_values[arc_coef(currentNode, nearestNode, inst)] = 1.0;
		visited[nearestNode] = true;
		currentNode = nearestNode;
	}
	arc_values[arc_coef(currentNode, 0, inst)] = 1.0;

	pathSolver.addMIPStart(arcs, arc_values, IloCplex::MIPStartAuto, "withfixedpath");
	
	// Solve the model and report some statistics
	clock_t time_start = clock();
	bool feasible = pathSolver.solve();
	clock_t time_end = clock();

	cout << "Open nodes: " << endl;
	for (const auto& pair : inst->tree) {
		if (pair.second.currentpath.size() > 0) {
			cout << pair.second.LP << " : ";
			for (auto v : pair.second.currentpath) {
				cout << v << "  ";
			}
			cout << endl;
		}
	}

	cout << "Model solved, status = " << pathSolver.getStatus()
		<< endl;
	inst->bestMIP = INFINITY;
	if (feasible) {
		cout << "Objective value: " << pathSolver.getObjValue() << endl;
		inst->bestMIP = pathSolver.getObjValue();
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {

				inst->active_arcs.push_back(pairs(i, j));
				std::cout << "Arc: " << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << endl;
			}

		}


	vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {
		solution[i] = vector<double>(inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << endl;
			solution[i][k] = pathSolver.getValue(z[i][k]);
		}
	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << "Customer " << k << " served by: ";
		for (int i = 0; i < inst->dimension_f; i++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << i << " ";
		}
		cout << endl;
	}

	/*for (int i = 0; i < inst->dimension_f; i++) {
		cout << t[i].getName() << " = " << pathSolver.getValue(t[i]) << endl;
	}*/


	for (int k = 0; k < inst->dimension_c; k++) {

		for (int i = 0; i < inst->dimension_f; i++) {

			for (int ii = 0; ii < inst->dimension_f; ii++) {

				if (pathSolver.getValue(f[k][i][ii]) != 0)
					cout << f[k][i][ii].getName() << " = " << pathSolver.getValue(f[k][i][ii]) << endl;

			}
		}

	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << s[k].getName() << " = " << pathSolver.getValue(s[k]) << endl;
	}

	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	write_log(inst, 1, inst->bestMIP, pathSolver.getBestObjValue(), solution_time, (int)pathSolver.getStatus());

	//generate_latex_model1(inst, solution, 1);
	inst->env.end();
}

void master1_bestbound(instance* inst) {

	//log files
	string filename = inst->input_file + "tree";
	inst->branchtree.open(filename, std::ios_base::app);

	//log files
	inst->model1 = IloModel(inst->env);

#ifdef impose_integrality
	IloNumVar::Type inttype = IloNumVar::Int;
#else
	IloNumVar::Type inttype = IloNumVar::Float;
#endif // 

	//Objective
	IloObjective  Lateness = IloAdd(inst->model1, IloMinimize(inst->env));

	//Variables
	//Variables: z
	IloArray<IloNumVarArray> z = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {

		z[i] = IloNumVarArray(inst->env, inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {

			z[i][k] = IloNumVar(inst->env, 0, 1, inttype);

			string name = "z." + to_string(i) + "." + to_string(k);
			z[i][k].setName(name.c_str());


		}

	}

	//Variable f
	IloArray<IloArray<IloNumVarArray>> f = IloArray<IloArray<IloNumVarArray>>(inst->env, inst->dimension_c);
	for (int k = 0; k < inst->dimension_c; k++) {

		f[k] = IloArray<IloNumVarArray>(inst->env, inst->dimension_f);
		for (int i = 0; i < inst->dimension_f; i++) {

			f[k][i] = IloNumVarArray(inst->env, inst->dimension_f, 0, IloInfinity, ILOFLOAT);
			for (int ii = 0; ii < inst->dimension_f; ii++) {

				string name = "f." + to_string(k) + "." + to_string(i) + "." + to_string(ii);
				f[k][i][ii].setName(name.c_str());

			}
		}

	}

	//Variable s
	IloNumVarArray s = IloNumVarArray(inst->env, inst->dimension_c, 0, IloInfinity, ILOFLOAT);
	for (int k = 0; k < inst->dimension_c; k++) {

		string name = "s." + to_string(k);
		s[k].setName(name.c_str());
		Lateness.setLinearCoef(s[k], 1.0);
	}

	//Variables: arc
	IloNumVarArray arcs = IloNumVarArray(inst->env, inst->dimension_f * inst->dimension_f, 0, 1, inttype);

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			string name = "arc." + to_string(i) + "." + to_string(j);

			if (i == j) {
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 0.0, inttype);
			}
			else
				arcs[arc_coef(i, j, inst)] = IloNumVar(inst->env, 0.0, 1.0, inttype);

			arcs[arc_coef(i, j, inst)].setName(name.c_str());

		}


	//Constraints
	for (int k = 0; k < inst->dimension_c; k++) {

		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += inst->c_f[i][ii] * f[k][i][ii];

			}

		for (IloInt i = 0; i < inst->dimension_f; i++) {

			v += z[i][k] * inst->c_fc[i][k];

		}

		inst->model1.add(v <= inst->u[k] + s[k]);
		v.end();
	}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				v += f[k][ii][i];
				v -= f[k][i][ii];
			}

			if (i == 0)
				inst->model1.add(v == -1 + z[0][k]);
			else
				inst->model1.add(v == z[i][k]);
			v.end();
		}

	for (int k = 0; k < inst->dimension_c; k++)
		for (IloInt i = 0; i < inst->dimension_f; i++)
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				inst->model1.add(f[k][i][ii] <= arcs[arc_coef(i, ii, inst)]);
			}

	for (int k = 0; k < inst->dimension_c; k++)
		for (int i = 0; i < inst->dimension_f; i++) {

			IloExpr v(inst->env);
			for (IloInt ii = 0; ii < inst->dimension_f; ii++) {

				//v += arcs[ii][i];
				v += arcs[arc_coef(ii, i, inst)];

			}
			inst->model1.add(z[i][k] <= v);
			v.end();

		}

	for (int k = 0; k < inst->dimension_c; k++) {
		IloExpr v(inst->env);
		for (IloInt i = 0; i < inst->dimension_f; i++) {
			v += z[i][k];
		}
		inst->model1.add(v == 1.0);
	}


	//Constraints:Capacity
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr v(inst->env);
		for (int k = 0; k < inst->dimension_c; k++) {

			v += inst->demand[k] * z[i][k];

		}
		inst->model1.add(v <= inst->capacity[i]);
		v.end();
	}



	//X constraints
	for (int i = 0; i < inst->dimension_f; i++) {

		IloExpr Xdeltaminus(inst->env);
		IloExpr Xdeltaplus(inst->env);
		for (IloInt ii = 0; ii < inst->dimension_f; ii++) {
			//Xdeltaminus += arcs[ii][i];
			Xdeltaminus += arcs[arc_coef(ii, i, inst)];
			//Xdeltaplus += arcs[i][ii];
			Xdeltaplus += arcs[arc_coef(i, ii, inst)];
		}

		inst->model1.add(Xdeltaminus <= 1);
		inst->model1.add(Xdeltaminus == Xdeltaplus);
		if (i == 0) {

			inst->model1.add(Xdeltaminus == 1);

		}
		Xdeltaminus.end();
		Xdeltaplus.end();

	}

	IloCplex pathSolver(inst->model1);
	pathSolver.setParam(IloCplex::Param::TimeLimit, inst->timelimit);

	pathSolver.setParam(IloCplex::Param::MIP::Display, 3);
	pathSolver.setParam(IloCplex::Param::RandomSeed, 100);
	pathSolver.setParam(IloCplex::Param::Emphasis::MIP, IloCplex::MIPEmphasisBestBound);
	
	//MIPstart
	IloNumArray arc_values(inst->env, inst->dimension_f* inst->dimension_f);
	for (int i = 0; i < inst->dimension_f * inst->dimension_f; ++i) {
		arc_values[i] = IloInfinity;  // Initialize with IloInfinity (invalid value for CPLEX to ignore)
	}
	//IloNumVarArray arc_vars(inst->env);
	vector<int> tour;
	vector<bool> visited(inst->dimension_f, false);

	int currentNode = 0;
	tour.push_back(currentNode);
	visited[currentNode] = true;

	for (int i = 1; i < inst->dimension_f; i++) {
		int nearestNode = findNearestNode(inst, currentNode, visited);
		if (nearestNode == -1) break; // If no nearest node is found (should not happen in a complete graph)
		tour.push_back(nearestNode);
		arc_values[arc_coef(currentNode, nearestNode, inst)] = 1.0;
		visited[nearestNode] = true;
		currentNode = nearestNode;
	}
	arc_values[arc_coef(currentNode, 0, inst)] = 1.0;

	pathSolver.addMIPStart(arcs, arc_values, IloCplex::MIPStartAuto, "withfixedpath");

	// Solve the model and report some statistics
	clock_t time_start = clock();
	bool feasible = pathSolver.solve();
	clock_t time_end = clock();

	cout << "Open nodes: " << endl;
	for (const auto& pair : inst->tree) {
		if (pair.second.currentpath.size() > 0) {
			cout << pair.second.LP << " : ";
			for (auto v : pair.second.currentpath) {
				cout << v << "  ";
			}
			cout << endl;
		}
	}

	cout << "Model solved, status = " << pathSolver.getStatus()
		<< endl;
	inst->bestMIP = INFINITY;
	if (feasible) {
		cout << "Objective value: " << pathSolver.getObjValue() << endl;
		inst->bestMIP = pathSolver.getObjValue();
	}

	for (int i = 0; i < inst->dimension_f; i++)
		for (int j = 0; j < inst->dimension_f; j++) {

			if (pathSolver.getValue(arcs[arc_coef(i, j, inst)]) > RC_EPS) {

				inst->active_arcs.push_back(pairs(i, j));
				std::cout << "Arc: " << arcs[arc_coef(i, j, inst)].getName() << " = " << pathSolver.getValue(arcs[arc_coef(i, j, inst)]) << endl;
			}

		}


	vector<vector<double>> solution = vector<vector<double>>(inst->dimension_f);
	for (int i = 0; i < inst->dimension_f; i++) {
		solution[i] = vector<double>(inst->dimension_c);
		for (int k = 0; k < inst->dimension_c; k++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << z[i][k].getName() << " = " << pathSolver.getValue(z[i][k]) << endl;
			solution[i][k] = pathSolver.getValue(z[i][k]);
		}
	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << "Customer " << k << " served by: ";
		for (int i = 0; i < inst->dimension_f; i++) {
			if (pathSolver.getValue(z[i][k]) != 0)
				cout << i << " ";
		}
		cout << endl;
	}

	/*for (int i = 0; i < inst->dimension_f; i++) {
		cout << t[i].getName() << " = " << pathSolver.getValue(t[i]) << endl;
	}*/


	for (int k = 0; k < inst->dimension_c; k++) {

		for (int i = 0; i < inst->dimension_f; i++) {

			for (int ii = 0; ii < inst->dimension_f; ii++) {

				if (pathSolver.getValue(f[k][i][ii]) != 0)
					cout << f[k][i][ii].getName() << " = " << pathSolver.getValue(f[k][i][ii]) << endl;

			}
		}

	}

	for (int k = 0; k < inst->dimension_c; k++) {
		cout << s[k].getName() << " = " << pathSolver.getValue(s[k]) << endl;
	}

	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;

	write_log(inst, 16, inst->bestMIP, pathSolver.getBestObjValue(), solution_time, (int)pathSolver.getStatus());

	//generate_latex_model1(inst, solution, 1);
	inst->env.end();
}