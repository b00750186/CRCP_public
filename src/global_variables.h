#pragma once

#include <ilcplex/ilocplex.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ostream>
#include <vector>
#include <string>
#include <queue>
#include <set>
#include <random>
#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <iterator>
#include <algorithm>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/r_c_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <fstream>
#include <iomanip>

#define RC_EPS 1.0e-6
#define truck_speed 30000 //speed in meters per hour
#define robot_speed 5000 //in meters per hour
#define impose_integrality

//model 2b
//time limit for generation
#define CGtimelimit 3600.0
#define addedpathlimit 10
#define CGsolvertimelimit 1000.0
//global for tree
#define timetosolvetree 500.0
#define timelimitfinalMIP 100.0

using namespace std;

static std::vector<std::vector<double>> g_min_f2c_via_f; //shortest path in terms of time from facility to customer
static std::vector<std::vector<double>> g_min_rc_f2c; //shortest path in terms of dual values from facility to customer

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloArray<IloRangeArray>> NumConstrMatrix3;
typedef IloArray<IloRangeArray> NumConstrMatrix2;
typedef IloArray<IloArray<IloNumVarArray>> NumVarMatrix3;

class path;
class path_and_score;
class subproblem_return;

class pairs{
	
public:
	pairs(int f, int s):first(f),second(s) {}
	int first;
	int second;
};

class node_store{
public:
	node_store(vector<int> cp, double LP_): currentpath(cp),LP(LP_){}
	vector<int> currentpath;
	double LP;
	node_store() = default;
};

static double currentu;
static double currentmostnegative;
static double currentdualCov;
static int numberoffacilities;
static int customernodenum;
static int currentcustomer_original_index;
static int dim_f, dim_c;

//checking for lost variables
static std::map<path,vector<int>> all_path;

typedef struct{

	ofstream branchtree;

	vector<vector<double>> c_f;
	vector<vector<double>> c_fc;
	vector<vector<double>> pos_f;
	vector<vector<double>> pos_c;
	string input_file;
	int dimension_f;
	int dimension_c;
	vector<double> capacity;
	vector<double> demand, u;
	IloEnv env, envLP;
	IloModel model1, modelLP;
	NumVarMatrix arcs, arcsLP;
	vector<double> arcsLPsol,arcsLPCGsol;
	IloObjective  Lateness;
	IloNumArray zerosArcs;
	NumConstrMatrix2 arcconstr;
	NumConstrMatrix3 arcconstr2;
	IloRangeArray  Coverage;
	IloRangeArray  Capacity;
	IloCplex CGlocalsolver,LPsolver,CGsolver;
	int path_counter;
	vector<vector<IloNumVar>> lambdatnode;
	vector<bool> capacity_constraint_initialized;
	IloNumVarArray la;
	double UB;
	double bestMIP;
	int CGcycle;
	vector<pairs> active_arcs;
	vector<path> paths;
	set<path> finalpaths;
	IloArray<IloNumVarArray> z;
	vector<int> sorted;
	vector<vector<int>> cust_to_lambda;
	vector<path> bestheur;
	IloObjective tLateness;
	IloNumArray tzerosArcs;
	NumConstrMatrix3 tarcconstr;
	NumConstrMatrix3 tarcconstr2;
	IloRangeArray  tCoverage;
	IloRangeArray  tCapacity;
	NumVarMatrix tarcs;
	vector<vector<IloNumVar>> tlambdatnode;
	vector<bool> tcapacity_constraint_initialized;
	int tpath_counter;
	vector<path> tpaths;
	IloEnv tenv;
	IloModel tempmodelCG;
	IloCplex tCGlocalsolver;
	int tCGcycle;

	//info for boost
	double most_negative_dualarcII;
	double dualCover;
	double currentu;

	//subproblem_return sr;
	vector<path> sr_return_paths;
	bool sr_return_has_negative;

	//global time counters
	double boost_time_global;
	double CG_LP_time_global;
	int pricing_number_of_cycles;
	double createCGLPtime;
	double MIPsolvetimeatnode;
	double finalMIPtime;
	double startglobal;
	double total_LS_time;

	std::map<path, int> path_quality;
	std::map<int, vector<int>> compact_node;
	std::map<int, node_store> tree;
	double bestMIPfromcompact;

	double totaldemand;
	double totalcapacity;

	double timelimit;

	bool use_LP_basedMIP;
	bool use_MIP_basedMIP;
	int numSol_limit;
	bool use_LS;
	int LS_depth_limit;
	int MIPthr;
	int LSthr;

	int boost_limit_solutions;

} instance;

class path {

	int customer;
	vector<int> visitedFacilities;
	double lateness;
public:
	path(int customer_init, vector<int> visitedFacilities_init, instance* inst) {

		customer = customer_init;
		visitedFacilities = visitedFacilities_init;

		double distance = 0.0;

		if (visitedFacilities.size() > 1) {
			int prev = visitedFacilities[0];
			int last;
			for (int i = 1; i < visitedFacilities.size(); i++) {

				last = visitedFacilities[i];
				distance += inst->c_f[prev][last];
				prev = last;
			}
			distance += inst->c_fc[last][customer];
		}
		else {
			distance = inst->c_fc[visitedFacilities[0]][customer];
		}

		lateness = max(distance - inst->u[customer], 0.0);

	}

	path() {

	}

	int const getWorkingFacility() {
		if (!visitedFacilities.empty()) {
			return visitedFacilities.back();
		}
		else {
			return -1;
		}
	}

	int const getCustomer() {
		return customer;
	}

	std::vector<int> getFacilityPath() {

		return visitedFacilities;

	}

	double getLateness() {
		return lateness;
	}

	bool operator<(const path& other) const {
		if (customer != other.customer) {
			return customer < other.customer;
		}

		return visitedFacilities < other.visitedFacilities;
	}

};

class path_and_score {
public:
	path_and_score(path p_, int score_) :p(p_), score(score_) {}
	path p;
	int score;
};

class pathvariableLess
{
public:
	bool operator()(const path_and_score a1, const path_and_score a2) const {

		return a1.score < a2.score;
		//return a1.currentXCG < a2.currentXCG;
		/*
		if (a1.currentXCG == 0.0 || a2.currentXCG == 0.0)
			return abs(a1.currentX - 0.5) > abs(a2.currentX - 0.5);
		else*/
		//return abs(a1.currentXCG - 0.5) > abs(a2.currentXCG - 0.5);

	}
};

struct nodeinfo {

	//info for boost
	double most_negative_dualarcII;
	double dualCover;
	double currentu;
	nodeinfo(double a1 = 0.0,double a2=0.0,double a3=0.0)
		:most_negative_dualarcII(a1),dualCover(a2),currentu(a3)
	{
	
	}

};

struct Lnode {
	vector<int> road;
	vector<path> p;
	double LP;
	double currentX;
	double currentXCG;
	bool requiredCG;
	vector<int> prohibited_next;
	//IloModel currentLP; 
	Lnode(vector<int> road_init, double LP_init, double currentX, double currentXCG, bool req,vector<path> p, vector<int> prohibited_next) : 
		road(road_init), LP(LP_init), currentX(currentX), currentXCG(currentXCG), requiredCG(req), p(p),prohibited_next(prohibited_next) {

	}
	Lnode() = default;
};



class subproblem_return{
public:
	bool has_negative;
	vector<path> paths;
	
	subproblem_return(bool has_negative, vector<path> ps):
	paths(ps),has_negative(has_negative){
		
	}
	subproblem_return() = default;
};

class subproblem_return_single {
public:
	double reduced_cost;
	path p;

	subproblem_return_single(double rc, path p_) {
		reduced_cost = rc;
		p = p_;
	}
	subproblem_return_single() = default;
};
class LP_return {
public:
	vector<path> p;
	double LP;
	bool status;
	LP_return(double LP, bool status, vector<path> p) :LP(LP),status(status),p(p){}
	LP_return() = default;
};



