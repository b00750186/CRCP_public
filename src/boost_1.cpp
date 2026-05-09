#include "boost_1.h"
using namespace boost;

inline double round6(double x) { return std::round(x * 1e6) / 1e6; }

void precompute_min_f2c_via_f(const instance* inst) { //precompute shortest path from a facility to a customer
    const int F = inst->dimension_f;
    const int C = inst->dimension_c;

    const double INF = std::numeric_limits<double>::infinity();

    g_min_f2c_via_f.assign(F, std::vector<double>(C, INF));

    for (int j = 0; j < F; ++j) {
        for (int k = 0; k < C; ++k) {
            double best = INF;
            for (int i = 0; i < F; ++i) {
                double a = inst->c_f[j][i];   // facility j -> facility i
                double b = inst->c_fc[i][k];  // facility i -> customer k
                if (std::isinf(a) || std::isinf(b)) continue; // unreachable leg
                double cand = a + b;
                if (cand < best) best = cand;
            }
            g_min_f2c_via_f[j][k] = best; // stays INF if no 2-hop path exists
        }
    }

    /*
    if (g_min_rc_f2c.size() != static_cast<size_t>(F) ||
        (F > 0 && g_min_rc_f2c[0].size() != static_cast<size_t>(C))) {
        g_min_rc_f2c.assign(F, std::vector<double>(C, INF));
    }*/
}

struct Candidate {
    double score; // the reduced-cost value
    int idx;      // index in opt_solutions
};

struct CandidateCmp {
    // min-heap by score
    bool operator()(const Candidate& a, const Candidate& b) const {
        return a.score > b.score;
    }
};

//boost structures
struct Graph_Node
{
    Graph_Node(int n = 0, double t_ = 0, double t = 0)
        : num(n), t_bar(t_), time(t){}
    int num;
    double t_bar;
    double time;
};

struct Graph_Arc {

    Graph_Arc(int n = 0, double t_ = 0, double t = 0)
        : num(n), t_bar(t_), time(t){}
    int num;
    double t_bar;
    double time;
};

typedef adjacency_list<vecS, vecS, directedS, Graph_Node, Graph_Arc> Subproblem_Graph;

//boost data structures for shortest path problem with time windows (spptw)

struct subproblem_graph_cont
{
    subproblem_graph_cont(double t_ = 0, double t = 0) : t_bar(t_), time(t){
    
    }
    subproblem_graph_cont& operator=(const subproblem_graph_cont& other)
    {
        if (this == &other)
            return *this;
        this->~subproblem_graph_cont();
        new(this) subproblem_graph_cont(other);
        return *this;
    }
    double t_bar;
    double time;
    int destination;
};


bool operator==(const subproblem_graph_cont& arg1, const subproblem_graph_cont& arg2)
{
    return (arg1.t_bar == arg2.t_bar && arg1.time == arg2.time);
}


bool operator<(const subproblem_graph_cont& arg1, const subproblem_graph_cont& arg2)
{
    if ((arg1.t_bar <= arg2.t_bar && arg1.time < arg2.time) || (arg1.t_bar < arg2.t_bar && arg1.time <= arg2.time))
        return true;
    if (arg1.time > currentu && arg2.time > currentu)
        return (arg1.time + arg1.t_bar < arg2.time + arg2.t_bar);
    return false;
}


// label extention

class label_extention
{
public:
    inline bool operator()(const Subproblem_Graph& g,
        subproblem_graph_cont& new_cont,
        const subproblem_graph_cont& old_cont,
        graph_traits
        <Subproblem_Graph>::edge_descriptor ed) const
    {
        const Graph_Arc& arc_prop = get(edge_bundle, g)[ed];
        const Graph_Node& vert_prop = get(vertex_bundle, g)[target(ed, g)];
        new_cont.t_bar = old_cont.t_bar + arc_prop.t_bar;
        new_cont.time = old_cont.time + arc_prop.time;
        new_cont.destination = vert_prop.num;
        return true;
    }
};


// DominanceFunction model
class dominance_rules_old
{
public:
    inline bool operator()(const subproblem_graph_cont& arg1, const subproblem_graph_cont& arg2) const
     {
       
        if (max(0.0, arg2.time - currentu) + arg2.t_bar - currentdualCov + currentmostnegative > RC_EPS)
            return true;

        if (arg1.time > currentu && arg2.time > currentu)
            return (arg1.time + arg1.t_bar <= arg2.time + arg2.t_bar);

        else
            return (arg1.time <= arg2.time && arg1.t_bar <= arg2.t_bar);

    }
};

class dominance_rules
{
public:
    inline bool operator()(const subproblem_graph_cont& arg1, const subproblem_graph_cont& arg2) const
    {
        const int v2 = arg2.destination;
        const int v1 = arg1.destination;

        double shortestfromheretocustomer = 0;

        if (v2 < dim_f && v2 >= 0 && v1 == v2) {
            shortestfromheretocustomer = g_min_f2c_via_f[v2][currentcustomer_original_index];
        }

        if (v2 < dim_f && v2 >= 0 && v1 == v2) {
            
            if (max(0.0, arg2.time + shortestfromheretocustomer - currentu) + arg2.t_bar - currentdualCov + currentmostnegative > -RC_EPS)
                return true;        
        }        
   
        if (arg1.time + shortestfromheretocustomer > currentu && arg2.time + shortestfromheretocustomer > currentu)
            return (arg1.time + arg1.t_bar <= arg2.time + arg2.t_bar);

        else
            return (arg1.time <= arg2.time && arg1.t_bar <= arg2.t_bar);

    }
};

class dominance_rules_late
{
public:
    inline bool operator()(const subproblem_graph_cont& arg1, const subproblem_graph_cont& arg2) const
    {

        if (max(0.0, arg2.time - currentu) + arg2.t_bar - currentdualCov + currentmostnegative > -RC_EPS)
            return true;

        return (arg1.time + arg1.t_bar <= arg2.time + arg2.t_bar);
        
    }
};

template <class Graph>
double shortest_time_value_dijkstra(
    const Graph& g,
    typename boost::graph_traits<Graph>::vertex_descriptor s,
    typename boost::graph_traits<Graph>::vertex_descriptor t)
{
    using boost::get;
    const std::size_t n = num_vertices(g);
    auto index = get(&Graph_Node::num, g);
    auto weight = get(&Graph_Arc::time, g);

    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pred(n);
    auto dist_map = boost::make_iterator_property_map(dist.begin(), index);
    auto pred_map = boost::make_iterator_property_map(pred.begin(), index);

    boost::dijkstra_shortest_paths(
        g, s,
        boost::weight_map(weight)
        .distance_map(dist_map)
        .predecessor_map(pred_map)
        .vertex_index_map(index)
    );
    return get(dist_map, t);
}

subproblem_return_single boost_shortest(int customer, instance* inst, IloNumArray dualCapacity, IloArray<NumMatrix> dualArcs, IloNumArray dualCoverage)
{
    Subproblem_Graph gg;

    double u = inst->u[customer];
    int counter = 0;

    for (int i = 0; i < inst->dimension_f; i++) {
        add_vertex(Graph_Node(counter, 0, 0), gg);
        counter++;
    }

    add_vertex(Graph_Node(counter, 0, 0), gg);

    //inter facility
    counter = 0;
    for (int i = 0; i < inst->dimension_f - 1; i++)
        for (int ii = i + 1; ii < inst->dimension_f; ii++) {
            add_edge(i, ii, Graph_Arc(counter, dualArcs[customer][i][ii], inst->c_f[i][ii]), gg);
            counter++;
            add_edge(ii, i, Graph_Arc(counter, dualArcs[customer][ii][i], inst->c_f[ii][i]), gg);
            counter++;
        }

    for (int i = 0; i < inst->dimension_f; i++)
        if(inst->capacity[i] >= inst->demand[customer])
    {
        add_edge(i, inst->dimension_f, Graph_Arc(counter, -dualCapacity[i] * inst->demand[customer], inst->c_fc[i][customer]), gg);
        counter++;
    }
    // spp without resource constraints
    graph_traits<Subproblem_Graph>::vertex_descriptor s = 0;
    graph_traits<Subproblem_Graph>::vertex_descriptor t = inst->dimension_f;
    customernodenum = inst->dimension_f;

    std::vector<std::vector<graph_traits<Subproblem_Graph>::edge_descriptor>> opt_solutions;
    std::vector<subproblem_graph_cont> pareto_opt_rcs_no_rc;

    currentu = inst->u[customer];
    currentdualCov = dualCoverage[customer];

    numberoffacilities = inst->dimension_f;
    currentmostnegative = 0;

    r_c_shortest_paths
    (gg,
        get(&Graph_Node::num, gg),
        get(&Graph_Arc::num, gg),
        s,
        t,
        opt_solutions,
        pareto_opt_rcs_no_rc,
        subproblem_graph_cont(),
        label_extention(),
        dominance_rules_old(),
        std::allocator
        <r_c_shortest_paths_label
        <Subproblem_Graph, subproblem_graph_cont> >(),
        default_r_c_shortest_paths_visitor());


    vector<int> p_;
    //std::cout << static_cast<int>(opt_solutions.size()) << std::endl;

    //optimal at the last node
    int optim = 0;
    double incumbent = INFINITY;
    //cout << "Shortest paths: " << customer << endl;
    for (int i = 0; i < static_cast<int>(opt_solutions.size()); ++i) {
        /*cout << "t: " << pareto_opt_rcs_no_rc[i].time << " bar t: " << pareto_opt_rcs_no_rc[i].t_bar << endl;
        for (int j = static_cast<int>(opt_solutions[i].size()) - 1; j >= 0; --j) {
            std::cout << source(opt_solutions[i][j], gg) << std::endl;
        }*/

        if (max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar < incumbent) {
            optim = i;
            incumbent = max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar;
        }
    }
    //cout << "Customer: " << customer << " Reduced cost: " << incumbent - inst->localsolver.getDual(inst->Coverage[customer]) << endl;
    //cout << "Path:" << endl;
    for (int j = static_cast<int>(opt_solutions[optim].size()) - 1; j >= 0; --j) {
        //std::cout << source(opt_solutions[optim][j], gg) << std::endl;
        p_.push_back(source(opt_solutions[optim][j], gg));
    }
    // std::cout << "Time: " << pareto_opt_rcs_no_rc[optim].time << std::endl;
     //std::cout << "Time_: " << pareto_opt_rcs_no_rc[optim].t_bar << std::endl;

    path p = path(customer, p_, inst);
    subproblem_return_single sr = subproblem_return_single(max(pareto_opt_rcs_no_rc[optim].time - u, 0.0) + pareto_opt_rcs_no_rc[optim].t_bar - dualCoverage[customer], p);
    return sr;
}

void boost_shortest11(instance * inst, vector<int> &road, const std::vector<char>& inRoad, const std::vector<char>& avoid)
{
	bool has_negative = false;
	inst->sr_return_paths.clear();

	Subproblem_Graph gg;

	int counter = 0;
	nodeinfo info;

	//facility nodes
	for (int i = 0; i < inst->dimension_f; i++) {
		add_vertex(Graph_Node(counter, 0, 0), gg);
		counter++;
	}
	//customer nodes start from index inst->dimension_f
	for (int i = 0; i < inst->dimension_c; i++) {
		add_vertex(Graph_Node(counter, 0, 0), gg);
		counter++;
	}

    const int F = inst->dimension_f;
    

    std::vector<double> capDual(F, 0.0);
    for (int i = 0; i < F; ++i) capDual[i] = round6(inst->CGsolver.getDual(inst->Capacity[i]));

    std::vector<std::vector<double>> arcDual(F, std::vector<double>(F, 0.0));
    for (int i = 0; i < F; ++i)
        for (int j = 0; j < F; ++j)
            if (i != j)
                arcDual[i][j] = round6(inst->CGsolver.getDual(inst->arcconstr[i][j]));

    std::vector<double> piarcII(F, 0.0);
    for (int i = 0; i < F; ++i) {
        if (i == 0) { piarcII[i] = 0.0; continue; }
        double s = 0.0;
        for (int ii = 0; ii < F; ++ii) {
            s += inst->CGsolver.getDual(inst->arcconstr2[0][i][ii])
                + inst->CGsolver.getDual(inst->arcconstr2[1][ii][i]);
        }
        piarcII[i] = round6(s);
    }

    

    //inter facility
    
	counter = 0;
	int prev = 0;
	for (auto v : road) if (v != 0) {
		add_edge(prev, v, Graph_Arc(counter, arcDual[prev][v], inst->c_f[prev][v]), gg);
		prev = v;
		counter++;
	}

	for (int ii = 0; ii < inst->dimension_f; ii++) {
		if (!inRoad[ii] && !avoid[ii]) {
			add_edge(road.back(), ii, Graph_Arc(counter, arcDual[road.back()][ii], inst->c_f[road.back()][ii]), gg);
			counter++;
		}
	}

    std::vector<int> freeF; freeF.reserve(F);
    for (int i = 0; i < F; ++i) if (!inRoad[i]) freeF.push_back(i);
    for (size_t a = 0; a + 1 < freeF.size(); ++a) {
        const int i = freeF[a];
        for (size_t b = a + 1; b < freeF.size(); ++b) {
            const int j = freeF[b];
            add_edge(i, j, Graph_Arc(counter, arcDual[i][j], inst->c_f[i][j]), gg);
            counter++;
            add_edge(j, i, Graph_Arc(counter, arcDual[j][i], inst->c_f[j][i]), gg);
            counter++;
        }
    }

	//for (int i = 0; i < inst->dimension_f - 1; i++)
	//	if (std::find(road.begin(), road.end(), i) == road.end()) {
	//		for (int ii = i + 1; ii < inst->dimension_f; ii++)
	//		{
	//			if (std::find(road.begin(), road.end(), ii) == road.end()) {
 //                   double dual = std::round(inst->CGsolver.getDual(inst->arcconstr[i][ii])*1000000)/1000000;
	//				add_edge(i, ii, Graph_Arc(counter, dual, inst->c_f[i][ii]), gg);
	//				counter++;
 //                   dual = std::round(inst->CGsolver.getDual(inst->arcconstr[ii][i])*1000000)/1000000;
 //                   add_edge(ii, i, Graph_Arc(counter, dual, inst->c_f[ii][i]), gg);
	//				counter++;
	//			}
	//		}
	//	}

	//to customers
	//for (int i = 0; i < inst->dimension_f; i++)
	//    for (int j = 0; j < inst->dimension_c; j++) {
	//        add_edge(i, j + inst->dimension_f, Graph_Arc(counter, -dualCapacity[i] * inst->demand[customer], inst->c_fc[i][j]), gg);
	//        counter++;
	//    }

	vector<double> mostnegative = vector<double>(inst->dimension_c, 0.0);
	for (int customer = 0; customer < inst->dimension_c; customer++)
		for (int i = 0; i < inst->dimension_f; i++)
		{
			if (inst->capacity[i] >= inst->demand[customer]) {
				//arc constraints II duals
                const double w_t = -capDual[i] * inst->demand[customer] + piarcII[i];
                add_edge(i, F + customer, Graph_Arc(counter, w_t, inst->c_fc[i][customer]), gg);
                counter++;
                if (w_t < mostnegative[customer]) mostnegative[customer] = w_t;
                
                //inst->dimension_f + customer - index of customer node
			}
		}
    // Print all vertices and edges
    //std::cout << "Vertices:" << std::endl;
    //auto vp = vertices(gg);
    //for (auto it = vp.first; it != vp.second; ++it) {
    //    std::cout << *it << std::endl;
    //}

   /* std::cout << "Edges:" << std::endl;
    auto ep = edges(gg);
    for (auto it = ep.first; it != ep.second; ++it) {
        auto edge = *it;
        std::cout << source(edge, gg) << " -> " << target(edge, gg)
            << " number: " << gg[edge].num << " t:" << gg[edge].time << " t_:" << gg[edge].t_bar << std::endl;
    }*/
    //END pring all vertices and edges

    graph_traits<Subproblem_Graph>::vertex_descriptor s = 0;

      

    for (int customer = 0; customer < inst->dimension_c; customer++) {
        // spp without resource constraints		

        //end point - the customer
        graph_traits<Subproblem_Graph>::vertex_descriptor t = inst->dimension_f + customer;
        customernodenum = inst->dimension_f + customer;

        std::vector<std::vector<graph_traits<Subproblem_Graph>::edge_descriptor>> opt_solutions;
        std::vector<subproblem_graph_cont> pareto_opt_rcs_no_rc;
        currentcustomer_original_index = customer;
        currentu = inst->u[customer];
        currentdualCov = inst->CGsolver.getDual(inst->Coverage[customer]);
        currentmostnegative = mostnegative[customer];
        numberoffacilities = inst->dimension_f;
        //cout << "customer: " << customer << endl;        

        //const double sp_time = shortest_time_value_dijkstra(gg, s, t);
        //constexpr double EPS = 1e-9;

        //const bool use_late = !std::isfinite(sp_time) || (sp_time > currentu + EPS);

        r_c_shortest_paths(
            gg,
            get(&Graph_Node::num, gg),
            get(&Graph_Arc::num, gg),
            s,
            t,
            opt_solutions,
            pareto_opt_rcs_no_rc,
            subproblem_graph_cont(),
            label_extention(),
            dominance_rules(),
            std::allocator< r_c_shortest_paths_label<Subproblem_Graph, subproblem_graph_cont> >(),
            default_r_c_shortest_paths_visitor()
        );


        //std::cout << static_cast<int>(opt_solutions.size()) << std::endl;
        vector<vector<int>> negative_paths;

        //optimal at the last node
        /*int optim = 0.0;
        double incumbent = INFINITY;
        for (int i = 0; i < static_cast<int>(opt_solutions.size()); ++i) {
            if (max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar < incumbent) {
                optim = i;
                incumbent = max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar;
            }
        }
        cout << "Customer: " << customer << " Reduced cost: " << incumbent - inst->localsolver.getDual(inst->Coverage[customer]) << endl;
        for (int j = static_cast<int>(opt_solutions[optim].size()) - 1; j >= 0; --j) {
            std::cout << source(opt_solutions[optim][j], gg) << std::endl;
            p_.push_back(source(opt_solutions[optim][j], gg));
        }

        std::cout << "Time: " << pareto_opt_rcs_no_rc[optim].time << std::endl;
        std::cout << "Time_: " << pareto_opt_rcs_no_rc[optim].t_bar << std::endl;


        */

        ////print all
        //for (int i = 0; i < static_cast<int>(opt_solutions.size()); ++i) {
        //        for (int j = static_cast<int>(opt_solutions[i].size()) - 1; j >= 0; --j) {
        //            std::cout << source(opt_solutions[i][j], gg) << " - ";
        //        }
        //        cout << endl;
        //    }       

        
        
        std::priority_queue<Candidate, std::vector<Candidate>, CandidateCmp> pq;
        double dualcov = inst->CGsolver.getDual(inst->Coverage[customer]);
        for (int i = 0; i < static_cast<int>(opt_solutions.size()); ++i) {
            const double rc =
                std::max(pareto_opt_rcs_no_rc[i].time - currentu, 0.0) +
                pareto_opt_rcs_no_rc[i].t_bar - dualcov;

            if (rc < -RC_EPS) {
                has_negative = true;
                pq.push(Candidate{ rc, i });
            }

        }

        const int take = std::min<int>(static_cast<int>(pq.size()), inst->boost_limit_solutions);

        for (int k = 0; k < take; ++k) {
            const int idx = pq.top().idx;
            pq.pop();

            std::vector<int> p_;
            p_.reserve(opt_solutions[idx].size());
            for (int j = static_cast<int>(opt_solutions[idx].size()) - 1; j >= 0; --j) {
                p_.push_back(source(opt_solutions[idx][j], gg));
            }
            inst->sr_return_paths.push_back(path(customer, p_, inst));
        }    

    }
	//subproblem_return sr = subproblem_return(has_negative, paths);   
	inst->sr_return_has_negative = has_negative;
	//return sr;
}



//subproblem_return tboost_shortest11(instance* inst)//, IloNumArray dualCapacity, IloArray<NumMatrix> dualArcs, IloNumArray dualCoverage, IloArray<NumMatrix> dualArcs2)
//{
//    bool has_negative = false;
//    vector<path> paths;
//    for (int customer = 0; customer < inst->dimension_c; customer++) {
//        Subproblem_Graph gg;
//
//        double u = inst->u[customer];
//        int counter = 0;
//        nodeinfo info;
//        for (int i = 0; i < inst->dimension_f; i++) {
//            add_vertex(Graph_Node(counter, 0, 0), gg);
//            counter++;
//        }
//
//        add_vertex(Graph_Node(counter, 0, 0), gg);
//
//        /*
//        for (int j = 0; j < inst->dimension_c; j++) {
//
//            add_vertex(Graph_Node(counter, 0, 0, u), gg);
//            counter++;
//        }
//        */
//        //inter facility
//        counter = 0;
//        for (int i = 0; i < inst->dimension_f - 1; i++)
//            for (int ii = i + 1; ii < inst->dimension_f; ii++) {
//
//                //add_edge(i, ii, Graph_Arc(counter, dualArcs[customer][i][ii], inst->c_f[i][ii]), gg);
//                double temp = inst->tCGsolver.getDual(inst->tarcconstr[customer][i][ii]);
//                add_edge(i, ii, Graph_Arc(counter, inst->tCGsolver.getDual(inst->tarcconstr[customer][i][ii]), inst->c_f[i][ii]), gg);
//                counter++;
//                add_edge(ii, i, Graph_Arc(counter, inst->tCGsolver.getDual(inst->tarcconstr[customer][ii][i]), inst->c_f[ii][i]), gg);
//                counter++;
//            }
//
//        //to customers
//        //for (int i = 0; i < inst->dimension_f; i++)
//        //    for (int j = 0; j < inst->dimension_c; j++) {
//        //        add_edge(i, j + inst->dimension_f, Graph_Arc(counter, -dualCapacity[i] * inst->demand[customer], inst->c_fc[i][j]), gg);
//        //        counter++;
//        //    }
//
//        double mostnegative = 0;
//        for (int i = 0; i < inst->dimension_f; i++)
//        {
//            //arc constraints II duals
//            double piarcII = 0.0;
//            if (i != 0) {
//                for (int ii = 0; ii < inst->dimension_f; ii++) {
//
//                    piarcII += inst->tCGsolver.getDual(inst->tarcconstr2[0][i][ii]) + inst->tCGsolver.getDual(inst->tarcconstr2[1][ii][i]);
//
//                }
//            }
//            if (piarcII < mostnegative)
//                mostnegative = piarcII;
//            if (inst->tcapacity_constraint_initialized[i])
//                add_edge(i, inst->dimension_f, Graph_Arc(counter, -inst->tCGsolver.getDual(inst->tCapacity[i]) * inst->demand[customer] + piarcII, inst->c_fc[i][customer]), gg);
//            else
//                add_edge(i, inst->dimension_f, Graph_Arc(counter, 0 + piarcII, inst->c_fc[i][customer]), gg);
//            counter++;
//        }
//        // spp without resource constraints
//        graph_traits<Subproblem_Graph>::vertex_descriptor s = 0;
//        graph_traits<Subproblem_Graph>::vertex_descriptor t = inst->dimension_f;
//        customernodenum = inst->dimension_f;
//
//        std::vector<std::vector<graph_traits<Subproblem_Graph>::edge_descriptor>> opt_solutions;
//        std::vector<subproblem_graph_cont> pareto_opt_rcs_no_rc;
//
//        currentu = inst->u[customer];
//        currentdualCov = inst->tCGsolver.getDual(inst->tCoverage[customer]);
//        currentmostnegative = mostnegative;
//        numberoffacilities = inst->dimension_f;
//        already_found_negative_path = false;
//        r_c_shortest_paths
//        (gg,
//            get(&Graph_Node::num, gg),
//            get(&Graph_Arc::num, gg),
//            s,
//            t,
//            opt_solutions,
//            pareto_opt_rcs_no_rc,
//            subproblem_graph_cont(),
//            label_extention(),
//            dominance_rules(),
//            std::allocator
//            <r_c_shortest_paths_label
//            <Subproblem_Graph, subproblem_graph_cont> >(),
//            default_r_c_shortest_paths_visitor());
//
//
//        //std::cout << static_cast<int>(opt_solutions.size()) << std::endl;
//        vector<vector<int>> negative_paths;
//
//        //optimal at the last node
//        /*int optim = 0.0;
//        double incumbent = INFINITY;
//        for (int i = 0; i < static_cast<int>(opt_solutions.size()); ++i) {
//            if (max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar < incumbent) {
//                optim = i;
//                incumbent = max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar;
//            }
//        }
//        cout << "Customer: " << customer << " Reduced cost: " << incumbent - inst->localsolver.getDual(inst->Coverage[customer]) << endl;
//        for (int j = static_cast<int>(opt_solutions[optim].size()) - 1; j >= 0; --j) {
//            std::cout << source(opt_solutions[optim][j], gg) << std::endl;
//            p_.push_back(source(opt_solutions[optim][j], gg));
//        }
//
//        std::cout << "Time: " << pareto_opt_rcs_no_rc[optim].time << std::endl;
//        std::cout << "Time_: " << pareto_opt_rcs_no_rc[optim].t_bar << std::endl;
//
//
//        */
//
//        ////print all
//        //for (int i = 0; i < static_cast<int>(opt_solutions.size()); ++i) {
//        //        for (int j = static_cast<int>(opt_solutions[i].size()) - 1; j >= 0; --j) {
//        //            std::cout << source(opt_solutions[i][j], gg) << " - ";
//        //        }
//        //        cout << endl;
//        //    }       
//
//        int n_considered = 0;
//        for (int i = 0; i < static_cast<int>(opt_solutions.size()); ++i) {
//            //cout << "Reduced cost: " << max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar - inst->localsolver.getDual(inst->Coverage[customer]) << endl;
//            if (max(pareto_opt_rcs_no_rc[i].time - u, 0.0) + pareto_opt_rcs_no_rc[i].t_bar - inst->tCGsolver.getDual(inst->tCoverage[customer]) < -RC_EPS) {
//                has_negative = true;
//                n_considered++;
//                vector<int> p_;
//                for (int j = static_cast<int>(opt_solutions[i].size()) - 1; j >= 0; --j) {
//                    //std::cout << source(opt_solutions[i][j], gg) << std::endl;
//                    p_.push_back(source(opt_solutions[i][j], gg));
//                }
//                path p = path(customer, p_, inst);
//                paths.push_back(p);
//            }
//            if (n_considered > addedpathlimit)
//                break;
//        }
//    }
//    subproblem_return sr = subproblem_return(has_negative, paths);
//    return sr;
//}
