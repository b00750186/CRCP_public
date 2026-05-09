#include "global_functions.h"

void readData(instance* inst)
{
    std::ifstream in(inst->input_file);
    string line;
    if (in) {
       
        inst->totalcapacity = 0.0;
        inst->totaldemand = 0.0;

        in >> inst->dimension_f;
        in >> inst->dimension_c;
        inst->capacity = vector<double>(inst->dimension_f);
        inst->demand = vector<double>(inst->dimension_c);
        inst->u = vector<double>(inst->dimension_c);

        inst->pos_f = vector<vector<double>>(inst->dimension_f);
        for (int i = 0; i < inst->dimension_f; i++) { 
            
            inst->pos_f[i] = vector<double>(2); 
        
        }

        inst->pos_c = vector<vector<double>>(inst->dimension_c);
        for (int j = 0; j < inst->dimension_c; j++) {

            inst->pos_c[j] = vector<double>(2);

        }

        for (int i = 0; i < inst->dimension_f; i++) {
        
            int dummy_int;
            in >> dummy_int;
            in >> inst->pos_f[i][0];
            in >> inst->pos_f[i][1];
            in >> inst->capacity[i];
            inst->totalcapacity += inst->capacity[i];
        }
        
        for (int j = 0; j < inst->dimension_c; j++) {

            int dummy_int;
            in >> dummy_int;
            in >> inst->pos_c[j][0];
            in >> inst->pos_c[j][1];
            in >> inst->demand[j];
            in >> inst->u[j];
            inst->totaldemand += inst->demand[j];
        }
        
        inst->c_f = vector<vector<double>>(inst->dimension_f, vector<double>(inst->dimension_f));
        inst->c_fc = vector<vector<double>>(inst->dimension_f, vector<double>(inst->dimension_c));

        for (int i = 0; i < inst->dimension_f; i++) {
            for (int ii = i + 1; ii < inst->dimension_f; ii++) {

                inst->c_f[i][ii] = sqrt((inst->pos_f[i][0] - inst->pos_f[ii][0]) * (inst->pos_f[i][0] - inst->pos_f[ii][0]) + (inst->pos_f[i][1] - inst->pos_f[ii][1]) * (inst->pos_f[i][1] - inst->pos_f[ii][1])) / truck_speed;
                inst->c_f[ii][i] = inst->c_f[i][ii];

            }

            for (int j = 0; j < inst->dimension_c; j++) {
            
                inst->c_fc[i][j] = sqrt((inst->pos_f[i][0] - inst->pos_c[j][0]) * (inst->pos_f[i][0] - inst->pos_c[j][0]) + (inst->pos_f[i][1] - inst->pos_c[j][1]) * (inst->pos_f[i][1] - inst->pos_c[j][1])) / robot_speed;
            
            }
        }

    }
    else {
        std::cerr << "No such file: " << inst->input_file << std::endl;
        throw(1);
    }
}

//generate latex for path model 1
void generate_latex_model2(instance* inst, IloNumArray solution) {

    ofstream compact_file;
    compact_file.open("network_solution2.tex");   
    compact_file
        << "\\documentclass[]{article}" << endl
        << "\\usepackage{tikz}" << endl
        << "\\pagenumbering{ gobble }" << endl
        << "\\begin{document}" << endl
        << "\\begin{tikzpicture}" << endl
        << "\\usetikzlibrary{shapes}" << endl
        << "\\tikzset{square/.style={regular polygon,fill=blue,regular polygon sides=4, inner sep=2}}" << endl
        << "\\tikzset{triangle/.style = {fill=blue, regular polygon, regular polygon sides=3 ,shape border rotate=120,scale = 0.7}}" << endl
        << "\\tikzset{trianglempty/.style = {fill=gray, regular polygon, regular polygon sides=3 ,shape border rotate=120,scale = 0.7}}" << endl
        << "\\tikzset{mycircle/.style={circle,scale = 0.8,fill = violet}}" << endl;

    double scale_factor = 500;
    double size = 5000;
    double point = size / scale_factor + 5;
    double shift = 2;

    compact_file << "\\clip(0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;
    
    for (int i = 0; i < inst->active_arcs.size(); i++) {
        int prev = inst->active_arcs[i].first;
        int next = inst->active_arcs[i].second;
        compact_file << "\\draw[line width=0.5mm] (" << inst->pos_f[prev][0] / scale_factor + shift << "," << inst->pos_f[prev][1] / scale_factor +  shift << ") -- (" << inst->pos_f[next][0] / scale_factor + shift << "," << inst->pos_f[next][1] / scale_factor + shift << ");" << endl;

    }

    compact_file <<
        "\\node[square, label={[font=\\small,text=blue]30:" << inst->capacity[0] << "}] at (" <<
        round((inst->pos_f[0][0]) / scale_factor * 100) / 100 + shift <<
        "," <<
        round((inst->pos_f[0][1]) / scale_factor * 100) / 100 + shift <<
        "){" << 0 << "};"
        << endl;

    for (int i = 1; i < inst->dimension_f; i++) {
    
        compact_file <<
            "\\node[trianglempty, label={[font=\\small,text=blue]30:" << inst->capacity[i] << "}] at (" <<
            round((inst->pos_f[i][0]) / scale_factor * 100) / 100 + shift <<
            "," <<
            round((inst->pos_f[i][1]) / scale_factor * 100) / 100 + shift <<
            "){" << i << "};"
            << endl;
    
    }

    for (int i = 0; i < solution.getSize(); i++) {
    
        if (solution[i] > RC_EPS) {
            
            int customer = inst->paths[i].getCustomer();
            int working_facility = inst->paths[i].getWorkingFacility();
            
            compact_file << "\\draw[green,dashed] (" << inst->pos_f[working_facility][0] / scale_factor + shift << "," << inst->pos_f[working_facility][1] / scale_factor + shift << ") -- (" << inst->pos_c[customer][0] / scale_factor + shift << "," << inst->pos_c[customer][1] / scale_factor + shift << ");" << endl;
            
            if (working_facility == 0) {
            
               
            }
            else {
                compact_file <<
                    "\\node[triangle, label={[font=\\small,text=blue]30:" << inst->capacity[working_facility] << "}] at (" <<
                    round((inst->pos_f[working_facility][0]) / scale_factor * 100) / 100 + shift <<
                    "," <<
                    round((inst->pos_f[working_facility][1]) / scale_factor * 100) / 100 + shift <<
                    "){" << working_facility << "};"
                    << endl;
            }
            
            compact_file <<
                "\\node[mycircle, label={[font=\\small,text=violet]30:" << inst->demand[customer] << "}] at (" <<
                round((inst->pos_c[customer][0]) / scale_factor * 100) / 100 + shift <<
                "," <<
                round((inst->pos_c[customer][1]) / scale_factor * 100) / 100 + shift <<
                "){" << customer << "};"
                << endl;

        }

    }
    

    compact_file << "\\draw[ line width=0.7mm ] (0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;
    compact_file
        << "\\end{tikzpicture}" << "\n"
        << "\\end{document}" << "\n"
        << endl;

    compact_file.close();

}

void generate_latex_model1(instance* inst, vector<vector<double>> solution) {

    ofstream compact_file;
    compact_file.open("network_solution1.tex");
    compact_file
        << "\\documentclass[]{article}" << endl
        << "\\usepackage{tikz}" << endl
        << "\\pagenumbering{gobble}" << endl
        << "\\begin{document}" << endl
        << "\\begin{tikzpicture}" << endl
        << "\\usetikzlibrary{shapes}" << endl
        << "\\tikzset{square/.style={regular polygon,fill=blue,regular polygon sides=4, inner sep=2}}" << endl
        << "\\tikzset{triangle/.style = {fill=blue, regular polygon, regular polygon sides=3 ,shape border rotate=120,scale = 0.7}}" << endl
        << "\\tikzset{trianglempty/.style = {fill=gray, regular polygon, regular polygon sides=3 ,shape border rotate=120,scale = 0.7}}" << endl
        << "\\tikzset{mycircle/.style={circle,scale = 0.8,fill = violet}}" << endl;

    double scale_factor = 500;
    double size = 5000;
    double point = size / scale_factor + 5;
    double shift = 2;

    compact_file << "\\clip(0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;

    for (int i = 0; i < inst->active_arcs.size(); i++) {
        int prev = inst->active_arcs[i].first;
        int next = inst->active_arcs[i].second;
        compact_file << "\\draw[line width=0.5mm] (" << inst->pos_f[prev][0] / scale_factor + shift << "," << inst->pos_f[prev][1] / scale_factor + shift << ") -- (" << inst->pos_f[next][0] / scale_factor + shift << "," << inst->pos_f[next][1] / scale_factor + shift << ");" << endl;

    }

    compact_file <<
        "\\node[square, label={[font=\\small,text=blue]30:"<<inst->capacity[0]<<"}] at (" <<
        round((inst->pos_f[0][0]) / scale_factor * 100) / 100 + shift <<
        "," <<
        round((inst->pos_f[0][1]) / scale_factor * 100) / 100 + shift <<
        "){" << 0 << "};"
        << endl;



    for (int i = 1; i < inst->dimension_f; i++) {

        compact_file <<
            "\\node[trianglempty, label={[font=\\small,text=blue]30:" << inst->capacity[i] << "}] at (" <<
            round((inst->pos_f[i][0]) / scale_factor * 100) / 100 + shift <<
            "," <<
            round((inst->pos_f[i][1]) / scale_factor * 100) / 100 + shift <<
            "){"<< i <<"};"
            << endl;

    }

    for (int k = 0; k < inst->dimension_c; k++) {
        for (int i = 0; i < inst->dimension_f; i++) {
            if (solution[i][k] > RC_EPS) {

                int customer = k;
                int working_facility = i;

                compact_file << "\\draw[green,dashed] (" << inst->pos_f[working_facility][0] / scale_factor + shift << "," << inst->pos_f[working_facility][1] / scale_factor + shift << ") -- (" << inst->pos_c[customer][0] / scale_factor + shift << "," << inst->pos_c[customer][1] / scale_factor + shift << ");" << endl;

                if (working_facility == 0) {

                    
                }
                else {
                    compact_file <<
                        "\\node[triangle, label={[font=\\small,text=blue]30:" << inst->capacity[working_facility] << "}] at (" <<
                        round((inst->pos_f[working_facility][0]) / scale_factor * 100) / 100 + shift <<
                        "," <<
                        round((inst->pos_f[working_facility][1]) / scale_factor * 100) / 100 + shift <<
                        "){" << working_facility << "};"
                        << endl;
                }


                compact_file <<
                    "\\node[mycircle, label={[font=\\small,text=violet]30:" << inst->demand[customer]<< "}] at (" <<
                    round((inst->pos_c[customer][0]) / scale_factor * 100) / 100 + shift <<
                    "," <<
                    round((inst->pos_c[customer][1]) / scale_factor * 100) / 100 + shift <<
                    "){" << customer << "};"
                    << endl;

                /*compact_file <<
                    "\\draw[line width=0.01mm,fill={rgb:red," << 100 << " ;green," << 50 << ";blue," << 100 << "}] (" <<
                    round((inst->pos_c[customer][0]) / scale_factor * 100) / 100 + shift <<
                    "," <<
                    round((inst->pos_c[customer][1]) / scale_factor * 100) / 100 + shift <<
                    ") circle (0.35);"
                    << endl;*/



            }
        }
    }


    compact_file << "\\draw[ line width=0.7mm ] (0,0) -- (" << point << ",0) -- (" << point << "," << point << ") -- (0," << point << ") -- cycle;" << endl;
    compact_file
        << "\\end{tikzpicture}" << "\n"
        << "\\end{document}" << "\n"
        << endl;

    compact_file.close();

}



void write_log(instance* inst,int method,double objective,double bestlp,double solution_time,int status) {
    std::ofstream outFile("log.txt", std::ios::app);

    if (outFile.is_open()) {
        outFile << inst->input_file << " \t " << method << " \t " << objective << " \t " << bestlp << " \t " << solution_time << " \t " << status << "\n";
    }
    else {
        std::cerr << "Error opening the file." << std::endl;
    }
}

void write_log_comp(instance* inst,
    int method,
    double objective,
    double bestlp,
    double rootlb,
    double solution_time,
    int status,
    long long nodes)   
{
    std::ofstream outFile("log.txt", std::ios::app);

    if (outFile.is_open()) {
        outFile << inst->input_file << " \t "
            << method << " \t "
            << objective << " \t "
            << bestlp << " \t "
            << rootlb << " \t "
            << solution_time << " \t "
            << status << " \t "
            << nodes << "\n";
    }
    else {
        std::cerr << "Error opening the file." << std::endl;
    }
}


void write_log1(instance* inst, int method, double objective, double solution_time, int numberofcycles,int numberoflambda,int status) {
    std::ofstream outFile("log.txt", std::ios::app);

    if (outFile.is_open()) {
        outFile << inst->input_file << " \t " << method << " \t " << objective << " \t " << solution_time << " \t " << numberofcycles << " \t " << numberoflambda << " \t " << status << "\n";
    }
    else {
        std::cerr << "Error opening the file." << std::endl;
    }
}


void write_log21(instance* inst, int method, double objective, double bestMIP, double LB, double solution_time, int numberofcycles, int numberoflambda, int status) {
    std::ofstream outFile("log.txt", std::ios::app);  // Open file in append mode

    if (outFile.is_open()) {
        // Use ostream for output, and string for file name
        outFile << inst->input_file << " \t " << method << " \t " << objective << " \t " << bestMIP << "\t" << LB << "\t" << solution_time << " \t " << numberofcycles << " \t " << numberoflambda << " \t " << status << "\n";

        // outFile automatically closed when it goes out of scope
    }
    else {
        // Handle the error, for example, print an error message
        std::cerr << "Error opening the file." << std::endl;
    }
}


int arc_coef(int i, int j, instance* inst) {
    return i * inst->dimension_f + j;
}

int findNearestNode(instance* inst, int currentNode, const vector<bool>& visited) {
    int nearestNode = -1;
    int minDistance = numeric_limits<int>::max();

    for (size_t i = 0; i < inst->dimension_f; ++i) {
        if (!visited[i] && inst->c_f[currentNode][i] < minDistance) {
            nearestNode = i;
            minDistance = inst->c_f[currentNode][i];
        }
    }

    return nearestNode;
}