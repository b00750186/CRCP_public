#include "global_functions.h"
#include "global_variables.h"
#include "boost_1.h"
#include "model2.h"
#include "model2b.h"
#include "model1.h"

void printUsage() {
    std::cout
        << "Usage: filename [-c] [-la] [-ap] [-ls] "
        << "[time=<seconds>] [nsol=<value>] [lsdepth=<value>] "
        << "[MIPthr=<value>] [LSthr=<value>] [boost_limit_solutions=<value>]\n";
}

int main(int argc, char * argv[]) {
	
    if (argc < 2) {
        printUsage();
        std::exit(EXIT_FAILURE);
    }

    instance inst;
    inst.input_file = argv[1];
    readData(&inst);

    dim_f = inst.dimension_f;
    dim_c = inst.dimension_c;

    inst.use_LP_basedMIP = false; //flag to use Heuristic 2
    inst.use_MIP_basedMIP = false; //flag to use Heuristic 1
    inst.use_LS = false; //flag to use local search
    inst.timelimit = 3600; 
    inst.numSol_limit = 1; //number of solutions to which we apply local search
    inst.MIPthr = 1; //skip factor for Heuristic 1 or 2
    inst.LSthr = 1; //skip factor for local search, local search cannot be performed more frequently than Heuristic 1 or 2
    inst.LS_depth_limit = 0; //BnP search tree depth from which local search is allowed
    inst.boost_limit_solutions = 10; //max number of paths with negative reduced cost for each customer we take from pricing subproblem
    

    // Parse command-line arguments
    for (int i = 2; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg == "-h") {
            printUsage();
            return 0;
        }

        if (arg == "-naive") {
            cout << "Naive solution:" << inst.input_file << "|" << dummy(&inst) << endl;
            return 0;
        }

        if (arg == "-c") {
            master14(&inst); //BnP initialized with heuristic
            return 0;
        }

        if (arg == "-la") {
            inst.use_MIP_basedMIP = true; 
        }
        else if (arg == "-ap") {
            inst.use_LP_basedMIP = true;           
        }
        else if (arg == "-ls") {
            inst.use_LS = true;
        }
        else if (arg.find("time=") == 0) {
            try {
                inst.timelimit = std::stoi(arg.substr(5));
            }
            catch (...) {
                std::cerr << "Invalid format for time.\n";
                printUsage();
                std::exit(EXIT_FAILURE);
            }
        }
        else if (arg.find("nsol=") == 0) {
            try {
                inst.numSol_limit = std::stoi(arg.substr(5));
            }
            catch (...) {
                std::cerr << "Invalid format for number_of_sol.\n";
                printUsage();
                std::exit(EXIT_FAILURE);
            }
        }
        else if (arg.find("lsdepth=") == 0) {
            try {
                inst.LS_depth_limit = std::stoi(arg.substr(8));
            }
            catch (...) {
                std::cerr << "Invalid format for LS depth limit.\n";
                printUsage();
                std::exit(EXIT_FAILURE);
            }
        }
        else if (arg.find("MIPthr=") == 0) {
            try {
                inst.MIPthr = std::stoi(arg.substr(7)); 
            }
            catch (...) {
                std::cerr << "Invalid format for MIP throttle.\n";
                printUsage();
                std::exit(EXIT_FAILURE);
            }
        }
        else if (arg.find("LSthr=") == 0) {
            try {
                inst.LSthr = std::stoi(arg.substr(6)); 
            }
            catch (...) {
                std::cerr << "Invalid format for LS throttle.\n";
                printUsage();
                std::exit(EXIT_FAILURE);
            }
        }
        else if (arg.rfind("boost_limit_solutions=", 0) == 0) { 
            try {
                const auto key = std::string("boost_limit_solutions=");
                inst.boost_limit_solutions = std::stoi(arg.substr(key.size()));
            }
            catch (...) {
                std::cerr << "Invalid format for boost_limit_solutions.\n";
                printUsage(); std::exit(EXIT_FAILURE);
            }
        }
        else {
            std::cerr << "Unknown or malformed argument: " << arg << std::endl;
            printUsage();
            std::exit(EXIT_FAILURE);
        }
    }

    precompute_min_f2c_via_f(&inst);

    // Output the results of argument parsing
    std::cout << "Flags set:\n";
    std::cout << "  Use lambda in {0,1}: " << (inst.use_MIP_basedMIP ? "Yes" : "No") << "\n";
    std::cout << "  Use x in A_p: " << (inst.use_LP_basedMIP ? "Yes" : "No") << "\n";
    std::cout << "  Use local search: " << (inst.use_LS ? "Yes" : "No") << "\n";
    std::cout << "Parameters:\n";
    std::cout << "  time: " << inst.timelimit << "\n";
    std::cout << "  number_of_sol: " << inst.numSol_limit << "\n";
    std::cout << "  LS depth limit: " << inst.LS_depth_limit << "\n";
    std::cout << "  MIP throttle: " << inst.MIPthr << "\n";
    std::cout << "  LS throttle: " << inst.LSthr << "\n";
    master2c(&inst); // with branching
    //master2(&inst); //no branching

    return 0;
}
