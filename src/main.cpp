#include <iostream>

#include "input.h"
#include "core/algorithms.h"
#include "core/deployment.h"

using namespace std;

void run_experiment(input &par) {
    if (par.experiment != 0) {
        par = load_parameters(par);
    }

    print_parameters(par);

    // Deployment creation with respect to the input parameters
    deployment dep(par);

    // Algorithms creation and invocation
    algorithms alg(&dep);
    alg.algorithm_1();
    alg.algorithm_2();
}

int main(int argc, char **argv) {

    input par;

    if (argc >= 3) {
        par.experiment = atoi(argv[1]);
        par.exp_name = argv[2];
    } else {
        cerr << "Error. Use: executable 0|1 filename" << endl;
        exit(1);
    }

    if (par.experiment != 0 && par.experiment != 1) {
        cerr << "\n[ERROR main] experiment to run not present\n";
        exit(1);
    }

    run_experiment(par);

    return 0;
}
