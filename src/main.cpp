#include <iostream>

#include "input.h"
#include "core/sensor.h"
#include "core/deployment.h"

using namespace std;

void run_experiment(input &par) {
//    int n = 100;
//    int length = 1000;
//    int width = 1000;
//    double radius = 150;
//    double doi = 1;
//
//    deployment dep(n, length, width, radius, doi);

    if (par.experiment != 0) {
        par = load_parameters(par);
    }

    // Print used parameters
    print_parameters(par);
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
