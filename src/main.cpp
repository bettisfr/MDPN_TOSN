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
    // cout << dep << endl;
    // Algorithms creation and invocation
    algorithms alg(&dep);
    
    // vector<sensor> sensors = dep.get_sensors();
    // sort(sensors.begin(), sensors.end(), [](const sensor& a, const sensor& b) {
    //     return a.get_pos_x() < b.get_pos_x();
    // });
  
    // auto [tspn_tour_im, tspn_cost_im] = alg.improved_tsp_neighbors(sensors, dep.get_sensor_radius());
    // auto [tspn_tour, tspn_cost] = alg.tsp_neighbors(sensors, dep.get_sensor_radius());

    // double total_cost_im;
    // for (auto c : tspn_cost_im){
    //     total_cost_im += c;
    // }
    // cout << "improve_tspn_cost " << total_cost_im << endl;

    // double total_cost;
    // for (auto c : tspn_cost){
    //     total_cost += c;
    // }
    // cout << "tspn_cost " << total_cost << endl;

    alg.approxTSPN_S();
    //alg.approxMPN_S();
    //alg.approxTSPN_M();
    //alg.approxMPN_M();
    
    // with DOI
    //alg.approxTSPN_S_doi();

}

int main(int argc, char** argv) {

    input par;

    if (argc >= 3) {
        try {
            par.experiment = stoi(argv[1]);
        } catch (const invalid_argument& e) {
            cerr << "Error. Experiment number must be an integer." << endl;
            exit(1);
        } catch (const out_of_range& e) {
            cerr << "Error. Experiment number is out of range." << endl;
            exit(1);
        }

        par.exp_name = argv[2];
    } else {
        cerr << "Error. Use: " << argv[0] << " 0|1 filename" << endl;
        exit(1);
    }

    if (par.experiment != 0 && par.experiment != 1) {
        cerr << "\n[ERROR main] Experiment to run not present\n";
        exit(1);
    }

    run_experiment(par);

    return 0;
}
