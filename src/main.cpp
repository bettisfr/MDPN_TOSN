#include <iostream>

#include "input.h"
#include "core/algorithms.h"
#include "core/deployment.h"
#include "model/energy.h"

void test();

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

void test() {
    energy en;

    // Drone's parameters
    double payload = 0.5;

    // Wind's parameters
    double wind_speed = 0;
    double wind_direction = 0;

    // en required for HOVERING (time in seconds)
    double time = 1;
    double e_hovering = en.get_energy_hovering(time, payload, wind_speed);
    cout << "en_Hovering=" << e_hovering << endl;

    // en required for MOVING (distance in meters)
    double speed = 1.75;
    // Points A and B (assuming to fly from A to B)
    // Distance between vertex A and vertex B
    double Ax = 0, Ay = 0, Bx = 1, By = 0;
    double distance = sqrt(pow(Ax - Bx, 2) + pow(Ay - By, 2));
    // Directions between vertex A and vertex B
    double direction = atan2(By - Ay, Bx - Ax) * 180 / M_PI;
    if (direction < 0) {
        direction += 360;
    }
    // Relative wind direction (difference between the two directions)
    double relative_wind_direction = fabs(direction - wind_direction);
    double e_moving = en.get_energy_movement(distance, payload, speed, wind_speed, relative_wind_direction);
    cout << "en_Moving=" << e_moving << endl;

    exit(0);
}

int main(int argc, char** argv) {

//    test();

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
