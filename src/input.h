#ifndef TOSN_INPUT_H
#define TOSN_INPUT_H

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>
#include <fstream>
#include <sstream>

using namespace std;

// Default values
struct input {
    // Application parameters
    int seed = 0;

    // 0: Default values; 1: From cfg file
    int experiment = 0;

    // 0: Regular; 1: With DOI; 2: With variable DTR
    int scenario = 0;

    // 0: TSPN_S; 1: MPN_S; 2: TSPN_M; 3: MPN_M
    int algorithm = 0;

    // Name of the experiment (just a name)
    string exp_name = "default";

    // Number of sensors
    int num_sensors = 10;

    // Number of depots
    int num_depots = 1;

    // Length of the area (in meters)
    int area_length = 100;

    // Width of the area (in meters)
    int area_width = 100;

    // Radius of the sensors (in case of no DOI)
    double sensor_radius = 20;

    // Radius of the sensors in case of DOI (reduced radius)
    double sensor_radius_doi = 10;

    // Degree of Irregularity (DOI)
    // -> https://www.sciencedirect.com/science/article/pii/S1574119218305406 (Section 4)
    // DOI=0 sphere, DOI>0 irregularities
    double doi = 0;

    // Max data (in MB)
    int max_data = 1024;

    // The following 4 parameters
    // -> https://www.sciencedirect.com/science/article/pii/S0022000023000806 (Section 6)
    // Energy budget of each drone (in kJ)
    int energy_budget = 5000;

    // Average energy consumption for flying for each meter traveled (in J/m)
    double energy_cons_fly = 200;

    // Average energy consumption for hovering for each second spent (in J/s = W)
    double energy_cons_hover = 700;

    // Maximum data transfer rate between drone and sensors (in MB/s)
    double data_transfer_rate = 50.0;

    // Budget violation constant (1 + \epsilon)
    // Impacts the performance of the algorithm: > 0 slow (no violation), 1 <= fast (2x budget required)
    double epsilon = 0.5;

};

void print_parameters(const input &);
void save_parameters(const input &);
input load_parameters(input &);

#endif //TOSN_INPUT_H
