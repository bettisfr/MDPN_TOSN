#ifndef TOSN_INPUT_H
#define TOSN_INPUT_H

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

// Default values
struct input {
    // 0: No save to CSV; 1: Save to CSV
    int save = 1;

    // Application parameters
    int seed = 0;

    // 0: Default values; 1: From cfg file; 2 From command line
    int experiment = 0;

    // 0: Regular; 1: With DOI; 2: With variable DTR
    int scenario = 0;

    // 0: TSPN_S; 1: MPN_S; 2: TSPN_M; 3: MPN_M
    int algorithm = 0;

    // Number or random instances to be performed
    int iterations = 30;

    // Name of the experiment (just a name)
    string exp_name = "default";

    // Number of sensors
    int num_sensors = 10;

    // Number of depots
    int num_depots = 1;

    // Length of the area (in meters)
    int area_length = 2000;

    // Width of the area (in meters)
    int area_width = 2000;

    // Radius of the sensors (in case of no DOI)
    double sensor_radius = 20;

    // Radius of the sensors in case of DOI (reduced radius in percentage, e.g., 0.8 means 80% of sensor_radius)
    double sensor_radius_doi_percentage = 1;

    // Degree of Irregularity (DOI)
    // -> https://www.sciencedirect.com/science/article/pii/S1574119218305406 (Section 4)
    // DOI=0 sphere, DOI>0 irregularities
    double doi = 0;

    // Max data (in MB)
    int max_data = 1024;

    // The following 4 parameters
    // -> https://www.sciencedirect.com/science/article/pii/S0022000023000806 (Section 6)
    // Energy budget of each drone (in J)
    int energy_budget = 5000000;

    // Average energy consumption for flying for each meter traveled (in J/m)
    double energy_cons_fly = 150;

    // Average energy consumption for hovering for each second spent (in J/s = W)
    double energy_cons_hover = 300;

    // 0: WiFi-5; 1: WiFi-4; 2: Zigbee; 3: Bluetooth
    int wireless_technology = 1;

    // Budget violation constant (1 + \epsilon)
    // Impacts the performance of the algorithm: > 0 slow (no violation), 1 <= fast (2x budget required)
    double epsilon = 0.1;

};

void print_parameters(const input &);
void save_parameters(const input &);
input load_parameters(input &);
input read_parameters(input &, int, char* []);

#endif //TOSN_INPUT_H
