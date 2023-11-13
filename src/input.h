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

    // Number of random scenarios
    int scenarios = 10;

    // 0: to do
    int experiment = 0;

    // Name of the experiment (just a name)
    string exp_name = "default";

    // Number of sensors
    int num_sensors = 100;

    // Number of depots
    int num_depots = 1;

    // Length of the area (in meters)
    int area_length = 1000;

    // Width of the area (in meters)
    int area_width = 1000;

    // Radius of the sensors (in case of no DOI)
    double sensor_radius = 150;

    // Degree of Irregularity (DOI)
    // https://www.sciencedirect.com/science/article/pii/S1574119218305406 (Section 4)
    // DOI=0 sphere, DOI>0 irregularities
    double doi = 0;

    // Max data (in MB)
    int max_data = 1024;

};

void print_parameters(const input &);
void save_parameters(const input &);
input load_parameters(input &);

#endif //TOSN_INPUT_H
