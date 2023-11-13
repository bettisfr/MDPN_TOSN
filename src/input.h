#ifndef TOSN_INPUT_H
#define TOSN_INPUT_H

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>

using namespace std;

// Default values
struct input {
    // Application parameters
    int seed = 0;

    // 0: to do
    int experiment = 0;

    // Name of the experiment (just a name)
    string exp_name = "default";

    // Number of sensors
    int num_sensors;

    // Length of the area (in meters)
    int area_length;

    // Width of the area (in meters)
    int area_width;

    // Radius of the sensors (in case of no DOI)
    double sensor_radius;

    // Degree of Irregularity (DOI)
    double doi;

};

#endif //TOSN_INPUT_H
