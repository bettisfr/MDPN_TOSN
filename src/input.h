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

    // 0: to do
    int experiment = 0;

    // Name of the experiment (just a name)
    string exp_name = "default";

    // Number of sensors
    int num_sensors = 100;

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

    // Number of random scenarios
    int scenarios = 10;

};

map<int, string> exp_str = {
    {0, "Default parameters"},
    {1, "From loaded parameters"},
};

void print_parameters(const input &);
void save_parameters(const input &);
input load_parameters(input &);

void print_parameters(const input &par) {
    cout << "Seed=" << par.seed << endl;
    cout << "Experiment name=" << par.exp_name << endl;
    cout << "Number of scenarios=" << par.scenarios << endl;

    cout << "Experiment=" << par.experiment << " (" << exp_str[par.experiment] << ")" << endl << endl;

    cout << "Area length=" << par.area_length << "m" << endl;
    cout << "Area width=" << par.area_width << "m" << endl;
    cout << "Number of sensors=" << par.num_sensors << endl;
    cout << "Radius of sensor=" << par.sensor_radius << endl;
    cout << "DOI=" << par.doi << endl;

    cout << endl << endl;
}

void save_parameters(const input &par) {
    string cfg_filename = "input/" + par.exp_name + ".cfg";
    ofstream file_cfg(cfg_filename);

    file_cfg << "seed=" << par.seed << endl;
    file_cfg << "area_length=" << par.area_length << endl;
    file_cfg << "area_width=" << par.area_width << endl;
    file_cfg << "num_sensors=" << par.num_sensors << endl;
    file_cfg << "sensor_radius=" << par.sensor_radius << endl;
    file_cfg << "doi=" << par.doi << endl;
    file_cfg << "scenarios=" << par.scenarios << endl;

    file_cfg << endl;

    file_cfg.close();
}

input load_parameters(input &par) {
    string cfg_filename = "input/" + par.exp_name + ".cfg";
    ifstream file_cfg(cfg_filename);

    string line;
    while (getline(file_cfg, line)) {
        stringstream lineStream(line);
        string key;
        string value;

        if (getline(lineStream, key, '=') && lineStream >> value) {
            if (key == "seed") {
                par.seed = stoi(value);
            } else if (key == "area_length") {
                par.area_length = stoi(value);
            } else if (key == "area_width") {
                par.area_width = stoi(value);
            } else if (key == "num_sensors") {
                par.num_sensors = stoi(value);
            } else if (key == "sensor_radius") {
                par.sensor_radius = stod(value);
            } else if (key == "doi") {
                par.doi = stod(value);
            } else if (key == "scenarios") {
                par.scenarios = stoi(value);
            }
        }
    }

    file_cfg.close();

    return par;
}




#endif //TOSN_INPUT_H
