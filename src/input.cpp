#include "input.h"

void print_parameters(const input &par) {
    cout << "Save=" << par.save << endl;
    cout << "Seed=" << par.seed << endl;
    cout << "Experiment name=" << par.exp_name << endl;

    map<int, string> experiment_str = {
            {0, "Default parameters"},
            {1, "Loaded parameters from config file"},
            {2, "Read parameters from command line"},
    };

    map<int, string> scenario_str = {
            {0, "Regular"},
            {1, "With DOI"},
            {2, "With variable DTR"},
    };

    map<int, string> wireless_technology_str = {
            {0, "WiFi-5"},
            {1, "WiFi-4"},
            {2, "Zigbee"},
            {3, "Bluetooth"},
    };

    map<int, string> algorithm_str = {
            {0, "TSPN_S"},
            {1, "MPN_S"},
            {2, "TSPN_M"},
            {3, "MPN_M"},
    };

    cout << "Experiment=" << par.experiment << " (" << experiment_str[par.experiment] << ")" << endl << endl;

    cout << "Area length=" << par.area_length << " m" << endl;
    cout << "Area width=" << par.area_width << " m" << endl;
    cout << "Number of sensors=" << par.num_sensors << endl;
    cout << "Number of depots=" << par.num_depots << endl;
    cout << "Radius of sensor=" << par.sensor_radius << endl;
    cout << "Radius DOI of sensor=" << par.sensor_radius_doi_percentage << " %" << endl;
    cout << "DOI=" << par.doi << endl;
    cout << "Max data=" << par.max_data << " MB" << endl;
    cout << "Energy budget=" << par.energy_budget << " J" << endl;
    cout << "Energy consumption for flying=" << par.energy_cons_fly << " J/m" << endl;
    cout << "Energy consumption for hovering=" << par.energy_cons_hover << " J/s" << endl;
    cout << "Wireless technology=" << wireless_technology_str[par.wireless_technology] << endl;
    cout << "Epsilon (budget violation)=" << par.epsilon << endl;

    cout << "Scenario=" << scenario_str[par.scenario] << endl;
    cout << "Algorithm=" << algorithm_str[par.algorithm] << endl;
    cout << "Iterations=" << par.iterations << endl;

    cout << endl << endl;
}

void save_parameters(const input &par) {
    string cfg_filename = "input/" + par.exp_name + ".cfg";
    ofstream file_cfg(cfg_filename);

    file_cfg << "save=" << par.save << endl;
    file_cfg << "seed=" << par.seed << endl;
    file_cfg << "area_length=" << par.area_length << endl;
    file_cfg << "area_width=" << par.area_width << endl;
    file_cfg << "num_sensors=" << par.num_sensors << endl;
    file_cfg << "num_depots=" << par.num_depots << endl;
    file_cfg << "sensor_radius=" << par.sensor_radius << endl;
    file_cfg << "sensor_radius_doi_percentage=" << par.sensor_radius_doi_percentage << endl;
    file_cfg << "doi=" << par.doi << endl;
    file_cfg << "max_data=" << par.max_data << endl;
    file_cfg << "energy_budget=" << par.energy_budget << endl;
    file_cfg << "energy_cons_fly=" << par.energy_cons_fly << endl;
    file_cfg << "energy_cons_hover=" << par.energy_cons_hover << endl;
    file_cfg << "wireless_technology=" << par.wireless_technology << endl;
    file_cfg << "epsilon=" << par.epsilon << endl;
    file_cfg << "scenario=" << par.scenario << endl;
    file_cfg << "algorithm=" << par.algorithm << endl;
    file_cfg << "iterations=" << par.iterations << endl;

    file_cfg << endl;

    file_cfg.close();
}

input load_parameters(input &par) {
    string cfg_filename = "../input/" + par.exp_name + ".cfg";
    ifstream file_cfg(cfg_filename);

    // Check if the file exists and can be opened
    if (!file_cfg.is_open()) {
        cerr << "Error opening config file: " << cfg_filename << endl;
        exit(-1);
    }

    string line;
    while (getline(file_cfg, line)) {
        stringstream lineStream(line);
        string key;
        string value;

        if (getline(lineStream, key, '=') && lineStream >> value) {
            if (key == "save") {
                par.save = stoi(value);
            } else if (key == "seed") {
                par.seed = stoi(value);
            } else if (key == "area_length") {
                par.area_length = stoi(value);
            } else if (key == "area_width") {
                par.area_width = stoi(value);
            } else if (key == "num_sensors") {
                par.num_sensors = stoi(value);
            } else if (key == "num_depots") {
                par.num_depots = stoi(value);
            } else if (key == "sensor_radius") {
                par.sensor_radius = stod(value);
            } else if (key == "sensor_radius_doi_percentage") {
                par.sensor_radius_doi_percentage = stod(value);
            } else if (key == "doi") {
                par.doi = stod(value);
            } else if (key == "max_data") {
                par.max_data = stoi(value);
            } else if (key == "energy_budget") {
                par.energy_budget = stoi(value);
            } else if (key == "energy_cons_fly") {
                par.energy_cons_fly = stod(value);
            } else if (key == "energy_cons_hover") {
                par.energy_cons_hover = stod(value);
            } else if (key == "wireless_technology") {
                par.wireless_technology = stod(value);
            } else if (key == "epsilon") {
                par.epsilon = stod(value);
            } else if (key == "scenario") {
                par.scenario = stoi(value);
            } else if (key == "algorithm") {
                par.algorithm = stoi(value);
            } else if (key == "iterations") {
                par.iterations = stoi(value);
            }
        }
    }

    file_cfg.close();

    return par;
}

input read_parameters(input &par, int argc, char* argv[]) {
    // Iterate through command line arguments
    // Start from 2 to skip program name and --params
    for (int i = 2; i < argc; ++i) {
        string arg = argv[i];

        // Check if the argument is a flag (e.g., starts with '-')
        if (!arg.empty() && arg[0] == '-') {
            // Process the flag or option accordingly
            if (arg == "-exp_name") {
                par.exp_name = argv[i + 1];
            } else if (arg == "-save") {
                par.save = stoi(argv[i + 1]);
            } else if (arg == "-seed") {
                par.seed = stoi(argv[i + 1]);
            } else if (arg == "-area_length") {
                par.area_length = stoi(argv[i + 1]);
            } else if (arg == "-area_width") {
                par.area_width = stoi(argv[i + 1]);
            } else if (arg == "-num_sensors") {
                par.num_sensors = stoi(argv[i + 1]);
            } else if (arg == "-num_depots") {
                par.num_depots = stoi(argv[i + 1]);
            } else if (arg == "-sensor_radius") {
                par.sensor_radius = stoi(argv[i + 1]);
            } else if (arg == "-sensor_radius_doi_percentage") {
                par.sensor_radius_doi_percentage = stod(argv[i + 1]);
            } else if (arg == "-doi") {
                par.doi = stod(argv[i + 1]);
            } else if (arg == "-max_data") {
                par.max_data = stoi(argv[i + 1]);
            } else if (arg == "-energy_budget") {
                par.energy_budget = stoi(argv[i + 1]);
            } else if (arg == "-energy_cons_fly") {
                par.energy_cons_fly = stoi(argv[i + 1]);
            } else if (arg == "-energy_cons_hover") {
                par.energy_cons_hover = stoi(argv[i + 1]);
            } else if (arg == "-wireless_technology") {
                par.wireless_technology = stoi(argv[i + 1]);
            } else if (arg == "-epsilon") {
                par.epsilon = stod(argv[i + 1]);
            } else if (arg == "-scenario") {
                par.scenario = stoi(argv[i + 1]);
            } else if (arg == "-algorithm") {
                par.algorithm = stoi(argv[i + 1]);
            } else if (arg == "-iterations") {
                par.iterations = stoi(argv[i + 1]);
            } else {
                cerr << "Unknown option: " << arg << endl;
            }
        }
    }

    return par;
}
