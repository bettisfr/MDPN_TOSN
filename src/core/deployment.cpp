#include <random>
#include "deployment.h"

deployment::deployment(const input &par) {
    num_sensors = par.num_sensors;
    num_depots = par.num_depots;
    area_length = par.area_length;
    area_width = par.area_width;
    max_data = par.max_data;
    sensor_radius = par.sensor_radius;
    doi = par.doi;

    mt19937 re(par.seed);
    uniform_real_distribution<double> length_rand(0, area_length);
    uniform_real_distribution<double> width_rand(0, area_width);
    uniform_int_distribution<int> data_rand(0, max_data);

    // Sensors creation
    for (int i = 0; i < num_sensors; i++) {
        double x = length_rand(re);
        double y = width_rand(re);
        double data = data_rand(re);

        // Suitably fix this according to Degree of Irregularity (DOI)
        // https://www.sciencedirect.com/science/article/pii/S1574119218305406 (Section 4)
        // DOI=0 sphere, DOI>0 irregularities
        vector<double> r_doi(360);
        for (int angle = 0; angle < 360; angle++) {
            // FIXME
            r_doi[angle] = sensor_radius;
        }

        sensors.emplace_back(x, y, data, r_doi);
    }

    // Depots creation
    for (int i = 0; i < num_depots; i++) {
        double x = length_rand(re);
        double y = width_rand(re);

        depots.emplace_back(x, y);
    }
}

int deployment::get_num_sensors() const {
    return num_sensors;
}

int deployment::get_area_length() const {
    return area_length;
}

int deployment::get_area_width() const {
    return area_width;
}

double deployment::get_sensor_radius() const {
    return sensor_radius;
}

vector<sensor> deployment::get_sensors() {
    return sensors;
}

vector<point> deployment::get_depots() {
    return depots;
}

