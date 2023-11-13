#include "deployment.h"

deployment::deployment(int n_sensors, int length, int width, double radius, double _doi) {
    num_sensors = n_sensors;
    area_length = length;
    area_width = width;
    sensor_radius = radius;
    doi = _doi;

    for (int i = 0; i < num_sensors; i++) {
        double x = 0;
        double y = 0;
        double data = 0;
        vector<double> r_doi(360);
        sensors.emplace_back(x, y, data, r_doi);
    }
}
