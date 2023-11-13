#ifndef TOSN_DEPLOYMENT_H
#define TOSN_DEPLOYMENT_H

#include <iostream>
#include <vector>

#include "../input.h"
#include "sensor.h"

using namespace std;

typedef tuple<double, double> point;

class deployment {
private:
    int num_sensors;
    int num_depots;
    int area_length;
    int area_width;
    int max_data;
    double sensor_radius;
    double doi;

    vector<sensor> sensors;
    vector<point> depots;

public:
    deployment(const input &);

    int get_num_sensors() const;
    int get_area_length() const;
    int get_area_width() const;
    double get_sensor_radius() const;
    vector<sensor> get_sensors();
    vector<point> get_depots();
};


#endif //TOSN_DEPLOYMENT_H
