#ifndef TOSN_DEPLOYMENT_H
#define TOSN_DEPLOYMENT_H

#include <iostream>
#include <vector>
#include <random>
#include <numeric>

#include "../input.h"
#include "sensor.h"

using namespace std;

typedef tuple<double, double> point;

class deployment {

    friend ostream &operator<<(ostream &os, const deployment &s);

private:
    int num_sensors;
    int num_depots;
    int area_length;
    int area_width;
    int max_data;
    double sensor_radius;
    double doi;

    const double c = 0.3e+9;
    const double f_c = 6.5e+9;
    const double shape = 0.67;
    const double scale = 0.16;

    vector<sensor> sensors;
    vector<point> depots;

public:
    explicit deployment(const input &);

    int get_num_sensors() const;
    int get_area_length() const;
    int get_area_width() const;
    double get_sensor_radius() const;
    vector<sensor> get_sensors();
    vector<point> get_depots();
};


#endif //TOSN_DEPLOYMENT_H
