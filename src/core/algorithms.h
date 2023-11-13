#ifndef TOSN_ALGORITHMS_H
#define TOSN_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <cmath>

#include "deployment.h"
#include "../input.h"

using namespace std;

class algorithms {

private:
    int num_sensors;
    double sensor_radius;
    vector<sensor> sensors;
    vector<point> depots;

public:
    explicit algorithms(deployment*);

    double get_distance(sensor, point);
    int get_angle(sensor, point);
    bool is_within_radius(const sensor&, point);
    bool is_within_radius_doi(const sensor&, point);

    // TODO: random names now
    void algorithm_1();
    void algorithm_2();
};


#endif //TOSN_ALGORITHMS_H
