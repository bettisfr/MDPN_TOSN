#ifndef TOSN_ALGORITHMS_H
#define TOSN_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <cmath>

#include "deployment.h"
#include "../input.h"

class algorithms {
private:
    int num_sensors;
    double sensor_radius;
    vector<sensor> sensors;
    vector<depot> depots;

    double get_distance(point, point);

public:
    explicit algorithms(deployment*);

    double get_distance(sensor, sensor);
    double get_distance(depot, sensor);

    // TODO: random names now
    void algorithm_1();
    void algorithm_2();
};


#endif //TOSN_ALGORITHMS_H
