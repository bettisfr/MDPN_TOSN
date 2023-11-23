#ifndef TOSN_ALGORITHMS_H
#define TOSN_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>

#include "tsp/tsp.h"

#include "deployment.h"
#include "../input.h"

using namespace std;

class algorithms {

private:
    deployment *dep;
    // change this in a way that you can also store the id of sensor
    vector<tuple<double, double, int>> tspn_result;

public:
    explicit algorithms(deployment *);

    double get_distance(sensor, point);

    int get_angle(sensor, point);

    bool is_within_radius(const sensor &, point);

    bool is_within_radius_doi(const sensor &, point);

    vector<point> get_intersection_points(point, point);

    void tsp_neighbors();

    vector<tuple<double, double, double, double>> sorted_sensors;

    // TODO: random names now
    void algorithm_1();

    void algorithm_2();

    void draw_result();
};


#endif //TOSN_ALGORITHMS_H
