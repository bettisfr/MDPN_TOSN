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
    int area_length;
    int area_width;
    double sensor_radius;
    vector<sensor> sensors;
    vector<point> depots;

    vector<int> tsp_result_id;
    vector<point_3d> tsp_result;

public:
    explicit algorithms(deployment*);

    double get_distance(sensor, point);
    int get_angle(sensor, point);
    bool is_within_radius(const sensor&, point);
    bool is_within_radius_doi(const sensor&, point);

    vector<point> get_intersection_points(point, point);
    void tsp_neighbors();

    // TODO: random names now
    void algorithm_1();
    void algorithm_2();

    void draw_result();
};


#endif //TOSN_ALGORITHMS_H
