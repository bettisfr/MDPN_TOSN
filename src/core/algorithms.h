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
    vector<tuple<double, double, int>> tspn_result;
    vector<tuple<double, double, double, double>> sorted_sensors;
    vector<double> tspn_cost;
    vector<vector<tuple<double, double, int>>> tspn_tours;

public:
    explicit algorithms(deployment *);

    double get_distance(sensor, point);

    int get_angle(sensor, point);

    bool is_within_radius(const sensor &, point);

    bool is_within_radius_doi(const sensor &, point);

    vector<point> get_intersection_points(point, point);

    double tour_cost(vector<tuple<double, double, int>>);

    void tsp_neighbors();
    void tsp_split();

    // TODO: random names now
    void ApproxTSPN_S();

    void algorithm_2();

    void draw_result();
};


#endif //TOSN_ALGORITHMS_H
