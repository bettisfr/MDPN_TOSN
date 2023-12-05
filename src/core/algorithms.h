#ifndef TOSN_ALGORITHMS_H
#define TOSN_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <unordered_set>

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

    double tour_cost(vector<tuple<double, double, int>>, int, int, point);

    double compute_energy_hovering(tuple<double, double, double, double>);

    void tsp_neighbors(vector<sensor>);
    void tsp_split(int, point);
    void ApproxMPN(point);
    vector<vector<tuple<double, double, int>>> approAlgNei(vector<tuple<double, double, double, double>>, int);
    void DFS(int, unordered_set<int>&, unordered_set<int>&, vector<vector<int>>);

    void ApproxTSPN_S();
    void ApproxMPN_S();

    void ApproxTSPN_M();
    void ApproxMPN_M();

    void draw_result();
};


#endif //TOSN_ALGORITHMS_H
