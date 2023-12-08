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

public:
    explicit algorithms(deployment *);

    static double get_distance(point, point);

    static double get_distance(const sensor&, point);

    static double get_distance(const sensor&, const sensor&);

    static int get_angle(const sensor&, point);

    bool is_within_radius(const sensor &, point);

    static bool is_within_radius_doi(const sensor &, point);

    vector<point> get_intersection_points(point, point);

    double tour_cost(vector<tuple<point, int>>, vector<double>, int, int, point);

    double compute_energy_hovering(tuple<double, double, double, double>);

    double compute_energy_hovering(sensor);

    tuple<vector<tuple<point, int>>, vector<double>> tsp_neighbors(const vector<sensor>&);

    tuple<vector<vector<tuple<point, int>>>, vector<double>> tsp_split(vector<tuple<point, int>>, const vector<double>&, point);

    void approxMPN(point depot);

    vector<vector<tuple<double, double, int>>> approAlgNei(vector<tuple<double, double, double, double>>, int);

    void DFS(int, unordered_set<int> &, unordered_set<int> &, vector<vector<int>>);

    void approxTSPN_S();

    void approxMPN_S();

    void approxTSPN_M();

    void approxMPN_M();

    void draw_result(vector<vector<tuple<point, int>>>);
};

#endif //TOSN_ALGORITHMS_H
