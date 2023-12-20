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

    solution approxTSPN_S(double);
    solution approxMPN_S(double);
    solution approxTSPN_M(double);
    solution approxMPN_M(double);

    void approxMPN_S_doi();
    void approxTSPN_S_doi();
    void approxTSPN_M_doi();
    void approxMPN_M_doi();

    void approxTSPN_S_xxx();

    //tuple<vector<tuple<point, int>>, vector<double>> tsp_neighbors(const vector<sensor>&, double);
    solution improved_tsp_neighbors(const vector<sensor>&, double);

    solution tsp_split(vector<tuple<point, int>>, const vector<double>&, point, const vector<sensor>&, bool);

    solution approxMPN(point, double);

    solution appro_alg_nei(vector<sensor> V, int jth, point depot, double radius);

    void DFS(int, unordered_set<int> &, unordered_set<int> &, vector<vector<int>>);

    static double get_distance(point, point);

    static double get_distance(const sensor&, point);

    static double get_distance(const sensor&, const sensor&);

    static int get_angle(const sensor&, point);

    bool is_within_radius(const sensor &, point);

    static bool is_within_radius_doi(const sensor &, point);

    vector<point> get_intersection_points(point, point, double);

    double tour_cost(vector<tuple<point, int>>, vector<double>, int, int, point, const vector<sensor>&);

    double compute_energy_hovering(sensor);
    double compute_hovering_time(sensor);

    void draw_result(vector<vector<tuple<point, int>>>, bool);
};

#endif //TOSN_ALGORITHMS_H
