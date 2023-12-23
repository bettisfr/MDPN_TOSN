#ifndef TOSN_ALGORITHMS_H
#define TOSN_ALGORITHMS_H

#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <cmath>
#include <unordered_set>

#include "tsp/tsp.h"

#include "deployment.h"
#include "../input.h"

using namespace std;

class algorithms {

private:
    deployment *dep;

    vector<function<void(algorithms&)>> algorithm_functions = {
            &algorithms::approxTSPN_S,
            &algorithms::approxMPN_S,
            &algorithms::approxTSPN_M,
            &algorithms::approxMPN_M,

            &algorithms::approxTSPN_S_DOI,
            &algorithms::approxMPN_S_DOI,
            &algorithms::approxTSPN_M_DOI,
            &algorithms::approxMPN_M_DOI,

            &algorithms::approxTSPN_S_DTR,
            &algorithms::approxMPN_S_DTR,
            &algorithms::approxTSPN_M_DTR,
            &algorithms::approxMPN_M_DTR
    };

    solution internal_approxTSPN_S(double);

    solution internal_approxMPN_S(double);

    solution internal_approxTSPN_M(double);

    solution internal_approxMPN_M(double);

public:
    explicit algorithms(deployment *);

    void run_experiment(int, int);

    void approxTSPN_S();

    void approxMPN_S();

    void approxTSPN_M();

    void approxMPN_M();

    void approxMPN_S_DOI();

    void approxTSPN_S_DOI();

    void approxTSPN_M_DOI();

    void approxMPN_M_DOI();

    void approxTSPN_S_DTR();

    void approxMPN_S_DTR();

    void approxTSPN_M_DTR();

    void approxMPN_M_DTR();

    solution improved_tsp_neighbors(const vector<sensor> &, double);

    solution tsp_split(vector<tuple<point, int>>, const vector<double> &, point, const vector<sensor> &, bool);

    solution approxMPN(point, double);

    solution appro_alg_nei(vector<sensor> V, int jth, point depot, double radius);

    void DFS(int, unordered_set<int> &, unordered_set<int> &, vector<vector<int>>);

    static double get_distance(point, point);

    static double get_distance(const sensor &, point);

    static double get_distance(const sensor &, const sensor &);

    static int get_angle(const sensor &, point);

    bool is_within_radius(const sensor &, point);

    static bool is_within_radius_doi(const sensor &, point);

    vector<point> get_intersection_points(point, point, double);

    double tour_cost(vector<tuple<point, int>>, vector<double>, int, int, point, const vector<sensor> &);

    double compute_energy_hovering(sensor);

    double compute_hovering_time(sensor);

    void draw_result(vector<vector<tuple<point, int>>>, bool);
};

#endif //TOSN_ALGORITHMS_H
