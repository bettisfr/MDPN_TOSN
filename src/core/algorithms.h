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

    vector<function<solution(algorithms &)>> algorithm_functions = {
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

    static vector<point> get_line_circle_intersections_helper(const point &, const point &, double, const point &);

public:
    explicit algorithms(deployment *);

    solution run_experiment(int, int);

    solution approxTSPN_S();

    solution approxMPN_S();

    solution approxTSPN_M();

    solution approxMPN_M();

    solution approxMPN_S_DOI();

    solution approxTSPN_S_DOI();

    solution approxTSPN_M_DOI();

    solution approxMPN_M_DOI();

    solution approxTSPN_S_DTR();

    solution approxMPN_S_DTR();

    solution approxTSPN_M_DTR();

    solution approxMPN_M_DTR();

    solution tsp_neighbors_v1(const vector<sensor> &, double);

    solution tsp_neighbors_v2(const vector<sensor> &, double);

    solution tsp_split(vector<tuple<point, int>>, const vector<double> &, point, const vector<sensor> &, bool);

    solution approxMPN(point, double);

    solution appro_alg_nei(vector<sensor>, int, point, double);

    void DFS(int, unordered_set<int> &, unordered_set<int> &, vector<vector<int>>);

    static double get_distance(point, point);

    static double get_distance(const sensor &, point);

    static double get_distance(const sensor &, const sensor &);

    static int get_angle(const sensor &, point);

    bool is_within_radius(const sensor &, point);

    static bool is_within_radius_doi(const sensor &, point);

    static vector<point> get_circles_intersections(const point &, const point &, double);

    static vector<point> get_line_circle_intersections(const point &, const point &, double);

    static vector<point> get_line_circle_intersections(const point &, const point &, double, const point &);

    static vector<point> get_line_segment_circle_intersections(const point &, const point &, double, const point &);

    double tour_cost(vector<tuple<point, int>>, vector<double>, int, int, point, const vector<sensor> &);

    double compute_energy_hovering(sensor);

    double compute_hovering_time(sensor);

    int compute_uncovered_sensors(const solution &);

    tuple<double, double> compute_lost_data(const solution &);

    void draw_result(vector<vector<tuple<point, int>>>, bool, bool);
};

#endif //TOSN_ALGORITHMS_H
