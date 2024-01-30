#include "algorithms.h"
#include <vector>
#include <unordered_set>

algorithms::algorithms(deployment *m_dep) {
    dep = m_dep;
}

// Useless (and stupid) method, but it was nice to use a vector of pointers to methods :-D
solution algorithms::run_experiment(int scenario, int algorithm) {
    int index = scenario * 4 + algorithm;
    solution out;

    if (index >= 0 && index < 12) {
        out = algorithm_functions[index](*this);
    } else {
        cerr << "Invalid scenario or algorithm index." << endl;
    }

    return out;
}

solution algorithms::approxTSPN_S() {
    solution sol = internal_approxTSPN_S(dep->get_sensor_radius());

    draw_result(sol.tours, true, false);

    sol.uncovered_sensors = 0; // no DOI
    sol.lost_data = 0; // no DTR

    return sol;
}

solution algorithms::internal_approxTSPN_S(double radius) {
    //solution sol = tsp_neighbors_v1(dep->get_sensors(), radius);    //////////////
    solution sol = tsp_neighbors_v2(dep->get_sensors(), radius);


    sol = tsp_split(sol.tours[0], sol.tours_costs, dep->get_depots()[0], dep->get_sensors(), false);

    sol.tours_number = static_cast<int>(sol.tours.size());
    sol.tours_costs = sol.tours_costs;
    sol.total_sensors = dep->get_num_sensors();
    sol.uncovered_sensors = compute_uncovered_sensors(sol);
    auto [lost_data, total_data] = compute_lost_data(sol);
    sol.total_data = total_data;
    sol.lost_data = lost_data;

    return sol;
}

solution algorithms::approxMPN_S() {
    solution sol = internal_approxMPN_S(dep->get_sensor_radius());

//    draw_result(sol.tours, true, false);

    sol.uncovered_sensors = 0; // no DOI
    sol.lost_data = 0; // no DTR

    return sol;
}

solution algorithms::internal_approxMPN_S(double radius) {
    solution sol = approxMPN(dep->get_depots()[0], dep->get_sensor_radius());

    sol.tours_number = static_cast<int>(sol.tours.size());
    sol.tours_costs = sol.tours_costs;
    sol.total_sensors = dep->get_num_sensors();
    sol.uncovered_sensors = compute_uncovered_sensors(sol);
    auto [lost_data, total_data] = compute_lost_data(sol);
    sol.total_data = total_data;
    sol.lost_data = lost_data;

    return sol;
}

solution algorithms::approxTSPN_M() {
    solution sol = internal_approxTSPN_M(dep->get_sensor_radius());

//    draw_result(sol.tours, false, false);

    sol.uncovered_sensors = 0; // no DOI
    sol.lost_data = 0; // no DTR

    return sol;
}

solution algorithms::internal_approxTSPN_M(double radius) {
    double ecf = dep->get_energy_cons_fly();
    //double radius = dep->get_sensor_radius();
    double R_0_f = radius * ecf;
    vector<point> depots = dep->get_depots();
    vector<double> A_d(depots.size(), 0.0);

    // find the depot
    for (int i = 0; i < depots.size(); i++) {
        for (const auto &sensor: dep->get_sensors()) {
            double dist = get_distance(sensor, depots[i]);
            //auto pos = sensor.get_position();
            double energy_hovering = compute_energy_hovering(sensor);
            if (dist <= radius) {
                if (energy_hovering > A_d[i]) {
                    A_d[i] = energy_hovering;
                }
            } else {
                double a = 2 * (dist * ecf + (energy_hovering / 2.) + R_0_f);
                if (a > A_d[i]) {
                    A_d[i] = a;
                }
            }
        }
    }

    auto index = distance(A_d.begin(), min_element(A_d.begin(), A_d.end()));
    point depot = depots[index];

    solution sol = tsp_neighbors_v1(dep->get_sensors(), radius);
    sol = tsp_split(sol.tours[0], sol.tours_costs, depot, dep->get_sensors(), false);

    sol.tours_number = static_cast<int>(sol.tours.size());
    sol.total_sensors = dep->get_num_sensors();
    sol.uncovered_sensors = compute_uncovered_sensors(sol);
    auto [lost_data, total_data] = compute_lost_data(sol);
    sol.total_data = total_data;
    sol.lost_data = lost_data;

    return sol;
}

solution algorithms::approxMPN_M() {
    solution sol = internal_approxMPN_M(dep->get_sensor_radius());

//    draw_result(sol.tours, false, false);

    sol.uncovered_sensors = 0; // no DOI
    sol.lost_data = 0; // no DTR

    return sol;
}

solution algorithms::internal_approxMPN_M(double radius) {
    solution sol;
    vector<point> depots = dep->get_depots();
    vector<tuple<point, int>> vec;
    vec.emplace_back(make_tuple(0.0, 0.0), 0);

    for (int i = 0; i < 2 * dep->get_sensors().size(); i++) {
        sol.tours.push_back(vec);
    }

    for (auto depot: depots) {
        solution tmp = approxMPN(depot, radius);
        if (tmp.tours.size() < sol.tours.size()) {
            sol.tours = tmp.tours;
            sol.tours_costs = tmp.tours_costs;
        }
    }

    sol.tours_number = static_cast<int>(sol.tours.size());
    sol.total_sensors = dep->get_num_sensors();
    sol.uncovered_sensors = compute_uncovered_sensors(sol);
    auto [lost_data, total_data] = compute_lost_data(sol);
    sol.total_data = total_data;
    sol.lost_data = lost_data;

    return sol;
}

int algorithms::compute_uncovered_sensors(const solution &sol) {
    vector<sensor> uncovered_sensors;
    vector<int> uncovered_ids;
    for (auto tour: sol.tours) {
        for (int i = 1; i < tour.size() - 1; i++) {
            auto point_id = get<0>(tour[i]);
            int id = get<1>(tour[i]);
            if (!is_within_radius_doi(dep->get_sensors()[id], point_id)) {
                uncovered_sensors.push_back(dep->get_sensors()[id]);
                uncovered_ids.push_back(get<1>(tour[i]));
            }
        }
    }

    return static_cast<int>(uncovered_sensors.size());
}

solution algorithms::approxTSPN_S_DOI() {
    solution sol = internal_approxTSPN_S(dep->get_sensor_radius_doi());

//    draw_result(sol.tours, true, true);

    sol.lost_data = 0; // no DTR

    return sol;
}


solution algorithms::approxMPN_S_DOI() {
    solution sol = internal_approxMPN_S(dep->get_sensor_radius_doi());

//    draw_result(sol.tours, true, true);

    sol.lost_data = 0; // no DTR

    return sol;
}

solution algorithms::approxTSPN_M_DOI() {
    solution sol = internal_approxTSPN_M(dep->get_sensor_radius_doi());

//    draw_result(sol.tours, false, true);

    sol.lost_data = 0; // no DTR

    return sol;
}

solution algorithms::approxMPN_M_DOI() {
    solution sol = internal_approxMPN_M(dep->get_sensor_radius_doi());

//    draw_result(sol.tours, false, true);

    sol.lost_data = 0; // no DTR

    return sol;
}

tuple<double, double> algorithms::compute_lost_data(const solution& sol) {
    double lost_data = 0.;
    double total_data = 0.;

    for (auto tour: sol.tours) {
        for (int i = 1; i < tour.size() - 1; i++) {
            auto point_id = get<0>(tour[i]);
            int id = get<1>(tour[i]);
            total_data += dep->get_sensors()[id].get_data_size();
            double hovering_time = compute_hovering_time(dep->get_sensors()[id]);
            // distance of the point from sensor
            double dist = get_distance(dep->get_sensors()[id], point_id);
            double collected_data = hovering_time * dep->get_DTR(dist);
//            cout << "dist: " << dist << " DTR:" << dep->get_DTR(dist) << endl;
            double lost_data_i = dep->get_sensors()[id].get_data_size() - collected_data;
            lost_data += lost_data_i;
        }
    }

    return make_tuple(lost_data, total_data);
}

solution algorithms::approxTSPN_S_DTR() {
    solution sol = internal_approxTSPN_S(dep->get_sensor_radius());

//    draw_result(sol.tours, true, false);

    sol.uncovered_sensors = 0; // no DOI

    return sol;
}

solution algorithms::approxMPN_S_DTR() {
    solution sol = internal_approxMPN_S(dep->get_sensor_radius());

//    draw_result(sol.tours, true, false);

    sol.uncovered_sensors = 0; // no DOI

    return sol;
}

solution algorithms::approxTSPN_M_DTR() {
    solution sol = internal_approxTSPN_M(dep->get_sensor_radius());

//    draw_result(sol.tours, false, false);

    sol.uncovered_sensors = 0; // no DOI

    return sol;
}

solution algorithms::approxMPN_M_DTR() {
    solution sol = internal_approxMPN_M(dep->get_sensor_radius());

//    draw_result(sol.tours, false, false);

    sol.uncovered_sensors = 0; // no DOI

    return sol;
}


double algorithms::get_distance(point p1, point p2) {
    double x1, y1, x2, y2;
    tie(x1, y1) = p1;
    tie(x2, y2) = p2;

    double delta_x = x2 - x1;
    double delta_y = y2 - y1;

    return sqrt(delta_x * delta_x + delta_y * delta_y);
}

double algorithms::get_distance(const sensor &s, point p) {
    point pos_s = s.get_position();

    double x1, y1, x2, y2;
    tie(x1, y1) = p;
    tie(x2, y2) = pos_s;

    double delta_x = x2 - x1;
    double delta_y = y2 - y1;

    return sqrt(delta_x * delta_x + delta_y * delta_y);
}

double algorithms::get_distance(const sensor &s1, const sensor &s2) {
    point pos_s1 = s1.get_position();
    point pos_s2 = s2.get_position();

    double x1, y1, x2, y2;
    tie(x1, y1) = pos_s1;
    tie(x2, y2) = pos_s2;

    double delta_x = x2 - x1;
    double delta_y = y2 - y1;

    return sqrt(delta_x * delta_x + delta_y * delta_y);
}

int algorithms::get_angle(const sensor &s, point p) {
    point pos_s = s.get_position();

    double x1, y1, x2, y2;
    tie(x1, y1) = p;
    tie(x2, y2) = pos_s;

    double delta_x = x2 - x1;
    double delta_y = y2 - y1;

    double angle_rad = atan2(delta_y, delta_x);
    double angle_deg = angle_rad * (180.0 / M_PI);

    int rounded_angle = static_cast<int>(round(angle_deg));

    if (rounded_angle < 0) {
        rounded_angle += 360;
    }

    return rounded_angle;
}

bool algorithms::is_within_radius(const sensor &s, point p) {
    double dist = get_distance(s, p);
    return (dist <= dep->get_sensor_radius());
}

bool algorithms::is_within_radius_doi(const sensor &s, point p) {
    double dist = get_distance(s, p);

    int angle = get_angle(s, p);
    double actual_radius = s.get_radius_doi(angle);

    return (dist <= actual_radius);
}

double algorithms::compute_energy_hovering(sensor s) {
    double dtr = dep->get_DTR(0);
    double ech = dep->get_energy_cons_hover();
    double sensor_data = s.get_data_size();
    double required_time = sensor_data / dtr;
    return required_time * ech;
}

double algorithms::compute_hovering_time(sensor s) {
    double dtr = dep->get_DTR(0);
    double ech = dep->get_energy_cons_hover();
    double sensor_data = s.get_data_size();
    return double(sensor_data / dtr);
}


double algorithms::tour_cost(vector<tuple<point, int>> T, vector<double> tspn_cost, int start, int end, point depot,
                             const vector<sensor> &sensors) {
    vector<sensor> deployed_sensors = sensors; //dep->get_sensors();

    vector<sensor> orig_sensors = dep->get_sensors();

    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    double ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_DTR(0); // constant DTR (in MB/s)

    double cost_T_k = 0;

    auto lastIndex = T.size() - 2;

    point p0 = get<0>(T[0]);
    point p1 = get<0>(T[1]);
    point pLast = get<0>(T[lastIndex]);

    double dist1 = get_distance(p0, p1);
    double energy1_flying = dist1 * ecf;

    sensor sensor1 = orig_sensors[get<1>(T[1])];
    double energy1_hovering = compute_energy_hovering(sensor1);

    double dist2 = get_distance(p0, pLast);
    double energy2_flying = dist2 * ecf;

    sensor sensor2 = orig_sensors[get<1>(T[lastIndex])];
    double energy2_hovering = compute_energy_hovering(sensor2);

    for (int p = start; p < end; p++) {
        cost_T_k = cost_T_k + tspn_cost[p];
    }
    cost_T_k = cost_T_k + energy1_flying + energy2_flying + energy1_hovering / 2 + energy2_hovering / 2;

    return cost_T_k;
}

solution algorithms::approxMPN(point depot, double radius) {
    vector<sensor> deployed_sensors = dep->get_sensors();

    // get eps and compute t
    double epsilon = dep->get_epsilon();
    int t = ceil(log2(1 / epsilon));

    double energy_budget = dep->get_energy_budget();

    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    double ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_DTR(0); // constant DTR (in MB/s)
    //double radius = dep->get_sensor_radius();
    double R_0_f = radius * ecf;  // R_0

    vector<vector<sensor>> V(t + 1);

    // for each s in deployed_sensors compute lambd_v,d and j
    for (const auto& s: deployed_sensors) {
        double dist = get_distance(s, depot);
        double energy_flying = dist * ecf;

        double energy_hovering = compute_energy_hovering(s);
        double lambda = energy_flying + energy_hovering / 2. + R_0_f;

        int j = floor(log2((energy_budget - 2. * lambda) / (epsilon * energy_budget)) + 1.);

        // find V_0,..., V_t
        if ((static_cast<double>((1. - epsilon) * (energy_budget / 2.)) < lambda) &&
            (lambda <= static_cast<double>(energy_budget / 2.))) {
            V[0].push_back(s);
        } else {
            V[j].push_back(s);
        }
    }

    solution sol;

    for (int i = 0; i < V.size(); i++) {
        // for each V[i], run Minimum UAV Deployment Problem with Neighborhoods with budget 2^{j-1} epsilon B
        solution tmp;

        if (V[i].size() == 1) {
            sensor s = V[i][0];
            double dist = get_distance(s, depot);
            double energy_flying = dist * ecf;

            double energy_hovering = compute_energy_hovering(s);

            double total_energy = energy_flying * 2 + energy_hovering;
            tmp.tours_costs.push_back(total_energy);

            vector<tuple<point, int>> vec;
            vec.emplace_back(depot, -1);
            vec.emplace_back(make_tuple(s.get_pos_x(), s.get_pos_y()), 0);
            vec.emplace_back(depot, -1);

            tmp.tours.emplace_back(vec);

        } else {
            auto nei_sol = appro_alg_nei(V[i], i, depot, radius);
            tmp.tours = nei_sol.tours;
            tmp.tours_costs = nei_sol.tours_costs;
        }

        for (int k = 0; k < tmp.tours.size(); k++) {
            sol.tours.push_back(tmp.tours[k]);
            sol.tours_costs.push_back(tmp.tours_costs[k]);
        }
    }

    return sol;
}

solution algorithms::appro_alg_nei(vector<sensor> V, int jth, point depot, double radius) {
    double budget = dep->get_energy_budget();
    double epsilon = dep->get_epsilon();
    double energy_budget = pow(2, jth - 1) * epsilon * budget;

    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    double ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_DTR(0); // constant DTR (in MB/s)
    //double radius = dep->get_sensor_radius();

    vector<vector<tuple<point, int>>> sol_tours;
    vector<double> sol_costs;
    for (const auto& s: V) {
        vector<tuple<point, int>> A;
        A.emplace_back(make_tuple(s.get_pos_x(), s.get_pos_y()), -1);
        sol_tours.emplace_back(A);
        sol_costs.emplace_back(energy_budget);
    }

    auto n = V.size();

    // G', n x n, pairwise weight between points
    vector<vector<double>> G1(n, vector<double>(n));
    for (int i = 2; i <= n; i++) {
        double delta = energy_budget / i;
        // Algorithm 4 in [27]
        for (int j = 0; j < n; j++) {
            sensor sensor_j = V[j];
            double energy_j_hovering = compute_energy_hovering(sensor_j);

            for (int k = 0; k < n; k++) {
                // compute C(j, k)
                sensor sensor_k = V[k];
                double dist = get_distance(sensor_j, sensor_k);
                double energy_k_hovering = compute_energy_hovering(sensor_k);

                if (dist <= 2 * radius) {   // R_0
                    double w1 = energy_j_hovering + energy_k_hovering;
                    // G''
                    if (w1 <= delta) {
                        G1[j][k] = w1;
                    } else {
                        G1[j][k] = 0.;
                    }
                } else {
                    double w2 = (dist - 2 * radius) * ecf + energy_j_hovering + energy_k_hovering;  // R_0
                    // G''
                    if (w2 <= delta) {
                        G1[j][k] = w2;
                    } else {
                        G1[j][k] = 0.;
                    }
                }
            }
        }
        // find connected components
        vector<vector<int>> adjacencyList(n);

        for (int i_adj = 0; i_adj < n; i_adj++) {
            for (int j_adj = 0; j_adj < n; j_adj++) {
                if (G1[i_adj][j_adj] > 0) {
                    adjacencyList[i_adj].push_back(j_adj);
                }
            }
        }

        vector<vector<sensor>> components;
        unordered_set<int> visited;
        int p = 0;
        for (int v = 0; v < n; ++v) {
            if (visited.find(v) == visited.end()) {
                unordered_set<int> connectedComponent;
                DFS(v, visited, connectedComponent, adjacencyList);
                components.emplace_back();
                for (int u: connectedComponent) {
                    components[p].push_back(V[u]);
                }
                p++;
            }
        }

        vector<vector<tuple<point, int>>> tours1;
        vector<double> costs1;
        // for each component run tspn
        for (auto &component: components) {
            auto res = tsp_neighbors_v1(component, radius);
            res = tsp_split(res.tours[0], res.tours_costs, depot, component, true);

            for (int w = 0; w < res.tours.size(); w++) {
                tours1.push_back(res.tours[w]);
                costs1.push_back(res.tours_costs[w]);
            }
        }
        if (tours1.size() < sol_tours.size()) {
            sol_tours = tours1;
            sol_costs = costs1;
        }
    }

    solution sol;
    sol.tours = sol_tours;
    sol.tours_costs = sol_costs;

    return sol;
}

void algorithms::DFS(int v, unordered_set<int> &visited, unordered_set<int> &connectedComponent,
                     vector<vector<int>> adjacencyList) {
    visited.insert(v);
    connectedComponent.insert(v);
    for (int neighbor: adjacencyList[v]) {
        if (visited.find(neighbor) == visited.end()) {
            DFS(neighbor, visited, connectedComponent, adjacencyList);
        }
    }
}

solution algorithms::tsp_neighbors_v1(const vector<sensor> &sensors, double radius) {
    vector<tuple<point, int>> tspn_result;
    vector<double> tspn_cost;

    vector<sensor> deployed_sensors = sensors;

    vector<sensor> orig_sensors = dep->get_sensors();
    // Sort sensors by x-coordinate
    sort(deployed_sensors.begin(), deployed_sensors.end(), [](const sensor &a, const sensor &b) {
        return a.get_pos_x() < b.get_pos_x();
    });

    map<int, map<int, vector<point>>> intersections;
    map<int, int> checked_sensors;

    for (int i = 0; i < deployed_sensors.size(); i++) {
        if (checked_sensors[deployed_sensors[i].get_id()] == 0) {
            point p1 = deployed_sensors[i].get_position();
            int intersect = 0;
            for (int j = i + 1; j < deployed_sensors.size(); j++) {
                if (checked_sensors[deployed_sensors[j].get_id()] == 0) {
                    // find intersect of i and j
                    point p2 = deployed_sensors[j].get_position();
                    vector<point> intersec_points = get_circles_intersections(p1, p2, radius);
                    if (!intersec_points.empty()) {
                        intersections[deployed_sensors[i].get_id()][deployed_sensors[j].get_id()] = intersec_points;
                        intersect = 1;
                        checked_sensors[deployed_sensors[j].get_id()] = 1;
                    }
                }
            }
            checked_sensors[deployed_sensors[i].get_id()] = 1;
            if (intersect == 0) {
                vector<point> intersec_points = {make_tuple(-1, -1)};
                intersections[deployed_sensors[i].get_id()][deployed_sensors[i].get_id()] = intersec_points;
            }
        }
    }

    // cout << "******* " << endl;
    // for (auto t : intersections){
    //     cout << t.first << endl;
    //     for (auto p : t.second){
    //         cout << p.first << " : " ;
    //         for (auto q : p.second){
    //             cout << get<0>(q) << " " << get<1>(q) << endl;
    //         }
    //     }
    //     cout << "_________" << endl;
    // }

    // get the keys of map (independent set) and compute tsp
    vector<point_2d> points;
    for (const auto &p: intersections) {
        point_2d new_point = {orig_sensors[p.first].get_pos_x(), orig_sensors[p.first].get_pos_y(), p.first};
        points.push_back(new_point);
    }

    TSP tsp(points);
    tsp.solve();

    vector<int> tsp_result_id = tsp.get_path_id();
    //cout << "tsp_result_id " << tsp_result_id.size() << endl;
    // tsp_result = tsp.get_path();
    // cout << "TSP path: ";
    // for (auto p : tsp_result_id) {
    //     cout << p << ", ";
    // }
    // cout << endl;

    // for each tsp_result_id, solve tsp for intersections[orig_sensprs[tsp_result_id[i]]]
    // and create the final tour

    int counter = 0;
    for (int i = 0; i < tsp_result_id.size(); i++) {
        int id = points[tsp_result_id[i]].id;
        sensor s1 = orig_sensors[id];

        counter++;
        if (counter == tsp_result_id.size()) {
            tspn_result.emplace_back(make_tuple(s1.get_pos_x(), s1.get_pos_y()), id);
            break;
        }

        tspn_result.emplace_back(make_tuple(s1.get_pos_x(), s1.get_pos_y()), id);

        int next_id = points[tsp_result_id[i + 1]].id;
        sensor s2 = orig_sensors[next_id];

        // for intersection points in intersections[id] solve tsp
        vector<point_2d> intersection_points;
        for (const auto &p: intersections[id]) {
            for (const auto &pp: p.second) {
                if (p.first != id) {
                    point_2d new_point = {get<0>(pp), get<1>(pp), p.first};
                    intersection_points.push_back(new_point);
                }
            }
        }

        if (!intersection_points.empty()) {
            TSP tsp_new(intersection_points);
            tsp_new.solve();

            vector<int> tsp_result_id_intersection = tsp_new.get_path_id();

            // add new points using tsp_result_id_intersection to tspn_result
            map<int, int> checked_ids;
            for (int j: tsp_result_id_intersection) {
                point_2d p = intersection_points[j];
                if (checked_ids[p.id] == 0) {
                    tspn_result.emplace_back(make_tuple(p.x, p.y), p.id);
                    checked_ids[p.id] = 1;
                }
            }
        }
    }

    // cout << "tsp : " << endl;
    // for (int i = 0; i < tspn_result.size(); i++) {
    //     auto [x, y] = get<0>(tspn_result[i]);  // Structured binding
    //     cout << "(" << x << ", " << y << ")" << " : " << get<1>(tspn_result[i]) << endl;
    // }
    // cout << endl;

    // compute the cost
    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)

    for (int i = 0; i < tspn_result.size() - 1; i++) {
        sensor s1 = orig_sensors[get<1>(tspn_result[i])];
        sensor s2 = orig_sensors[get<1>(tspn_result[i + 1])];
        //double dist = get_distance(s1, s2);
        double dist = get_distance(get<0>(tspn_result[i]), get<0>(tspn_result[i+1]));
        double energy_flying = dist * ecf;

        double energy_hovering_s1 = compute_energy_hovering(s1);
        double energy_hovering_s2 = compute_energy_hovering(s2);

        double total_energy = energy_flying + energy_hovering_s1 / 2. + energy_hovering_s2 / 2.;
        tspn_cost.push_back(total_energy);
    }

    // vector<vector<tuple<point, int>>> tspn_tours;
    // tspn_tours.push_back(tspn_result);
    // draw_result(tspn_tours, true);

    solution sol;
    sol.tours.push_back(tspn_result);
    sol.tours_costs = tspn_cost;

    return sol;
}


//////////////////////////////////////
solution algorithms::tsp_neighbors_v2(const vector<sensor> &sensors, double radius) {
    vector<tuple<point, int>> tspn_result;
    vector<double> tspn_cost;

    vector<sensor> deployed_sensors = sensors;

    vector<sensor> orig_sensors = dep->get_sensors();
    // Sort sensors by x-coordinate
    sort(deployed_sensors.begin(), deployed_sensors.end(), [](const sensor &a, const sensor &b) {
        return a.get_pos_x() < b.get_pos_x();
    });

    // for (auto s:deployed_sensors)
    // {
    //     cout << s.get_pos_x() << ", " << s.get_pos_y() << endl;
    // }
    

    // Step 1: TSP Christofides on all vertices(deployed_sensors)
    vector<point_2d> points;
    for (const auto &s: deployed_sensors) {
        // find id of every s in orig_sensors
        int id = distance(orig_sensors.begin(), find(orig_sensors.begin(), orig_sensors.end(), s));
        //cout << "id   " << id << endl;
        point_2d new_point = {orig_sensors[id].get_pos_x(), orig_sensors[id].get_pos_y(), id};
        points.push_back(new_point);
    }

    TSP tsp(points);
    tsp.solve();

    vector<int> tsp_result_id = tsp.get_path_id();

    // for (auto t : tsp_result_id)
    // {
    //     cout << "id  " << t << endl;
    // }
    

    // Step 2: for any consecutive vertices u->v, see the intersections among line u->v, and the circle centered in v
    //    point u, v;
    //    vector<point> result = get_line_circle_intersections(u, v, radius);

    if (radius > 0){
       
        point point_u;
        int id_u;
        int id_v;
        for (int i = 0; i < tsp_result_id.size(); i++) {

            if (i == tsp_result_id.size() - 1){
                // add last point to the begining of tspn_result
                tuple<point, int> p = {point_u, id_v};
                tspn_result.insert(tspn_result.begin(), p);
                break;
            }

            // update u to the new obtained point from its intersection with v
            // v is sensor, u is a point with id_u

            // if (i == tsp_result_id.size() - 2){  // no need
            //     id_u = points[tsp_result_id[i]].id;
            //     id_v = points[tsp_result_id[0]].id;
            // } else {
            id_u = points[tsp_result_id[i]].id;
            id_v = points[tsp_result_id[i+1]].id;
            // }
            
            if (i == 0){
                point_u = {orig_sensors[id_u].get_pos_x(), orig_sensors[id_u].get_pos_y()};
            }
            
            sensor v = orig_sensors[id_v];
            
            vector<point> result = get_line_circle_intersections(point_u, make_tuple(v.get_pos_x(), v.get_pos_y()), radius);
            // if u is inside v, select the center of v
            double dist_u_v = get_distance(v, point_u); 

            if (dist_u_v <= radius){
                tspn_result.emplace_back(make_tuple(v.get_pos_x(), v.get_pos_y()), id_v);
                point_u = make_tuple(v.get_pos_x(), v.get_pos_y());
            } else {
                // select the intersection point closest to u
                double dist1 = get_distance(point_u, result[0]);
                double dist2 = get_distance(point_u, result[1]);
                if (dist1 < dist2){
                    tspn_result.emplace_back(result[0], id_v);
                    point_u = result[0];
                } else {
                    tspn_result.emplace_back(result[1], id_v);
                    point_u = result[1];
                }
            }
        }

        // Step 3: for any consecutive vertices u->v->w, see the intersections among line u->w, and the circle centered in v
        //    point u, v, w;
        //    vector<point> result = get_line_circle_intersections(u, v, radius, w);

        for (int j = 0; j < tspn_result.size(); j++){    // j < tspn_result.size() - ...      ???
            point a = get<0>(tspn_result[j]);
            point b = get<0>(tspn_result[j+1]);
            point c = get<0>(tspn_result[j+2]);

            vector<point> result2 = get_line_circle_intersections(a, b, radius, c);
            
            double dist_ab = get_distance(a, b);
            double dist_bc = get_distance(b, c);
            double dist_abc = dist_ab + dist_bc;
            double min_dist = dist_abc;
            point best_p;
            // select best p, from a to c
            for (auto p: result2){
                double dist_ap = get_distance(a, p);
                double dist_pc = get_distance(p, c);
                double dist_apc = dist_ap + dist_pc;
                if (dist_apc < min_dist){
                    min_dist = dist_apc;
                    best_p = p;
                }
            }
            
            /////// CONTINUE
        }
        


     

    } else {
        for (int i = 0; i < tsp_result_id.size(); i++) {
            int id = points[tsp_result_id[i]].id;
            sensor s = orig_sensors[id];
            tspn_result.emplace_back(make_tuple(s.get_pos_x(), s.get_pos_y()), id);
        }

    }

    // cout << "tspn_result.size()  " << tspn_result.size() << endl;
    // for (size_t j = 0; j < tspn_result.size(); j++){
    //     cout << "tsp " << get<0>(get<0>(tspn_result[j])) << ", " << get<1>(get<0>(tspn_result[j])) << " : " << get<1>(tspn_result[j]) <<  endl;
    // }
    


   
    

    // compute the cost
    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    double tot_flying = 0.;                ////////////////////////////////////
    double tot_hovering = 0.;        ////////////////////////////
    for (int i = 0; i < tspn_result.size() - 1; i++) {
        sensor s1 = orig_sensors[get<1>(tspn_result[i])];
        sensor s2 = orig_sensors[get<1>(tspn_result[i + 1])];
        double dist = get_distance(get<0>(tspn_result[i]), get<0>(tspn_result[i+1]));
        double energy_flying = dist * ecf;
        tot_flying = tot_flying + energy_flying;

        double energy_hovering_s1 = compute_energy_hovering(s1);
        tot_hovering = tot_hovering + energy_hovering_s1;
        double energy_hovering_s2 = compute_energy_hovering(s2);

        double total_energy = energy_flying + energy_hovering_s1 / 2. + energy_hovering_s2 / 2.;
        tspn_cost.push_back(total_energy);
    }
     //cout << "Tot flyingg  " << tot_flying << endl;
     //cout << "Tot hovering  " << tot_hovering << endl;



    // Step 3: for any consecutive vertices u->v->w, see the intersections among line u->w, and the circle centered in v
    //    point u, v, w;
    //    vector<point> result = get_line_circle_intersections(u, v, radius, w);


//        // Example usage
//        point pa = {-6, 2};
//        point pb = {-2.5, 2.5};
//        point pc = {1, 0.5};
//        double radius = 2;
//
//        auto points = get_line_circle_intersections(pa, pb, radius);
//        // Output the intersection points
//        for (const auto& intersection : points) {
//            double x = get<0>(intersection);
//            double y = get<1>(intersection);
//            cout << "Intersection Point: (" << x << ", " << y << ")\n";
//        }
//        cout << "----" << endl;
//        points = get_line_circle_intersections(pa, pb, radius, pc);
//        // Output the intersection points
//        for (const auto& intersection : points) {
//            double x = get<0>(intersection);
//            double y = get<1>(intersection);
//            cout << "Intersection Point: (" << x << ", " << y << ")\n";
//        }
//        exit(1);

    solution sol;
    sol.tours.push_back(tspn_result);
    sol.tours_costs = tspn_cost;

    return sol;
}

solution algorithms::tsp_split(vector<tuple<point, int>> tspn_result, const vector<double> &tspn_cost, point depot,
                               const vector<sensor> &sensors, bool violation) {
    vector<vector<tuple<point, int>>> sol_tours;
    vector<double> sol_costs;
    double energy_budget = dep->get_energy_budget();
    double epsilon = dep->get_epsilon();
    if (violation) {
        energy_budget = (1 + epsilon) * energy_budget;
    }

    int i = 0;
    int j = 0;
    auto n = sensors.size(); // dep->get_sensors().size();

    double within_budget = -1;
    while (j <= n - 1) {
        vector<tuple<point, int>> T_k;

        T_k.emplace_back(depot, -1);
        for (int p = i; p <= j; p++) {
            T_k.push_back(tspn_result[p]);
        }
        T_k.emplace_back(depot, -1);
        double cost = tour_cost(T_k, tspn_cost, i, j, depot, sensors);
        if (cost <= energy_budget) {
            within_budget = cost;
            if (j == n - 1) {
                sol_tours.push_back(T_k);
                sol_costs.push_back(within_budget);
                break;
            } else {
                j++;
            }
        } else {
            j--;
            // remove j from T_k
            T_k.erase(T_k.end() - 2);  // Remove the second-to-last element

            sol_tours.push_back(T_k);
            sol_costs.push_back(within_budget);
            T_k.clear();

            j++;
            i = j;
            within_budget = -1;
        }
    }

    solution sol;
    sol.tours = sol_tours;
    sol.tours_costs = sol_costs;

    return sol;
}

// given two points pa and pb, returns the intersection points
vector<point> algorithms::get_circles_intersections(const point &pa, const point &pb, const double radius) {
    point int_ab_1 = {-1, -1};
    point int_ab_2 = {-1, -1};
    vector<point> int_points;

    double x1 = get<0>(pa);
    double y1 = get<1>(pa);
    double x2 = get<0>(pb);
    double y2 = get<1>(pb);

    // Calculate distance between the centers of the circles
    double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

    // Check if circles are too far apart or coincide
    if (distance > 2 * radius || distance < fabs(radius - radius)) {
        // cout << "No intersection points found." << endl;
        return int_points;
    }

    // Calculate the distance from the center of circle A to the line joining the intersection points
    double a = (pow(radius, 2) - pow(radius, 2) + pow(distance, 2)) / (2 * distance);

    // Calculate the coordinates of the intersection points
    double x3 = x1 + a * (x2 - x1) / distance;
    double y3 = y1 + a * (y2 - y1) / distance;

    double h = sqrt(pow(radius, 2) - pow(a, 2));

    // Calculate the coordinates of the two intersection points
    double int_x1 = x3 + h * (y2 - y1) / distance;
    double int_y1 = y3 - h * (x2 - x1) / distance;

    double int_x2 = x3 - h * (y2 - y1) / distance;
    double int_y2 = y3 + h * (x2 - x1) / distance;

    int_ab_1 = {int_x1, int_y1};
    int_ab_2 = {int_x2, int_y2};

    int_points.push_back(int_ab_1);
    int_points.push_back(int_ab_2);

    return int_points;
}

// "pa" is a point, "pb" is the center of a circle of radius "radius".
// Return the intersection points between the line that passes through "pa" and "pb" which intersects the circle centered in "pb"
vector<point> algorithms::get_line_circle_intersections(const point& pa, const point& pb, double radius) {
    return get_line_circle_intersections_helper(pa, pb, radius, pb);
}

// "pa" is a point, "pb" is the center of a circle of radius "radius", "pc" is another point
// Return the intersection points between the line that passes through "pa" and "pc" which intersects the circle centered in "pb"
vector<point> algorithms::get_line_circle_intersections(const point& pa, const point& pb, double radius, const point& pc) {
    return get_line_circle_intersections_helper(pa, pb, radius, pc);
}

vector<point> algorithms::get_line_circle_intersections_helper(const point &pa, const point &pb, double radius, const point &pc) {
    // Extract coordinates from tuples
    double xa, ya, xb, yb, xc, yc;
    tie(xa, ya) = pa;
    tie(xb, yb) = pb;
    tie(xc, yc) = pc;

    // Helper function to check equality with epsilon
    auto are_equal = [](double a, double b) {
        return fabs(a - b) < numeric_limits<double>::epsilon();
    };

    // Check if the line PaPc is vertical
    if (are_equal(xa, xc)) {
        // For a vertical line, calculate the x-coordinate of the line
        double x = xa;

        // Calculate the y-coordinates of the intersection points
        double y1 = yb + sqrt(radius * radius - (x - xb) * (x - xb));
        double y2 = yb - sqrt(radius * radius - (x - xb) * (x - xb));

        // Add the intersection points to the vector
        return {make_tuple(x, y1), make_tuple(x, y2)};
    }

    // Calculate the slope of the line PaPc
    double slope_pc = (yc - ya) / (xc - xa);

    // Calculate the y-intercept of the line PaPc
    double intercept_pc = ya - slope_pc * xa;

    // Calculate coefficients for the quadratic equation
    double a = 1 + slope_pc * slope_pc;
    double b = -2 * xb + 2 * slope_pc * intercept_pc - 2 * yb * slope_pc;
    double c = xb * xb + yb * yb + intercept_pc * intercept_pc - 2 * intercept_pc * yb - radius * radius;

    // Calculate the discriminant
    double discriminant = b * b - 4 * a * c;

    vector<point> intersections;

    if (discriminant >= 0) {
        // Calculate the x-coordinates of intersection points
        double x1 = (-b + sqrt(discriminant)) / (2 * a);
        double x2 = (-b - sqrt(discriminant)) / (2 * a);

        // Calculate the corresponding y-coordinates
        double y1 = slope_pc * x1 + intercept_pc;
        double y2 = slope_pc * x2 + intercept_pc;

        // Add the intersection points to the vector
        intersections.emplace_back(x1, y1);
        intersections.emplace_back(x2, y2);
    }

    return intersections;
}

void algorithms::draw_result(vector<vector<tuple<point, int>>> tspn_tours, bool single, bool doi) {
    ofstream htmlFile("html/sensor_deployment.html");
    //ofstream htmlFile("../output/sensor_deployment.html");



    htmlFile << "<!DOCTYPE html>\n<html>\n<head>\n";
    htmlFile << "<title>Sensor Deployment</title>\n";
    htmlFile << "</head>\n<body>\n";

    htmlFile << "<script>\n";
    htmlFile << "function drawSensorsAndDepots() {\n";
    htmlFile << "var canvas = document.getElementById('sensorCanvas');\n";
    htmlFile << "var ctx = canvas.getContext('2d');\n";

    // Draw label for point (0, 0) with an offset of +15, +15
    htmlFile << "ctx.font = 'bold 15px Arial';\n";
    htmlFile << "ctx.fillText('(0, 0)', 15, " << dep->get_area_width() - 15 << ");\n";

    // Draw (0, 0)
    htmlFile << "ctx.beginPath();\n";
    htmlFile << "ctx.arc(" << 0 << ", " << dep->get_area_width() << ", 5, 0, 2 * Math.PI);\n";
    htmlFile << "ctx.fillStyle = 'black';\n";
    htmlFile << "ctx.fill();\n";
    htmlFile << "ctx.stroke();\n";

    // Draw grid lines
    htmlFile << "ctx.strokeStyle = 'rgba(0, 0, 0, 0.1)';\n";
    htmlFile << "ctx.lineWidth = 1;\n";
    for (int x = 100; x < dep->get_area_length(); x += 100) {
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.moveTo(" << x << ", 0);\n";
        htmlFile << "ctx.lineTo(" << x << ", " << dep->get_area_width() << ");\n";
        htmlFile << "ctx.stroke();\n";
    }

    for (int y = 100; y < dep->get_area_width(); y += 100) {
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.moveTo(0, " << y << ");\n";
        htmlFile << "ctx.lineTo(" << dep->get_area_length() << ", " << y << ");\n";
        htmlFile << "ctx.stroke();\n";
    }

    // Draw sensors
    vector<sensor> sensors = dep->get_sensors();
    for (size_t i = 0; i < sensors.size(); ++i) {
        auto s = sensors[i];
        auto pos = s.get_position();

        double sensor_x = get<0>(pos);
        double sensor_y = get<1>(pos);

        if (doi) {
            // Draw sensor shape based on radius_doi
            htmlFile << "ctx.beginPath();\n";
            for (int angle = 0; angle < 360; ++angle) {
                double radius = s.get_radius_doi(angle);
                double x = sensor_x + radius * cos(angle * M_PI / 180.0);
                double y = dep->get_area_width() - sensor_y + radius * sin(angle * M_PI / 180.0);

                if (angle == 0) {
                    htmlFile << "ctx.moveTo(" << x << ", " << y << ");\n";
                } else {
                    htmlFile << "ctx.lineTo(" << x << ", " << y << ");\n";
                }
            }
            htmlFile << "ctx.closePath();\n";
            htmlFile << "ctx.fillStyle = 'rgba(255, 165, 0, 0.15)';\n";
            htmlFile << "ctx.fill();\n";
            htmlFile << "ctx.stroke();\n";
        } else {
            // Draw sensor circle
            htmlFile << "ctx.beginPath();\n";
            htmlFile << "ctx.arc(" << sensor_x << ", " << dep->get_area_width() - sensor_y << ", "
                     << dep->get_sensor_radius()
                     << ", 0, 2 * Math.PI);\n";
            htmlFile << "ctx.fillStyle = 'rgba(255, 165, 0, 0.15)';\n";
            htmlFile << "ctx.fill();\n";
            htmlFile << "ctx.stroke();\n";
        }

        // Draw sensor center dot
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.arc(" << sensor_x << ", " << dep->get_area_width() - sensor_y << ", 2, 0, 2 * Math.PI);\n";
        htmlFile << "ctx.fillStyle = 'black';\n";
        htmlFile << "ctx.fill();\n";
        htmlFile << "ctx.stroke();\n";

//        // Draw sensor label with coordinates (bold)
//        htmlFile << "ctx.fillStyle = 'black';\n";
//        htmlFile << "ctx.font = 'bold 15px Arial';\n";
//        htmlFile << "ctx.fillText('" << i << " (" << sensor_x << ", " << sensor_y << ")', "
//                 << sensor_x + 10 << ", " << dep->get_area_width() - sensor_y + 10 << ");\n";
    }

    // Draw depots
    for (size_t i = 0; i < dep->get_depots().size(); ++i) {
        auto d = dep->get_depots()[i];
        auto pos = d;
        double depot_x = get<0>(pos);
        double depot_y = get<1>(pos);

        // Draw depot square
        htmlFile << "ctx.fillStyle = 'green';\n";
        htmlFile << "ctx.fillRect(" << depot_x - 7.5 << ", " << dep->get_area_width() - depot_y - 7.5 << ", 15, 15);\n";

//        // Draw depot label with coordinates
//        htmlFile << "ctx.fillStyle = 'black';\n";
//        htmlFile << "ctx.font = '15px Arial';\n";
//        htmlFile << "ctx.fillText('D" << i << " (" << depot_x << ", " << depot_y << ")', "
//                 << depot_x + 10 << ", " << dep->get_area_width() - depot_y + 10 << ");\n";

        if (single) {
            break;
        }
    }

    // Draw TSP circuit connecting points in the order specified by tspn_tours
    vector<string> path_colors = {"red", "blue", "green", "orange", "purple", "cyan", "magenta",
                                  "yellow", "pink", "brown", "teal", "olive", "maroon", "navy", "lavender", "turquoise", "indigo"};
    int col = 0;
    for (auto &tspn_tour: tspn_tours) {
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.lineWidth = 2;\n";
        htmlFile << "ctx.strokeStyle = '" << path_colors[col % path_colors.size()] << "';\n";
        htmlFile << "ctx.setLineDash([5, 5]);\n"; // 5-pixel dashes with 5-pixel gaps
        for (auto j: tspn_tour) {
            double path_x = get<0>(get<0>(j));
            double path_y = get<1>(get<0>(j));
            htmlFile << "ctx.lineTo(" << path_x << ", " << dep->get_area_width() - path_y << ");\n";
        }
        htmlFile << "ctx.stroke();\n";

        for (auto j: tspn_tour) {
            double path_x = get<0>(get<0>(j));
            double path_y = get<1>(get<0>(j));

            // Draw a circle at each coordinate
            htmlFile << "ctx.beginPath();\n";
            htmlFile << "ctx.arc(" << path_x << ", " << dep->get_area_width() - path_y << ", 3, 0, 2 * Math.PI);\n";
            htmlFile << "ctx.fillStyle = '" << path_colors[col % path_colors.size()] << "';\n";
            htmlFile << "ctx.fill();\n";
            htmlFile << "ctx.closePath();\n";
            htmlFile << "ctx.stroke();\n";
        }

        col++;
    }

    htmlFile << "}\n";
    htmlFile << "</script>\n";

    htmlFile << "<canvas id='sensorCanvas' width='" << dep->get_area_length() << "' height='" << dep->get_area_width()
             << "' style='border:1px solid #000;'></canvas>\n";

    htmlFile << "<script>\n";
    htmlFile << "drawSensorsAndDepots();\n";
    htmlFile << "</script>\n";

    htmlFile << "</body>\n</html>";

    htmlFile.close();
}

