#include "algorithms.h"

#include <unordered_set>

algorithms::algorithms(deployment *m_dep) {
    dep = m_dep;
}

void algorithms::approxTSPN_S() {
    cout << "approxTSPN_S: " << endl;

    auto [tspn_result, tspn_cost] = tsp_neighbors(dep->get_sensors());
    point depot = dep->get_depots()[0];
    auto our_tours = tsp_split(tspn_result, tspn_cost, dep->get_energy_budget(), depot);

    cout << "Number of drones: " << our_tours.size() << endl << endl;
//    for (const auto& r : our_tours) {
//        for (auto azz : r) {
//            cout << get<0>(azz) << ", " << get<1>(azz) << ", " << get<2>(azz) << endl;
//        }
//        cout << endl;
//    }
    draw_result(our_tours);
}

/*void algorithms::approxMPN_S() {
    cout << "approxMPN_S: " << endl;
    point depot = dep->get_depots()[0];
    approxMPN(depot);
    cout << "Number of drones: " << tspn_tours.size() << endl;
    draw_result();
}*/

/*void algorithms::approxTSPN_M() {
    cout << "approxTSPN_M: " << endl;
    double ecf = dep->get_energy_cons_fly();
    double radius = dep->get_sensor_radius();
    double R_0_f = radius * ecf;
    vector<point> depots = dep->get_depots();
    vector<double> A_d(depots.size(), 0.0);

    // find the depot
    for (int i = 0; i < depots.size(); i++) {
        for (const auto &sensor: dep->get_sensors()) {
            double dist = get_distance(sensor, depots[i]);
            auto pos = sensor.get_position();
            double energy_hovering = compute_energy_hovering(
                    make_tuple(get<0>(pos), get<1>(pos), sensor.get_data_size(), 0));
            if (dist <= radius) {
                if (energy_hovering > A_d[i]) {
                    A_d[i] = energy_hovering;
                }
            } else {
                double a = 2 * (dist * ecf + (energy_hovering / 2) + R_0_f);
                if (a > A_d[i]) {
                    A_d[i] = a;
                }
            }
        }
    }

    auto index = distance(A_d.begin(), min_element(A_d.begin(), A_d.end()));
    point depot = depots[index];

    tsp_neighbors(dep->get_sensors());
    tsp_split(dep->get_energy_budget(), depot);
    cout << "Number of drones: " << tspn_tours.size() << endl;
    draw_result();
}*/

/*void algorithms::approxMPN_M() {
    cout << "approxMPN_M: " << endl;
    vector<point> depots = dep->get_depots();
    vector<sensor> sensors = dep->get_sensors();
    vector<vector<tuple<double, double, int>>> tours;
    vector<tuple<double, double, int>> vec;
    vec.emplace_back(0.0, 0.0, 0);

    for (int i = 0; i < 2 * sensors.size(); i++) {
        tours.push_back(vec);
    }

    for (auto depot: depots) {
        approxMPN(depot);
        if (tspn_tours.size() < tours.size()) {
            tours = tspn_tours;
        }
    }

    tspn_tours = tours;
    cout << "Number of drones: " << tspn_tours.size() << endl;
    draw_result();
}*/

double algorithms::get_distance(point p1, point p2) {
    double x1, y1, x2, y2;
    tie(x1, y1) = p1;
    tie(x2, y2) = p2;

    double delta_x = x2 - x1;
    double delta_y = y2 - y1;

    return sqrt(delta_x * delta_x + delta_y * delta_y);
}

double algorithms::get_distance(const sensor& s, point p) {
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

int algorithms::get_angle(const sensor& s, point p) {
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

double algorithms::compute_energy_hovering(tuple<double, double, double, double> sensor) {
    double dtr = dep->get_data_transfer_rate();
    double ech = dep->get_energy_cons_hover();
    double sensor_data = get<2>(sensor);
    double required_time = sensor_data / dtr;
    return required_time * ech;
}

double algorithms::compute_energy_hovering(sensor s) {
    double dtr = dep->get_data_transfer_rate();
    double ech = dep->get_energy_cons_hover();
    double sensor_data = s.get_data_size();
    double required_time = sensor_data / dtr;
    return required_time * ech;
}

double algorithms::tour_cost(vector<tuple<point, int>> T, vector<double> tspn_cost, int start, int end, point depot) {
    vector<sensor> deployed_sensors = dep->get_sensors();

    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    double ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)

    double cost_T_k = 0;

    if (get<0>(depot) != -1) {
        auto lastIndex = T.size() - 2;

        point p0 = get<0>(T[0]);
        point p1 = get<0>(T[1]);
        point pLast = get<0>(T[lastIndex]);

        double dist1 = get_distance(p0, p1);
        double energy1_flying = dist1 * ecf;

        sensor sensor1 = deployed_sensors[get<1>(T[1])];
        double energy1_hovering = compute_energy_hovering(sensor1);

        double dist2 = get_distance(p0, pLast);
        double energy2_flying = dist2 * ecf;

        sensor sensor2 = deployed_sensors[get<1>(T[lastIndex])];
        double energy2_hovering = compute_energy_hovering(sensor2);

        for (int p = start; p < end; p++) {
            cost_T_k = cost_T_k + tspn_cost[p];
        }
        cost_T_k = cost_T_k + energy1_flying + energy2_flying + energy1_hovering / 2 + energy2_hovering / 2;

    } else {
        for (int p = start; p < end; p++) {
            cost_T_k += tspn_cost[p];
        }
    }

    return cost_T_k;
}

/*void algorithms::approxMPN(point depot) {
    sorted_sensors.clear();
    for (const auto &sensor: dep->get_sensors()) {       // it can be removed
        auto pos = sensor.get_position();
        sorted_sensors.emplace_back(get<0>(pos), get<1>(pos), sensor.get_data_size(), 0);
    }

    // get eps and compute t
    double epsilon = dep->get_epsilon();
    int t = ceil(log2(1 / epsilon));

    double energy_budget = dep->get_energy_budget();

    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    double ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)
    double radius = dep->get_sensor_radius();
    double R_0_f = radius * ecf;                // R_0

    vector<vector<tuple<double, double, double, double>>> V(t + 1);

    // for each s in sorted_sensors compute lambd_v,d and j
    for (auto s: sorted_sensors) {
        double dist = sqrt(pow(get<0>(s) - get<0>(depot), 2) + pow(get<1>(s) - get<1>(depot), 2));
        double energy_flying = dist * ecf;

        double energy_hovering = compute_energy_hovering(s);
        double lambda = energy_flying + energy_hovering / 2. + R_0_f;

        int j = floor(log2((energy_budget - 2. * lambda) / (epsilon * energy_budget)) + 1.);

        // find V_0,..., V_t
        if ((static_cast<double>((1. - epsilon) * (energy_budget / 2.)) < lambda) && (lambda <= static_cast<double>(energy_budget / 2.))) {
            V[0].push_back(s);
        } else {
            V[j].push_back(s);
        }
    }

    vector<vector<tuple<double, double, int>>> tours;
    for (int i = 0; i < V.size(); i++) {
        // for each V[i], run Minimum UAV Deployment Problem with Neighborhoods with budget 2^{j-1} epsilon B
        vector<vector<tuple<double, double, int>>> C;
        C = approAlgNei(V[i], i);

        for (auto & j : C) {
            j.insert(j.begin(), make_tuple(get<0>(depot), get<1>(depot), -1));
            j.emplace_back(get<0>(depot), get<1>(depot), -1);
        }
        for (const auto & j : C) {
            tours.push_back(j);
        }
    }
    tspn_tours.clear();
    tspn_tours = tours;

    // for (auto t : tspn_tours){
    //     for (auto i : t){
    //         cout << get<0>(i) << ", " << get<1>(i) << endl;
    //     }
    //     cout << "______________ " << endl;
    // }
}*/

/*vector<vector<tuple<double, double, int>>> algorithms::approAlgNei(vector<tuple<double, double, double, double>> V, int jth) {
    int total_budget = dep->get_energy_budget();
    double epsilon = dep->get_epsilon();
    double energy_budget = pow(2, jth - 1) * epsilon * total_budget;

    int ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    int ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)
    double radius = dep->get_sensor_radius();

    vector<vector<tuple<double, double, int>>> C;
    for (auto s: V) {
        vector<tuple<double, double, int>> A;
        A.emplace_back(get<0>(s), get<1>(s), -1);
        C.emplace_back(A);
    }

    sort(V.begin(), V.end());
    auto n = V.size();

    // G', n x n, pairwise weight between points
    vector<vector<double>> G1(n, vector<double>(n));

    for (int i = 2; i <= n; i++) {
        double delta = energy_budget / i;
        // Algorithm 4 in [27]
        for (int j = 0; j < n; j++) {
            tuple<double, double, double, double> sensor_j;
            sensor_j = V[j];
            double energy_j_hovering = compute_energy_hovering(sensor_j);

            for (int k = 0; k < n; k++) {
                // compute C(j, k)
                tuple<double, double, double, double> sensor_k;
                sensor_k = V[k];
                double dist = sqrt(
                        pow(get<0>(sensor_j) - get<0>(sensor_k), 2) + pow(get<1>(sensor_j) - get<1>(sensor_k), 2));

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

        vector<vector<tuple<double, double, double, double>>> components;
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

        vector<vector<tuple<double, double, int>>> C1;
        // for each component run tspn
        for (auto & component : components) {
            vector<sensor> comp_sensors;
            for (auto s: component) {
                comp_sensors.emplace_back(get<0>(s), get<1>(s), get<2>(s), vector<double>());
            }
            tsp_neighbors(comp_sensors);
            tsp_split(energy_budget, make_tuple(-1, -1));
            for (const auto & tspn_tour : tspn_tours) {
                C1.push_back(tspn_tour);
            }
        }
        if (C1.size() < C.size()) {
            C = C1;
        }
    }
    return C;
}*/

void algorithms::DFS(int v, unordered_set<int> &visited, unordered_set<int> &connectedComponent, vector<vector<int>> adjacencyList) {
    visited.insert(v);
    connectedComponent.insert(v);
    for (int neighbor: adjacencyList[v]) {
        if (visited.find(neighbor) == visited.end()) {
            DFS(neighbor, visited, connectedComponent, adjacencyList);
        }
    }
}

tuple<vector<tuple<point, int>>, vector<double>> algorithms::tsp_neighbors(const vector<sensor>& sensors) {
    vector<tuple<point, int>> tspn_result;
    vector<double> tspn_cost;

    vector<sensor> deployed_sensors = dep->get_sensors();

    map<int, map<int, vector<point>>> intersections;
    vector<int> checked_sensors(deployed_sensors.size());

    for (int i = 0; i < deployed_sensors.size(); i++) {
        if (checked_sensors[i] == 0) {
            point p1 = deployed_sensors[i].get_position();
            int intersect = 0;
            for (int j = i + 1; j < deployed_sensors.size(); j++) {
                if (checked_sensors[j] == 0) {
                    // add another if to check whether the distance between sensor i and j is less than 2R
                    // find intersect of i and j
                    point p2 = deployed_sensors[j].get_position();
                    vector<point> intersec_points = get_intersection_points(p1, p2);
                    if (!intersec_points.empty()) {
                        intersections[i][j] = intersec_points;
                        intersect = 1;
                        checked_sensors[j] = 1;
                    }
                }
            }
            checked_sensors[i] = 1;
            if (intersect == 0) {
                vector<point> intersec_points = {make_tuple(-1, -1)};
                intersections[i][i] = intersec_points;
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
    vector<int> sensor_ids;
    vector<point_3d> points;
    for (const auto &p: intersections) {
        point_3d new_point = {deployed_sensors[p.first].get_pos_x(), deployed_sensors[p.first].get_pos_y(), 0};
        sensor_ids.push_back(deployed_sensors[p.first].get_id());
        points.push_back(new_point);
    }

    TSP tsp(points);
    tsp.solve();

    vector<int> tsp_result_id = tsp.get_path_id();
    // tsp_result = tsp.get_path();
    // cout << "TSP path: ";
    // for (auto p : tsp_result_id) {
    //     cout << p << ", ";
    // }
    // cout << endl;

    double ecf = dep->get_energy_cons_fly(); // every meter (in J/m)

    int counter = 0;
    for (int i = 0; i < tsp_result_id.size(); i++) {
        int id = sensor_ids[tsp_result_id[i]];
        sensor s1 = deployed_sensors[id];

        counter++;
        if (counter == tsp_result_id.size()) {
            tspn_result.emplace_back(make_tuple(s1.get_pos_x(), s1.get_pos_y()), id);
            break;
        }

        tspn_result.emplace_back(make_tuple(s1.get_pos_x(), s1.get_pos_y()), id);

        int next_id = sensor_ids[tsp_result_id[i+1]];
        sensor s2 = deployed_sensors[next_id];

        int s1_to_s2 = 0;
        auto len_sen1 = intersections[id].size();
        for (const auto& j_intersects : intersections[id]) {
            s1_to_s2++;

            double energy1_hovering = compute_energy_hovering(s1);
            double energy2_hovering = compute_energy_hovering(s2);
            double energy_j_hovering = compute_energy_hovering(deployed_sensors[j_intersects.first]);
            double energy_hovering = energy1_hovering / 2. + energy_j_hovering + energy2_hovering / 2.;

            int skip = 0;
            vector<tuple<double, point>> dist_points;
            for (point p: j_intersects.second) {  // 2 points
                if (get<0>(p) == -1) {
                    skip = 1;
                    break;
                } else {
                    double dist_sen1_p = get_distance(s1, p);
                    double dist_p_sen2 = get_distance(s2, p);

                    double distance = dist_sen1_p + dist_p_sen2;
                    double energy_flying = distance * ecf;
                    double total_energy = energy_hovering + energy_flying;

                    dist_points.emplace_back(total_energy, p);
                }
            }

            // add center of j_intersects.first
            point p_j = deployed_sensors[j_intersects.first].get_position();

            double dist_s1_p = get_distance(s1, p_j);
            double dist_p_s2 = get_distance(s2, p_j);
            double distance = dist_s1_p + dist_p_s2;
            double energy_flying = distance * ecf;

            double total_energy = energy_hovering + energy_flying;

            dist_points.emplace_back(total_energy, p_j);

            if (skip == 1) {
                double dist_s1_s2 = get_distance(s1, s2);
                energy_flying = dist_s1_s2 * ecf;
                energy_hovering = energy1_hovering / 2 + energy2_hovering / 2;
                total_energy = energy_hovering + energy_flying;
                tspn_cost.push_back(total_energy);
                break;
            } else {
                // select the one that has the minimum distance
                sort(dist_points.begin(), dist_points.end());
                point pos = get<1>(dist_points[0]);
                tspn_result.emplace_back(pos, j_intersects.first);

                double dist_s1_pos = get_distance(s1, pos);
                energy_flying = dist_s1_pos * ecf;
                energy_hovering = (energy1_hovering + energy_j_hovering) / 2;
                total_energy = energy_hovering + energy_flying;
                tspn_cost.push_back(total_energy);

                if (s1_to_s2 == len_sen1) {   // or counter == tsp_result_id.size()-1
                    // energy from pos to sen2
                    double dist_pos_s2 = get_distance(s2, pos);
                    energy_flying = dist_pos_s2 * ecf;
                    energy_hovering = (energy2_hovering + energy_j_hovering) / 2;
                    total_energy = energy_hovering + energy_flying;
                    tspn_cost.push_back(total_energy);
                }

                // fake sensor to be replaced to s1 (only to compute distances and energy costs)
                sensor fs1(get<0>(pos), get<1>(pos), deployed_sensors[j_intersects.first].get_data_size(), {});
                s1 = fs1;
            }
        }
    }

    cout << "tsp : " << endl;
    for (int i = 0; i < tspn_result.size(); i++) {
        auto [x, y] = get<0>(tspn_result[i]);  // Structured binding
        cout << "(" << x << ", " << y << ")" << " : " << get<1>(tspn_result[i]) << " : " << tspn_cost[i] << endl;
    }
    cout << endl;

    return make_tuple(tspn_result, tspn_cost);
}

vector<vector<tuple<point, int>>> algorithms::tsp_split(vector<tuple<point, int>> tspn_result, vector<double> tspn_cost, double energy_budget, point depot) {
    vector<vector<tuple<point, int>>> tspn_tours;

    vector<sensor> deployed_sensors = dep->get_sensors();

    int i = 0;
    int j = 0;
    auto n = deployed_sensors.size();

    while (j <= n - 1) {
        vector<tuple<point, int>> T_k;

        if (get<0>(depot) != -1) {
            T_k.emplace_back(depot, -1);
        }
        for (int p = i; p <= j; p++) {
            T_k.push_back(tspn_result[p]);
        }
        if (get<0>(depot) != -1) {
            T_k.emplace_back(depot, -1);
        }

        double cost = tour_cost(T_k, tspn_cost, i, j, depot);

        if (cost <= energy_budget) {
//            cout << "Cost ok!: " << cost << endl;
            if (j == n - 1) {
                tspn_tours.push_back(T_k);
                break;
            } else {
                j++;
            }
        } else {
//            cout << "Cost KO: " << cost << endl;
            j--;
            // remove j from T_k
            if (get<0>(depot) != -1) {
//                T_k.erase(T_k.begin() + T_k.size() - 2);
                T_k.erase(T_k.end() - 2);  // Remove the second-to-last element
            } else {
//                T_k.erase(T_k.begin() + T_k.size() - 1);
                T_k.pop_back();  // Remove the last element
            }

            tspn_tours.push_back(T_k);
            T_k.clear();

            j++;
            i = j;
        }
    }

//     for (int i = 0; i < tspn_tours.size(); i++){
//         for (auto j: tspn_tours[i]){
//             cout<< get<0>(j) << ", " << get<1>(j) << endl;
//         }
//         cout << "-------" << endl;
//     }

    return tspn_tours;
}


// void algorithms::tsp_neighbors() {
//     // sort sensors from left to right

//     point p1 = {0, 0};
//     point p2 = {10, 10};
//     vector<point> aaa = get_intersection_points(p1, p2);
//     for (point p : aaa) {
//         cout << get<0>(p) << ", " << get<1>(p) << endl;
//     }
//     cout << endl;

//     // create independent set

//     // run TSP on the previous set

//     vector<point_3d> points;

//     // it would be a subset of "sensors"
//     for (auto s : sensors) {
//         auto pos = s.get_position();
//         point_3d new_point = {get<0>(pos), get<1>(pos), 0};
//         points.push_back(new_point);
//     }

//     // Also depot is inserted
//     auto pos_depot = depots[0];
//     point_3d depot = {get<0>(pos_depot), get<1>(pos_depot), 0};
//     points.push_back(depot);

//     TSP tsp(points);
//     tsp.solve();

//     tsp_result_id = tsp.get_path_id();
//     tsp_result = tsp.get_path();
//     cout << "TSP path: ";
//     for (auto p : tsp_result_id) {
//         cout << p << ", ";
//     }
//     cout << endl;

//     double len = tsp.get_length();
//     cout << "TSP len: " << len << endl;
// }

vector<point> algorithms::get_intersection_points(point pa, point pb) {
    // given two points a and b, returns the intersection points
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
    if (distance > 2 * dep->get_sensor_radius() ||
        distance < fabs(dep->get_sensor_radius() - dep->get_sensor_radius())) {
        // cout << "No intersection points found." << endl;
        return int_points;
    }

    // Calculate the distance from the center of circle A to the line joining the intersection points
    double a = (pow(dep->get_sensor_radius(), 2) - pow(dep->get_sensor_radius(), 2) + pow(distance, 2)) / (2 * distance);

    // Calculate the coordinates of the intersection points
    double x3 = x1 + a * (x2 - x1) / distance;
    double y3 = y1 + a * (y2 - y1) / distance;

    double h = sqrt(pow(dep->get_sensor_radius(), 2) - pow(a, 2));

    // Calculate the coordinates of the two intersection points
    double int_x1 = x3 + h * (y2 - y1) / distance;
    double int_y1 = y3 - h * (x2 - x1) / distance;

    double int_x2 = x3 - h * (y2 - y1) / distance;
    double int_y2 = y3 + h * (x2 - x1) / distance;

    // Output the intersection points
    // cout << "Intersection Point 1: (" << int_x1 << ", " << int_y1 << ")" << endl;
    // cout << "Intersection Point 2: (" << int_x2 << ", " << int_y2 << ")" << endl;

    int_ab_1 = {int_x1, int_y1};
    int_ab_2 = {int_x2, int_y2};

    int_points.push_back(int_ab_1);
    int_points.push_back(int_ab_2);

    return int_points;
}

void algorithms::draw_result(vector<vector<tuple<point, int>>> tspn_tours) {
    ofstream htmlFile("output/sensor_deployment.html");

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

    // Draw sensors and depots as before
    for (size_t i = 0; i < dep->get_sensors().size(); ++i) {
        auto s = dep->get_sensors()[i];
        auto pos = s.get_position();

        // for (size_t i = 0; i < sorted_sensors.size(); ++i) {
        //     auto s = sorted_sensors[i];
        //     auto pos = make_tuple(get<0>(s), get<1>(s));

        // Draw sensor circle
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.arc(" << get<0>(pos) << ", " << dep->get_area_width() - get<1>(pos) << ", "
                 << dep->get_sensor_radius()
                 << ", 0, 2 * Math.PI);\n";
        htmlFile << "ctx.fillStyle = 'rgba(0, 0, 255, 0.25)';\n";
        htmlFile << "ctx.fill();\n";
        htmlFile << "ctx.stroke();\n";

//        // Draw sensor shape based on radius_doi
//        htmlFile << "ctx.beginPath();\n";
//        for (int angle = 0; angle < 360; ++angle) {
//            double radius = s.get_radius_doi(angle);
//            double x = get<0>(pos) + radius * cos(angle * M_PI / 180.0);
//            double y = dep->get_area_width() - get<1>(pos) + radius * sin(angle * M_PI / 180.0);
//
//            if (angle == 0) {
//                htmlFile << "ctx.moveTo(" << x << ", " << y << ");\n";
//            } else {
//                htmlFile << "ctx.lineTo(" << x << ", " << y << ");\n";
//            }
//        }
//        htmlFile << "ctx.closePath();\n";
//        htmlFile << "ctx.fillStyle = 'rgba(0, 0, 255, 0.25)';\n";
//        htmlFile << "ctx.fill();\n";
//        htmlFile << "ctx.stroke();\n";

        // Draw sensor center dot
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.arc(" << get<0>(pos) << ", " << dep->get_area_width() - get<1>(pos)
                 << ", 2, 0, 2 * Math.PI);\n";
        htmlFile << "ctx.fillStyle = 'black';\n";
        htmlFile << "ctx.fill();\n";
        htmlFile << "ctx.stroke();\n";

        // Draw sensor label with coordinates (bold)
        htmlFile << "ctx.fillStyle = 'black';\n";
        htmlFile << "ctx.font = 'bold 15px Arial';\n";
        htmlFile << "ctx.fillText('" << i << " (" << get<0>(pos) << ", " << get<1>(pos) << ")', "
                 << get<0>(pos) + 10 << ", " << dep->get_area_width() - get<1>(pos) + 10 << ");\n";
    }

    for (size_t i = 0; i < dep->get_depots().size(); ++i) {
        auto d = dep->get_depots()[i];
        auto pos = d;

        // Draw depot square
        htmlFile << "ctx.fillStyle = 'brown';\n";
        htmlFile << "ctx.fillRect(" << get<0>(pos) - 7.5 << ", " << dep->get_area_width() - get<1>(pos) - 7.5
                 << ", 15, 15);\n";

        // Draw depot label with coordinates
        htmlFile << "ctx.fillStyle = 'black';\n";
        htmlFile << "ctx.font = '15px Arial';\n";
        htmlFile << "ctx.fillText('D" << i << " (" << get<0>(pos) << ", " << get<1>(pos) << ")', "
                 << get<0>(pos) + 10 << ", " << dep->get_area_width() - get<1>(pos) + 10 << ");\n";
    }

    // // Draw TSP circuit connecting points in the order specified by tsp_result
    // htmlFile << "ctx.strokeStyle = 'red';\n";
    // htmlFile << "ctx.lineWidth = 2;\n";
    // htmlFile << "ctx.beginPath();\n";
    // htmlFile << "ctx.moveTo(" << get<0>(tspn_result[0]) << ", " << dep->get_area_width() - get<1>(tspn_result[0]) << ");\n";
    // for (auto &i: tspn_result) {
    //     htmlFile << "ctx.lineTo(" << get<0>(i) << ", " << dep->get_area_width() - get<1>(i) << ");\n";
    // }
    // htmlFile << "ctx.lineTo(" << get<0>(tspn_result[0]) << ", " << dep->get_area_width() - get<1>(tspn_result[0]) << ");\n";
    // htmlFile << "ctx.stroke();\n";

//////////////////////////////////
    // Draw TSP circuit connecting points in the order specified by tspn_tours
    htmlFile << "ctx.strokeStyle = 'red';\n";
    htmlFile << "ctx.lineWidth = 2;\n";
    htmlFile << "ctx.beginPath();\n";
    htmlFile << "ctx.moveTo(" << get<0>(get<0>(tspn_tours[0][0])) << ", " << dep->get_area_width() - get<1>(get<0>(tspn_tours[0][0])) << ");\n";
    for (auto & tspn_tour : tspn_tours) {
        for (auto j: tspn_tour) {
            htmlFile << "ctx.lineTo(" << get<0>(get<0>(j)) << ", " << dep->get_area_width() - get<1>(get<0>(j)) << ");\n";
        }

    }
    htmlFile << "ctx.lineTo(" << get<0>(get<0>(tspn_tours[0][0])) << ", " << dep->get_area_width() - get<1>(get<0>(tspn_tours[0][0])) << ");\n";

    htmlFile << "ctx.stroke();\n";
////////////////////////////////////

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




