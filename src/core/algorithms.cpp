#include "algorithms.h"

#include <utility>
#include <unordered_set>

algorithms::algorithms(deployment *m_dep) {
    dep = m_dep;
}

void algorithms::ApproxTSPN_S() {
    cout << "ApproxTSPN_S: " << endl;
    tsp_neighbors(dep->get_sensors());
    point depot;
    depot = dep->get_depots()[0];
    tsp_split(dep->get_energy_budget(), depot);
    cout << "Number of drones: " << tspn_tours.size() << endl;
    draw_result();
}

void algorithms::ApproxMPN_S() {
    cout << "ApproxMPN_S: " << endl;
    point depot;
    depot = dep->get_depots()[0];
    ApproxMPN(depot);
    cout << "Number of drones: " << tspn_tours.size() << endl;
    draw_result();
}

void algorithms::ApproxTSPN_M() {
    cout << "ApproxTSPN_M: " << endl;
    int ecf = dep->get_energy_cons_fly();
    double radius = dep-> get_sensor_radius();
    double R_0_f = radius * ecf; 
    vector<point> depots = dep->get_depots();
    point depot;
    double A;
    vector<double> A_d(depots.size(), 0.0);
    // find the depot
    for (int i = 0; i < depots.size(); i++){
        for (const auto &sensor: dep->get_sensors()){
            double dist = get_distance(sensor, depots[i]);
            auto pos = sensor.get_position();
            double energy_hovering = compute_energy_hovering(make_tuple(get<0>(pos), get<1>(pos), sensor.get_data_size(), 0));
            if (dist <= radius){
                if (energy_hovering > A_d[i]){
                    A_d[i] = energy_hovering;
                }
            } else {
                double a = 2 * (dist * ecf + (energy_hovering / 2) + R_0_f);
                if (a > A_d[i]){
                    A_d[i] = a;
                }
            }
        } 
    }
    int indesx = distance(A_d.begin(), min_element(A_d.begin(), A_d.end()));
    depot = depots[indesx];

    tsp_neighbors(dep->get_sensors());
    tsp_split(dep->get_energy_budget(), depot);
    cout << "Number of drones: " << tspn_tours.size() << endl;
    draw_result();
}

void algorithms::ApproxMPN_M() {
    cout << "ApproxMPN_M: " << endl;
    vector<point> depots = dep->get_depots();
    vector<sensor> sensors = dep->get_sensors();
    vector<vector<tuple<double, double, int>>> tours;
    vector<tuple<double, double, int>> vec;
    vec.push_back(make_tuple(0.0, 0.0, 0));

    for (int i = 0; i < 2 * sensors.size(); i++){
        tours.push_back(vec);
    }
    
    for (auto depot : depots){
        ApproxMPN(depot);
        if (tspn_tours.size() < tours.size()){
            tours = tspn_tours;
        }
    }

    tspn_tours = tours;
    cout << "Number of drones: " << tspn_tours.size() << endl;
    draw_result();
    
}


double algorithms::get_distance(sensor s, point p) {
    point pos_s = s.get_position();

    double x1, y1, x2, y2;
    tie(x1, y1) = p;
    tie(x2, y2) = pos_s;

    double delta_x = x2 - x1;
    double delta_y = y2 - y1;

    return sqrt(delta_x * delta_x + delta_y * delta_y);
}

int algorithms::get_angle(sensor s, point p) {
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

double algorithms::compute_energy_hovering(tuple<double, double, double, double>sensor) {
    double dtr = dep->get_data_transfer_rate();
    int ech = dep->get_energy_cons_hover(); 
    double sensor_data = get<2>(sensor);
    double required_time = sensor_data / dtr;
    return required_time * ech;
}

double algorithms::tour_cost(vector<tuple<double, double, int>> T, int start, int end, point depot) {
    int ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    int ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)

    double cost_T_k = 0;

    if (get<0>(depot) != -1){
        tuple<double, double, double, double> sensor1;
        tuple<double, double, double, double> sensor2;

        double dist1 = sqrt(pow(get<0>(T[0]) - get<0>(T[1]), 2) + pow(get<1>(T[0]) - get<1>(T[1]), 2));
        double energy1_flying = dist1 * ecf;

        sensor1 = sorted_sensors[get<2>(T[1])];
        double energy1_hovering = compute_energy_hovering(sensor1);

        double dist2 = sqrt(pow(get<0>(T[0]) - get<0>(T[T.size()-2]), 2) + pow(get<1>(T[0]) - get<1>(T[T.size()-2]), 2));
        double energy2_flying = dist2 * ecf;

        sensor2 = sorted_sensors[get<2>(T[T.size()-2])];
        double energy2_hovering = compute_energy_hovering(sensor2);

        for (int p = start; p < end; p++){
            cost_T_k = cost_T_k + tspn_cost[p];
        }
        cost_T_k = cost_T_k + energy1_flying + energy2_flying + energy1_hovering / 2 + energy2_hovering / 2;

    } else {
        for (int p = start; p < end; p++){
            cost_T_k += tspn_cost[p]; 
        }
    }

    return cost_T_k;
}

void algorithms::ApproxMPN(point depot) {
    sorted_sensors.clear();
    for (const auto &sensor: dep->get_sensors()) {       // it can be removed
        auto pos = sensor.get_position();
        sorted_sensors.emplace_back(get<0>(pos), get<1>(pos), sensor.get_data_size(), 0);
    }

    // get eps and compute t
    double epsilon = dep->get_epsilon();
    int t = ceil(log2(1/epsilon));


    int energy_budget = dep->get_energy_budget();
    energy_budget = energy_budget * 1000;

    int ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    int ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)
    double radius = dep-> get_sensor_radius();
    double R_0_f = radius * ecf;                // R_0

    vector<vector<tuple<double, double, double, double>>> V;
    for (int i = 0; i <= t; i++){
       V.push_back(vector<tuple<double, double, double, double>>());
    }    

    // for each s in sorted_sensors compute lambd_v,d and j
    for (auto s: sorted_sensors){
        double dist = sqrt(pow(get<0>(s) - get<0>(depot), 2) + pow(get<1>(s) - get<1>(depot), 2));
        double energy_flying = dist * ecf;

        double energy_hovering = compute_energy_hovering(s);

        double lamdba = energy_flying + energy_hovering / 2 + R_0_f;

        int j = floor(log2((energy_budget - 2 * lamdba) / (epsilon * energy_budget)) + 1);

        // find V_0,..., V_t
        if (double((1-epsilon) * (energy_budget / 2)) < lamdba && lamdba <= double(energy_budget / 2)){
            V[0].push_back(s);
        } else {
            V[j].push_back(s);
        }
    }

    vector<vector<tuple<double, double, int>>> tours;
    for (int i = 0; i < V.size(); i++){
        // for each V[i], run Minimum UAV Deployment Problem with Neighborhoods with budget 2^{j-1} epsilon B
        vector<vector<tuple<double, double, int>>> C;
        C = approAlgNei(V[i], i);

        for (int j = 0 ; j < C.size(); j++){
            C[j].insert(C[j].begin(), make_tuple(get<0>(depot), get<1>(depot), -1));
            C[j].push_back(make_tuple(get<0>(depot), get<1>(depot), -1));
        }
        for (int j = 0 ; j < C.size(); j++){
            tours.push_back(C[j]);
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
    
}


vector<vector<tuple<double, double, int>>> algorithms::approAlgNei(vector<tuple<double, double, double, double>> V, int j) {    
    int total_budget = dep->get_energy_budget();
    double epsilon = dep->get_epsilon();
    double energy_budget = pow(2, j-1) * epsilon * total_budget;
    energy_budget = energy_budget * 1000;

    int ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    int ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)
    double radius = dep-> get_sensor_radius();

    vector<vector<tuple<double, double, int>>> C;
    for (auto s: V){
        vector<tuple<double, double, int>> A;
        A.emplace_back(make_tuple(get<0>(s), get<1>(s), -1));
        C.emplace_back(A);
    }
    
    sort(V.begin(), V.end());
    int n = V.size();

    // G', n x n, pairwise weight between points
    double **G1;
    for (int i = 2; i <= n; i++){
        double delta = energy_budget / i;
        // Algorithm 4 in [27]
        G1 = new double *[n];
        for (int j = 0; j < n; j++) {
            G1[j] = new double[n];
            tuple<double, double, double, double> sensor_j;
            sensor_j = V[j];
            double energy_j_hovering = compute_energy_hovering(sensor_j);

            for (int k = 0; k < n; k++) {
                // compute C(j, k)
                tuple<double, double, double, double> sensor_k;
                sensor_k = V[k];
                double dist = sqrt(pow(get<0>(sensor_j) - get<0>(sensor_k), 2) + pow(get<1>(sensor_j) - get<1>(sensor_k), 2));

                double energy_k_hovering = compute_energy_hovering(sensor_k);

                if (dist <= 2 * radius){   // R_0
                    double w1 = energy_j_hovering + energy_k_hovering;
                    // G''
                    if (w1 <= delta){
                        G1[j][k] = w1;
                    } else {
                        G1[j][k] = 0.;
                    }
                } else {
                    double w2 = (dist - 2 * radius) * ecf + energy_j_hovering + energy_k_hovering;  // R_0
                    // G''
                    if (w2 <= delta){
                        G1[j][k] = w2;
                    } else {
                        G1[j][k] = 0.;
                    }
                }
            }
        }

        // find connected components
        vector<vector<int>> adjacencyList;
        for (int i = 0; i < n; i++){
            adjacencyList.push_back(vector<int>());
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (G1[i][j] > 0){
                    adjacencyList[i].push_back(j);
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
                components.push_back(vector<tuple<double, double, double, double>>());
                for (int u : connectedComponent) {
                    components[p].push_back(V[u]);
                }
                p++;
            }
        }

        vector<vector<tuple<double, double, int>>> C1;
        // for each component run tspn
        for (int j = 0; j < components.size(); j++){
            vector<sensor> comp_sensors;
            for (auto s : components[j]){
                comp_sensors.emplace_back(get<0>(s), get<1>(s), get<2>(s), vector<double>());
            }
            tsp_neighbors(comp_sensors);
            tsp_split(energy_budget / 1000, make_tuple(-1, -1));
            for (int i = 0; i < tspn_tours.size(); i++){
                C1.push_back(tspn_tours[i]);                
            }            
        }
        if (C1.size() < C.size()){
            C = C1;
        }
    }    
    return C;
}

void algorithms::DFS(int v, unordered_set<int>& visited, unordered_set<int>& connectedComponent, vector<vector<int>> adjacencyList) {
    visited.insert(v);
    connectedComponent.insert(v);
    for (int neighbor : adjacencyList[v]) {
        if (visited.find(neighbor) == visited.end()) {
            DFS(neighbor, visited, connectedComponent, adjacencyList);
        }
    }
}


void algorithms::tsp_neighbors(vector<sensor> sensors) {
    tspn_result.clear();
    tspn_cost.clear();
    sorted_sensors.clear();
    tspn_tours.clear();
    
    for (const auto &sensor: sensors) {      // dep->get_sensors()     
        auto pos = sensor.get_position();
        sorted_sensors.emplace_back(get<0>(pos), get<1>(pos), sensor.get_data_size(), 0);
    }

    // sort sensors based on x-coordinates
    sort(sorted_sensors.begin(), sorted_sensors.end());
    // cout << "sensors : " << endl;
    // for (int i = 0; i < sorted_sensors.size(); i++) {
    //     cout << get<0>(sorted_sensors[i]) << ", " << get<1>(sorted_sensors[i]) << endl;
    // }
    // cout << "**********" << endl;

    map<int, map<int, vector<point>>> intersections;
    vector<int> checked_sensors(sorted_sensors.size());

    for (int i = 0; i < sorted_sensors.size(); i++) {
        if (checked_sensors[i] == 0) {
            point p1 = {get<0>(sorted_sensors[i]), get<1>(sorted_sensors[i])};
            int intersec = 0;
            for (int j = i + 1; j < sorted_sensors.size(); j++) {
                if (checked_sensors[j] == 0) {
                    // add another if to check whether the distance between i and j is less than 2R
                    // find intersec of i and j
                    point p2 = {get<0>(sorted_sensors[j]), get<1>(sorted_sensors[j])};
                    vector<point> intersec_points = get_intersection_points(p1, p2);
                    if (!intersec_points.empty()) {
                        intersections[i][j] = intersec_points;
                        intersec = 1;
                        checked_sensors[j] = 1;
                    }
                }
            }
            checked_sensors[i] = 1;
            if (intersec == 0) {
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
    vector<point_3d> points;
    //vector<int> I;
    for (const auto &p: intersections) {
        //I.push_back(p.first);
        point_3d new_point = {get<0>(sorted_sensors[p.first]), get<1>(sorted_sensors[p.first]), get<2>(sorted_sensors[p.first])};
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

    int ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    int ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)

    //vector<double> tspn_cost;
    int counter = 0;
    for (int i = 0; i < tsp_result_id.size(); i++) {
        // use points to access sorted_node and then intersections
        tuple<double, double, double, double> sensor1;
        tuple<double, double, double, double> sensor2;

        sensor1 = make_tuple(points[tsp_result_id[i]].x, points[tsp_result_id[i]].y, points[tsp_result_id[i]].z, 0);
        auto it = find_if(sorted_sensors.begin(), sorted_sensors.end(),
                          [&sensor1](const auto &element) { return element == sensor1; });

        int index1 = distance(sorted_sensors.begin(), it);

        counter++;
        if (counter == tsp_result_id.size()) {
            tspn_result.push_back(make_tuple(get<0>(sorted_sensors[index1]), get<1>(sorted_sensors[index1]), index1));
            break;
        }

        tspn_result.push_back(make_tuple(points[tsp_result_id[i]].x, points[tsp_result_id[i]].y, index1));
        
        sensor2 = make_tuple(points[tsp_result_id[i + 1]].x, points[tsp_result_id[i + 1]].y, points[tsp_result_id[i + 1]].z, 0);

        int sen1_to_sen2 = 0;
        int len_sen1 = intersections[index1].size();
        for (auto j_intersecs: intersections[index1]) {
            sen1_to_sen2++;

            double energy1_hovering = compute_energy_hovering(sensor1);
            double energy2_hovering = compute_energy_hovering(sensor2);
            double energy_j_hovering = compute_energy_hovering(sorted_sensors[j_intersecs.first]);
            double energy_hovering = energy1_hovering / 2 + energy_j_hovering + energy2_hovering / 2;

            int skip = 0;
            vector<tuple<double, point>> dist_points;
            for (point p: j_intersecs.second) {  // 2 points
                if (get<0>(p) == -1) {
                    skip = 1;
                    break;
                } else {
                    double dist_sen1_p = sqrt(
                            pow(get<0>(p) - get<0>(sensor1), 2) + pow(get<1>(p) - get<1>(sensor1), 2));
                    double dist_p_sen2 = sqrt(
                            pow(get<0>(sensor2) - get<0>(p), 2) + pow(get<1>(sensor2) - get<1>(p), 2));
                    
                    double distance = dist_sen1_p + dist_p_sen2;
                    double energy_flying = distance * ecf;
                    double total_energy = energy_hovering + energy_flying;

                    dist_points.push_back(make_tuple(total_energy, p));
                }
            }
            // add center of j_intersecs.first
            double x_j = get<0>(sorted_sensors[j_intersecs.first]);
            double y_j = get<1>(sorted_sensors[j_intersecs.first]);

            double dist_sen1_p = sqrt(pow(x_j - get<0>(sensor1), 2) + pow(y_j - get<1>(sensor1), 2));
            double dist_p_sen2 = sqrt(pow(get<0>(sensor2) - x_j, 2) + pow(get<1>(sensor2) - y_j, 2));
            double distance = dist_sen1_p + dist_p_sen2;
            double energy_flying = distance * ecf;

            double total_energy = energy_hovering + energy_flying;

            dist_points.push_back(make_tuple(total_energy, make_tuple(x_j, y_j)));

            if (skip == 1) {
                double dist_sen1_sen2 = sqrt(pow(get<0>(sensor2) - get<0>(sensor1), 2) + pow(get<1>(sensor2) - get<1>(sensor1), 2));
                double energy_flying = dist_sen1_sen2 * ecf;
                double energy_hovering = energy1_hovering / 2 + energy2_hovering / 2;
                double total_energy = energy_hovering + energy_flying;
                tspn_cost.push_back(total_energy);
                break;
            } else {
                // select the one that has the minimum distance
                sort(dist_points.begin(), dist_points.end());
                point pos = get<1>(dist_points[0]);
                tspn_result.push_back(make_tuple(get<0>(pos), get<1>(pos), j_intersecs.first));

                double dist_sen1_pos = sqrt(pow(get<0>(pos) - get<0>(sensor1), 2) + pow(get<1>(pos) - get<1>(sensor1), 2));
                double energy_flying = dist_sen1_pos * ecf;
                double energy_hovering = (energy1_hovering + energy_j_hovering) / 2;
                double total_energy = energy_hovering + energy_flying;
                tspn_cost.push_back(total_energy);

                if (sen1_to_sen2 == len_sen1){   // or counter == tsp_result_id.size()-1
                    // energy from pos to sen2
                    double dist_pos_sen2 = sqrt(pow(get<0>(pos) - get<0>(sensor2), 2) + pow(get<1>(pos) - get<1>(sensor2), 2));
                    double energy_flying = dist_pos_sen2 * ecf;
                    double energy_hovering = (energy2_hovering + energy_j_hovering) / 2;
                    double total_energy = energy_hovering + energy_flying;
                    tspn_cost.push_back(total_energy);
                }

                sensor1 = make_tuple(get<0>(pos), get<1>(pos), get<2>(sorted_sensors[j_intersecs.first]), 0);

            }
        }
    }   

    // cout << "tsp : " << endl;
    // for (int i=0 ; i < tspn_result.size(); i++) {
    //     cout << "(" << get<0>(tspn_result[i]) << ", " << get<1>(tspn_result[i]) << ")" << " : " << get<2>(tspn_result[i]) << " : " << tspn_cost[i] << endl;
    // }
    // cout << endl;
}


void algorithms::tsp_split(int energy_budget, point depot) {
    tspn_tours.clear();                   
    energy_budget = energy_budget * 1000;

    int i = 0;
    int j = 0;
    int n = sorted_sensors.size();
    
    while (j <= n-1){
        vector<tuple<double, double, int>> T_k; 
        if (get<0>(depot) != -1){
            T_k.push_back(make_tuple(get<0>(depot), get<1>(depot), -1));
        }
        for (int p = i; p <= j; p++){
            T_k.push_back(tspn_result[p]);
        }
        if (get<0>(depot) != -1){
            T_k.push_back(make_tuple(get<0>(depot), get<1>(depot), -1));
        }

        double cost = tour_cost(T_k, i, j, depot);

        if (cost <= energy_budget){
            if (j == n-1){
                tspn_tours.push_back(T_k);
                break;
            } else {
                j++;
            }
        } else {
            j--;
            // remove j from T_k
            if (get<0>(depot) != -1){
                T_k.erase(T_k.begin() + T_k.size()-2);
            } else {
                T_k.erase(T_k.begin() + T_k.size()-1);
            }
            
            tspn_tours.push_back(T_k);
            T_k.clear();

            j++;
            i = j;
        }
    }

    // for (int i = 0; i < tspn_tours.size(); i++){
    //     for (auto j: tspn_tours[i]){
    //         cout<< get<0>(j) << ", " << get<1>(j) << endl;
    //     }
    //     cout << "-------" << endl;
    // }

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
    double a =
            (pow(dep->get_sensor_radius(), 2) - pow(dep->get_sensor_radius(), 2) + pow(distance, 2)) / (2 * distance);

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

void algorithms::draw_result() {
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
        htmlFile << "ctx.arc(" << get<0>(pos) << ", " << dep->get_area_width() - get<1>(pos) << ", " << dep->get_sensor_radius()
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
        htmlFile << "ctx.arc(" << get<0>(pos) << ", " << dep->get_area_width() - get<1>(pos) << ", 2, 0, 2 * Math.PI);\n";
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
        htmlFile << "ctx.fillRect(" << get<0>(pos) - 7.5 << ", " << dep->get_area_width() - get<1>(pos) - 7.5 << ", 15, 15);\n";

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
    htmlFile << "ctx.moveTo(" << get<0>(tspn_tours[0][0]) << ", " << dep->get_area_width() - get<1>(tspn_tours[0][0]) << ");\n";
    for (int i = 0; i < tspn_tours.size(); i++){
        for (auto j: tspn_tours[i]){
            htmlFile << "ctx.lineTo(" << get<0>(j) << ", " << dep->get_area_width() - get<1>(j) << ");\n";
        }
        
    }
    htmlFile << "ctx.lineTo(" << get<0>(tspn_tours[0][0]) << ", " << dep->get_area_width() - get<1>(tspn_tours[0][0]) << ");\n";
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



