#include "algorithms.h"

#include <utility>

algorithms::algorithms(deployment *m_dep) {
    dep = m_dep;
}

void algorithms::algorithm_1() {
    cout << "alg1" << endl;
    tsp_neighbors();
}

void algorithms::algorithm_2() {
    cout << "alg2" << endl;
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


void algorithms::tsp_neighbors() {
    vector<tuple<double, double, double, double>> sorted_sensors;
    for (const auto &sensor: dep->get_sensors()) {
        auto pos = sensor.get_position();
        sorted_sensors.emplace_back(get<0>(pos), get<1>(pos), 0, 0);
    }

    // sort sensors based on x-coordinates
    sort(sorted_sensors.begin(), sorted_sensors.end());
    int s = sorted_sensors.size();
    cout << "sensors : " << endl;
    for (int i = 0; i < s; i++) {
        cout << get<0>(sorted_sensors[i]) << ", " << get<1>(sorted_sensors[i]) << endl;
    }

    map<int, map<int, vector<point>>> intersections;
    vector<int> checked_sensors(sorted_sensors.size());

    for (int i = 0; i < sorted_sensors.size(); i++) {
        if (checked_sensors[i] == 0) {
            point p1 = {get<0>(sorted_sensors[i]), get<1>(sorted_sensors[i])};
            int intersec = 0;
            for (int j = i + 1; j < sorted_sensors.size(); j++) {
                if (checked_sensors[j] == 0) {
                    // add another if to check whether the distance of i and j is less than 2R
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
                intersections[i][-1] = intersec_points;
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
    vector<int> I;
    for (const auto &p: intersections) {
        I.push_back(p.first);
        point_3d new_point = {get<0>(sorted_sensors[p.first]), get<1>(sorted_sensors[p.first]), 0};
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


    // ####################
    int energy_budget = dep->get_energy_budget();

    int ecf = dep->get_energy_cons_fly(); // every meter (in J/m)
    int ech = dep->get_energy_cons_hover(); // every second (in J/s)
    double dtr = dep->get_data_transfer_rate(); // constant DTR (in MB/s)

    sensor first_sensor = dep->get_sensors()[0];
    double sensor_data = first_sensor.get_data_size();
    double required_time = sensor_data / dtr; // total seconds
    double energy_hovering = required_time * ech;

    double test_distance = 100.0; // total meters
    double energy_flying = test_distance * ecf;

    double total_energy = energy_hovering + energy_flying;
//    cout << ((total_energy <= energy_budget) ? ":)" : ":((") << endl;
    // ####################


    // vector<tuple<double, double>> tspn_result;
    int counter = 0;
    for (int i = 0; i < tsp_result_id.size(); i++) {
        // use points to access sorted_node and then intersections
        counter++;
        if (counter == tsp_result_id.size()) {
            tspn_result.push_back(make_tuple(points[tsp_result_id[tsp_result_id.size() - 1]].x,
                                             points[tsp_result_id[tsp_result_id.size() - 1]].y));
            break;
        }

        tspn_result.push_back(make_tuple(points[tsp_result_id[i]].x, points[tsp_result_id[i]].y));

        tuple<double, double, double, double> sensor1;
        tuple<double, double, double, double> sensor2;

        sensor1 = make_tuple(points[tsp_result_id[i]].x, points[tsp_result_id[i]].y, 0, 0);
        sensor2 = make_tuple(points[tsp_result_id[i + 1]].x, points[tsp_result_id[i + 1]].y, 0, 0);

        auto it = find_if(sorted_sensors.begin(), sorted_sensors.end(),
                          [&sensor1](const auto &element) { return element == sensor1; });

        int index1 = distance(sorted_sensors.begin(), it);

        for (auto j_intersecs: intersections[index1]) {
            int skip = 0;
            //vector<tuple<point, double>> points_dist;
            vector<tuple<double, point>> dist_points;
            for (point p: j_intersecs.second) {  // 2 points
                if (get<0>(p) == -1) {
                    skip = 1;
                    break;
                } else {
                    // compute the distance from sensor1 to p and from p to sensor2
                    double dist_sen1_p = sqrt(
                            pow(get<0>(p) - get<0>(sensor1), 2) + pow(get<1>(p) - get<1>(sensor1), 2));
                    double dist_p_sen2 = sqrt(
                            pow(get<0>(sensor2) - get<0>(p), 2) + pow(get<1>(sensor2) - get<1>(p), 2));

                    double distance = dist_sen1_p + dist_p_sen2;
                    dist_points.push_back(make_tuple(distance, p));
                }
            }
            // add center of j_intersecs.first
            double x_j = get<0>(sorted_sensors[j_intersecs.first]);
            double y_j = get<1>(sorted_sensors[j_intersecs.first]);

            double dist_sen1_p = sqrt(pow(x_j - get<0>(sensor1), 2) + pow(y_j - get<1>(sensor1), 2));
            double dist_p_sen2 = sqrt(pow(get<0>(sensor2) - x_j, 2) + pow(get<1>(sensor2) - y_j, 2));

            double distance = dist_sen1_p + dist_p_sen2;
            dist_points.push_back(make_tuple(distance, make_tuple(x_j, y_j)));

            if (skip == 1) {
                break;
            } else {
                // select the one that has the minimum distance
                sort(dist_points.begin(), dist_points.end());
                point pos = get<1>(dist_points[0]);
                tspn_result.push_back(pos);
                sensor1 = sensor1 = make_tuple(get<0>(pos), get<1>(pos), 0, 0);
            }
        }
    }

    cout << "tsp : " << endl;
    for (auto p: tspn_result) {

        cout << "(" << get<0>(p) << ", " << get<1>(p) << ")" << " ; ";
    }

    draw_result();

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

    // Draw TSP circuit connecting points in the order specified by tsp_result
    htmlFile << "ctx.strokeStyle = 'red';\n";
    htmlFile << "ctx.lineWidth = 2;\n";
    htmlFile << "ctx.beginPath();\n";
    htmlFile << "ctx.moveTo(" << get<0>(tspn_result[0]) << ", " << dep->get_area_width() - get<1>(tspn_result[0]) << ");\n";
    for (auto &i: tspn_result) {
        htmlFile << "ctx.lineTo(" << get<0>(i) << ", " << dep->get_area_width() - get<1>(i) << ");\n";
    }
    htmlFile << "ctx.lineTo(" << get<0>(tspn_result[0]) << ", " << dep->get_area_width() - get<1>(tspn_result[0]) << ");\n";
    htmlFile << "ctx.stroke();\n";


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



