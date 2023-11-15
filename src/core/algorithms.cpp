#include "algorithms.h"

#include <utility>

algorithms::algorithms(deployment *dep) {
    area_length = dep->get_area_length();
    area_width = dep->get_area_width();
    sensor_radius = dep->get_sensor_radius();

    sensors = dep->get_sensors();
    depots = dep->get_depots();
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

bool algorithms::is_within_radius(const sensor& s, point p) {
    double dist = get_distance(s, p);
    return (dist <= sensor_radius);
}

bool algorithms::is_within_radius_doi(const sensor& s, point p) {
    double dist = get_distance(s, p);

    int angle = get_angle(s, p);
    double actual_radius = s.get_radius_doi(angle);

    return (dist <= actual_radius);
}

void algorithms::tsp_neighbors() {
    // sort sensors from left to right

    point p1 = {0, 0};
    point p2 = {10, 10};
    vector<point> aaa = get_intersection_points(p1, p2);
    for (point p : aaa) {
        cout << get<0>(p) << ", " << get<1>(p) << endl;
    }
    cout << endl;

    // create independent set

    // run TSP on the previous set

    vector<point_3d> points;

    // it would be a subset of "sensors"
    for (auto s : sensors) {
        auto pos = s.get_position();
        point_3d newPoint = {pos.first, pos.second, 0};
        points.push_back(newPoint);
    }

    TSP tsp(points);
    tsp.solve();

    tsp_result = tsp.get_path();
    cout << "TSP path: ";
    for (auto p : tsp_result) {
        cout << p << ", ";
    }
    cout << endl;

    double len = tsp.get_length();
    cout << "TSP len: " << len << endl;
}

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
    if (distance > 2 * sensor_radius || distance < fabs(sensor_radius - sensor_radius)) {
//        cout << "No intersection points found." << endl;
        return int_points;
    }

    // Calculate the distance from the center of circle A to the line joining the intersection points
    double a = (pow(sensor_radius, 2) - pow(sensor_radius, 2) + pow(distance, 2)) / (2 * distance);

    // Calculate the coordinates of the intersection points
    double x3 = x1 + a * (x2 - x1) / distance;
    double y3 = y1 + a * (y2 - y1) / distance;

    double h = sqrt(pow(sensor_radius, 2) - pow(a, 2));

    // Calculate the coordinates of the two intersection points
    double int_x1 = x3 + h * (y2 - y1) / distance;
    double int_y1 = y3 - h * (x2 - x1) / distance;

    double int_x2 = x3 - h * (y2 - y1) / distance;
    double int_y2 = y3 + h * (x2 - x1) / distance;

    // Output the intersection points
//    cout << "Intersection Point 1: (" << int_x1 << ", " << int_y1 << ")" << endl;
//    cout << "Intersection Point 2: (" << int_x2 << ", " << int_y2 << ")" << endl;

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

    // Draw sensors and depots as before
    for (size_t i = 0; i < sensors.size(); ++i) {
        auto s = sensors[i];
        auto pos = s.get_position();

        // Draw sensor circle
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.arc(" << get<0>(pos) << ", " << get<1>(pos) << ", " << sensor_radius << ", 0, 2 * Math.PI);\n";
        htmlFile << "ctx.fillStyle = 'rgba(0, 0, 255, 0.25)';\n";
        htmlFile << "ctx.fill();\n";
        htmlFile << "ctx.stroke();\n";

        // Draw sensor center dot
        htmlFile << "ctx.beginPath();\n";
        htmlFile << "ctx.arc(" << get<0>(pos) << ", " << get<1>(pos) << ", 2, 0, 2 * Math.PI);\n";
        htmlFile << "ctx.fillStyle = 'black';\n";
        htmlFile << "ctx.fill();\n";
        htmlFile << "ctx.stroke();\n";

        // Draw sensor label
        htmlFile << "ctx.fillStyle = 'black';\n";
        htmlFile << "ctx.font = '15px Arial';\n";
        htmlFile << "ctx.fillText('" << i << "', " << get<0>(pos) + 10 << ", " << get<1>(pos) + 10 << ");\n";
    }

    for (size_t i = 0; i < depots.size(); ++i) {
        auto d = depots[i];
        auto pos = d;

        // Draw depot square
        htmlFile << "ctx.fillStyle = 'brown';\n";
        htmlFile << "ctx.fillRect(" << get<0>(pos) << ", " << get<1>(pos) << ", 15, 15);\n";

        // Draw depot label
        htmlFile << "ctx.fillStyle = 'black';\n";
        htmlFile << "ctx.font = '15px Arial';\n";
        htmlFile << "ctx.fillText('D" << i << "', " << get<0>(pos) + 10 << ", " << get<1>(pos) + 10 << ");\n";
    }

    // Draw TSP circuit connecting sensors in the order specified by tsp_result
    htmlFile << "ctx.strokeStyle = 'red';\n";
    htmlFile << "ctx.lineWidth = 2;\n";
    htmlFile << "ctx.beginPath();\n";
    htmlFile << "ctx.moveTo(" << get<0>(sensors[tsp_result[0]].get_position()) << ", " << get<1>(sensors[tsp_result[0]].get_position()) << ");\n";
    for (size_t i = 1; i < tsp_result.size(); ++i) {
        auto currentSensor = sensors[tsp_result[i]].get_position();
        htmlFile << "ctx.lineTo(" << get<0>(currentSensor) << ", " << get<1>(currentSensor) << ");\n";
    }
    htmlFile << "ctx.lineTo(" << get<0>(sensors[tsp_result[0]].get_position()) << ", " << get<1>(sensors[tsp_result[0]].get_position()) << ");\n";
    htmlFile << "ctx.stroke();\n";

    htmlFile << "}\n";
    htmlFile << "</script>\n";

    htmlFile << "<canvas id='sensorCanvas' width='" << area_length << "' height='" << area_width << "' style='border:1px solid #000;'></canvas>\n";

    htmlFile << "<script>\n";
    htmlFile << "drawSensorsAndDepots();\n";
    htmlFile << "</script>\n";

    htmlFile << "</body>\n</html>";

    htmlFile.close();
}



