#include "algorithms.h"

algorithms::algorithms(deployment *dep) {
    num_sensors = dep->get_num_sensors();
    sensor_radius = dep->get_sensor_radius();

    sensors = dep->get_sensors();
    depots = dep->get_depots();
}

double algorithms::get_distance(sensor s1, sensor s2) {
    pair<double, double> pos_s1 = s1.get_position();
    pair<double, double> pos_s2 = s2.get_position();
    return get_distance(pos_s1, pos_s2);
}

double algorithms::get_distance(depot d, sensor s) {
    pair<double, double> pos_s = s.get_position();
    return get_distance(d, pos_s);
}

double algorithms::get_distance(point p1, point p2) {
    double x1 = get<0>(p1);
    double y1 = get<1>(p1);
    double x2 = get<0>(p2);
    double y2 = get<1>(p2);
    return sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
}

void algorithms::algorithm_1() {
    cout << "alg1" << endl;
    cout << "distance s0-s1=" << get_distance(sensors[0], sensors[1]) << endl;
}

void algorithms::algorithm_2() {
    cout << "alg2" << endl;
    cout << "distance depot0-s5=" << get_distance(depots[0], sensors[5]) << endl;
}
