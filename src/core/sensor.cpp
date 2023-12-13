#include "sensor.h"

#include <utility>

sensor::sensor(double x, double y, double data, vector<double> r_doi) {
    id = -1;
    pos_x = x;
    pos_y = y;
    data_size = data;
    radius_doi = std::move(r_doi);
}

sensor::sensor(double x, double y) {
    pos_x = x;
    pos_y = y;
    data_size = -1;
}


tuple<double, double> sensor::get_position() const {
    return make_tuple(pos_x, pos_y);
}

double sensor::get_data_size() const {
    return data_size;
}

double sensor::get_radius_doi(int angle) const {
    return radius_doi[angle];
}

ostream &operator<<(ostream &os, const sensor &s) {
    cout << "pos=(" << s.pos_x << "," << s.pos_y << "); ";
    cout << "data_size=" << s.data_size;
    for (int i = 0; i < 360; i++) {
        cout << s.radius_doi[i] << ", ";
    }

    return os;
}

double sensor::get_pos_x() const {
    return pos_x;
}

double sensor::get_pos_y() const {
    return pos_y;
}

bool sensor::operator==(const sensor &other) const {
    return (pos_x == other.pos_x) && (pos_y == other.pos_y);
}

void sensor::set_id(int new_id) {
    id = new_id;
}

int sensor::get_id() const {
    return id;
}



