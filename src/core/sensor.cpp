#include "sensor.h"

#include <utility>

sensor::sensor(double x, double y, double data, vector<double> r_doi) {
    pos_x = x;
    pos_y = y;
    data_size = data;
    radius_doi = std::move(r_doi);
}

pair<double, double> sensor::get_position() {
    return make_pair(pos_x, pos_y);
}

double sensor::get_data_size() const {
    return data_size;
}

vector<double> sensor::get_radius_doi() const {
    return radius_doi;
}

ostream &operator<<(ostream &os, const sensor &s) {
    cout << "pos=(" << s.pos_x << "," << s.pos_y << ")" << endl;
    cout << "data_size=" << s.data_size << endl;

    return os;
}
