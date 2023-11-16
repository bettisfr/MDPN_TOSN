#include "sensor.h"

#include <utility>

sensor::sensor(double x, double y, double data, vector<double> r_doi) {
    pos_x = x;
    pos_y = y;
    data_size = data;
    radius_doi = std::move(r_doi);
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
