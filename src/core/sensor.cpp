#include "sensor.h"

sensor::sensor(double _x, double _y, double _data, double _r) {
    x = _x;
    y = _y;
    data = _data;
    r = _r;

    for (double & i : r_doi) {
        i = r;
    }
}

pair<double, double> sensor::get_position() {
    return make_pair(x, y);
}

double sensor::get_data() const {
    return data;
}

double sensor::get_radius() const {
    return r;
}

double *sensor::get_radius_doi() {
    return r_doi;
}

ostream &operator<<(ostream &os, const sensor &s) {
    cout << "pos=(" << s.x << "," << s.y << ")" << endl;
    cout << "data=" << s.data << endl;
    cout << "radius=" << s.r << endl;

    return os;
}
