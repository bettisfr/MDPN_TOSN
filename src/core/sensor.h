#ifndef TOSN_SENSOR_H
#define TOSN_SENSOR_H

#include <iostream>
#include <utility>
#include <vector>
#include <tuple>

using namespace std;

class sensor {

    friend ostream &operator<<(ostream &os, const sensor &s);

private:
    double pos_x;
    double pos_y;
    double data_size;
    vector<double> radius_doi;

public:
    sensor(double, double, double, vector<double>);

    tuple<double, double> get_position() const;

    double get_pos_x() const;

    double get_pos_y() const;

    double get_data_size() const;

    double get_radius_doi(int) const;
};

#endif //TOSN_SENSOR_H
