#ifndef TOSN_SENSOR_H
#define TOSN_SENSOR_H

#include <iostream>
#include <utility>
#include <vector>

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

    pair<double, double> get_position();
    double get_data_size() const;
    vector<double> get_radius_doi() const;
};


#endif //TOSN_SENSOR_H
