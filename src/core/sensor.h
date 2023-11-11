#ifndef TOSN_SENSOR_H
#define TOSN_SENSOR_H

#include <iostream>
#include <utility>

using namespace std;

class sensor {

    friend ostream &operator<<(ostream &os, const sensor &s);

private:
    double x;
    double y;
    double data;
    double r;
    double r_doi[360];

public:
    sensor(double, double, double, double);

    pair<double, double> get_position();
    double get_data() const;
    double get_radius() const;
    double* get_radius_doi();
};


#endif //TOSN_SENSOR_H
