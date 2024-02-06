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
    int id;
    double pos_x;
    double pos_y;
    double data_size;
    vector<double> radius_doi;

public:
    sensor(double, double, double, vector<double>);

    sensor(double, double);

    bool operator==(const sensor &other) const;

    tuple<double, double> get_position() const;

    double get_pos_x() const;

    double get_pos_y() const;

    void set_pos(double, double);

    double get_data_size() const;

    double get_radius_doi(int) const;

    int get_id() const;

    void set_id(int);
};

#endif //TOSN_SENSOR_H
