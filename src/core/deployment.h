#ifndef TOSN_DEPLOYMENT_H
#define TOSN_DEPLOYMENT_H

#include <iostream>
#include <vector>

#include "sensor.h"

using namespace std;

class deployment {
private:
    int num_sensors;
    int area_length;
    int area_width;
    double sensor_radius;
    double doi;

    vector<sensor> sensors;

public:
    deployment(int, int, int, double, double);
};


#endif //TOSN_DEPLOYMENT_H
