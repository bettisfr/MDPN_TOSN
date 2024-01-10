#ifndef TOSN_DEPLOYMENT_H
#define TOSN_DEPLOYMENT_H

#include <iostream>
#include <vector>
#include <random>
#include <numeric>

#include "../input.h"
#include "../definitions.h"
#include "sensor.h"

using namespace std;

class deployment {

    friend ostream &operator<<(ostream &os, const deployment &s);

private:
    int num_sensors;
    int num_depots;
    int area_length;
    int area_width;
    int max_data;
    double sensor_radius;
    double sensor_radius_doi_percentage;
    double doi;
    double epsilon;

    int energy_budget;
    double energy_cons_fly;
    double energy_cons_hover;

    int wireless_technology;

    // Speed of light in m/s
    const double c = 3e8;

    // Frequency in Hz
    double f_c;

    // Transmitting power in dBm
    double P_Tx;

    // Noise power in watts
    double N;

    // Channel bandwidth
    double B;

    // Weibull parameters
    const double shape = 0.67;
    const double scale = 0.16;

    vector<sensor> sensors;
    vector<point> depots;

public:
    explicit deployment(const input &);

    // Free Space Path Loss
    double get_FSPL() const;

    // Data Transfer Rate
    double get_DTR(double) const;

    int get_num_sensors() const;

    int get_area_length() const;

    int get_area_width() const;

    double get_sensor_radius() const;

    double get_sensor_radius_doi() const;

    double get_epsilon() const;

    double get_doi() const;

    vector<sensor> get_sensors();

    vector<point> get_depots();

    double get_energy_budget() const;

    double get_energy_cons_fly() const;

    double get_energy_cons_hover() const;

    bool check_feasibility() const;
};


#endif //TOSN_DEPLOYMENT_H
