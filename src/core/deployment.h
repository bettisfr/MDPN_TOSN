#ifndef TOSN_DEPLOYMENT_H
#define TOSN_DEPLOYMENT_H

#include <iostream>
#include <vector>
#include <random>
#include <numeric>

#include "../input.h"
#include "sensor.h"

using namespace std;

typedef tuple<double, double> point;

class deployment {

    friend ostream &operator<<(ostream &os, const deployment &s);

private:
    int num_sensors;
    int num_depots;
    int area_length;
    int area_width;
    int max_data;
    double sensor_radius;
    double doi;
    double epsilon;

    double data_transfer_rate;
    int energy_budget;
    double energy_cons_fly;
    double energy_cons_hover;

    const double f_c = 2.4e9;    // Frequency in Hz (for 2.4 GHz Wi-Fi)
    const double c = 3e8;        // Speed of light in m/s
    const double P_Tx = 20;      // Transmitting power in dBm
    const double N = 1e-8;       // Noise power in watts
    const double B = 20e6;       // 20 MHz

    // Weibull parameters
    const double shape = 0.67;
    const double scale = 0.16;

    vector<sensor> sensors;
    vector<point> depots;

public:
    explicit deployment(const input &);

    // Free Space Path Loss
    double get_FSPL();

    // Data Transfer Rate
    double get_DTR(double);

    int get_num_sensors() const;

    int get_area_length() const;

    int get_area_width() const;

    double get_sensor_radius() const;
    
    double get_epsilon() const;

    vector<sensor> get_sensors();

    vector<point> get_depots();

    int get_energy_budget() const;

    int get_energy_cons_fly() const;

    int get_energy_cons_hover() const;

    double get_data_transfer_rate() const;
};


#endif //TOSN_DEPLOYMENT_H
