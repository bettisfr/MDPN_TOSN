#ifndef TOSN_ENERGY_H
#define TOSN_ENERGY_H

#include <complex>
#include <iostream>

using namespace std;

class energy {
public:
    energy();

    double get_energy_hovering(double, double, double);

    double get_energy_movement(double, double, double, double, double);

    void evaluate();

private:
    double m_drone;
    double m_battery;

    int num_rotors;
    double diameter;

    double pressure;
    double R;
    double temperature;
    double rho;

    double g;

    double eta;

    double drag_coefficient_drone;
    double drag_coefficient_battery;
    double drag_coefficient_package;

    double projected_area_drone;
    double projected_area_battery;
    double projected_area_package;
};

#endif // TOSN_ENERGY_H
