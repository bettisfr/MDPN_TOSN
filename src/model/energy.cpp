#include "energy.h"
#include "quartic.h"

#include <cmath>

energy::energy() {
    m_drone = 3.6;
    m_battery = 2.7;

    num_rotors = 4;
    diameter = 0.895;

    pressure = 100726;
    R = 287.058;
    temperature = 15 + 273.15;
    rho = pressure / (R * temperature);

    g = 9.81;

    eta = 0.7;

    drag_coefficient_drone = 1.49;
    drag_coefficient_battery = 1;
    drag_coefficient_package = 2.2;

    projected_area_drone = 0.224;
    projected_area_battery = 0.015;
    projected_area_package = 0.0929;
}

double energy::get_energy_hovering(double time, double payload_weight, double wind_speed) {
    double m_package = payload_weight;

    double v_air = wind_speed;

    // Drag force
    double F_drag_drone = 0.5 * rho * pow(v_air, 2) * drag_coefficient_drone * projected_area_drone;
    double F_drag_battery = 0.5 * rho * pow(v_air, 2) * drag_coefficient_battery * projected_area_battery;
    double F_drag_package = 0.5 * rho * pow(v_air, 2) * drag_coefficient_package * projected_area_package;

    double F_drag = F_drag_drone + F_drag_battery + F_drag_package;

    // Thrust
    double T = (m_drone + m_battery + m_package) * g + F_drag;

    // Power min hover
    double P_min_hover = pow(T, 1.5) / (sqrt(0.5 * M_PI * num_rotors * pow(diameter, 2) * rho));

    // Expended power
    double P = P_min_hover / eta;

    // Energy consumed
    double E = P * time;

    return E;
}

double energy::get_energy_movement(double distance, double payload_weight, double drone_speed, double wind_speed,
                                   double relative_wind_direction) {
    double m_package = payload_weight;

    double v_north = drone_speed - wind_speed * cos(relative_wind_direction * M_PI / 180);
    double v_east = -wind_speed * sin(relative_wind_direction * M_PI / 180);
    double v_air = sqrt(pow(v_north, 2) + pow(v_east, 2));

    // Drag force
    double F_drag_drone = 0.5 * rho * pow(v_air, 2) * drag_coefficient_drone * projected_area_drone;
    double F_drag_battery = 0.5 * rho * pow(v_air, 2) * drag_coefficient_battery * projected_area_battery;
    double F_drag_package = 0.5 * rho * pow(v_air, 2) * drag_coefficient_package * projected_area_package;

    double F_drag = F_drag_drone + F_drag_battery + F_drag_package;

    double alpha = atan(F_drag / ((m_drone + m_battery + m_package) * g));

    // Thrust
    double T = (m_drone + m_battery + m_package) * g + F_drag;

    double tmp_a = 2 * T;
    double tmp_b = M_PI * num_rotors * pow(diameter, 2) * rho;
    double tmp_c = pow(drone_speed * cos(alpha), 2);
    double tmp_d = drone_speed * sin(alpha);
    double tmp_e = tmp_a / tmp_b;

    complex<double> coeff[5] = {
            complex<double>(1, 0),
            complex<double>(2 * tmp_d, 0),
            complex<double>(tmp_c + pow(tmp_d, 2), 0),
            complex<double>(0, 0),
            complex<double>(-pow(tmp_e, 2), 0),
    };
    complex<double> sol[4];

    solve_quartic(coeff, sol);

    double maxRealPart = sol[0].real();  // Initialize with the real part of the first element

    for (int i = 1; i < 4; ++i) {
        double currentRealPart = sol[i].real();
        if (currentRealPart > maxRealPart) {
            maxRealPart = currentRealPart;
        }
    }

    // Now maxRealPart contains the maximum real part among the elements of sol
    double induced_speed = maxRealPart;

    // Power min to go forward
    double P_min = T * (drone_speed * sin(alpha) + induced_speed);

    // Expended power
    double P = P_min / eta;

    // Energy efficiency of travel
    double mu = P / drone_speed;

    // Energy consumed
    double E = mu * distance;

    return E;
}
