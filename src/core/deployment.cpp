#include <random>
#include <algorithm>

#include "deployment.h"

deployment::deployment(const input &par) {
    num_sensors = par.num_sensors;
    num_depots = par.num_depots;
    area_length = par.area_length;
    area_width = par.area_width;
    max_data = par.max_data;
    sensor_radius = par.sensor_radius;
    sensor_radius_doi_percentage = par.sensor_radius_doi_percentage;
    doi = par.doi;
    epsilon = par.epsilon;

    energy_budget = par.energy_budget;
    energy_cons_fly = par.energy_cons_fly;
    energy_cons_hover = par.energy_cons_hover;

    wireless_technology = par.wireless_technology;
    if (wireless_technology == 0) {
        // WiFi-5
        f_c = 5e9;
        P_Tx = 20;
        N = 1e-9;
        B = 80e6;
    } else if (wireless_technology == 1) {
        // WiFi-4
        f_c = 2.4e9;
        P_Tx = 20;
        N = 1e-9;
        B = 20e6;
    } else if (wireless_technology == 2) {
        // Zigbee
        f_c = 2.4e9;
        P_Tx = 5;
        N = 1e-9;
        B = 2e6;
    } else if (wireless_technology == 3) {
        // Bluetooth
        f_c = 2.4e9;
        P_Tx = 0;
        N = 1e-9;
        B = 1e6;
    }

    if (!check_feasibility()) {
        cout << "Error: Feasibility error" << endl;
        exit(1);
    }

    static mt19937 re(par.seed);
    uniform_real_distribution<double> length_rand(0, area_length);
    uniform_real_distribution<double> width_rand(0, area_width);
    uniform_int_distribution<int> data_rand(0, max_data);

    double FSPL = get_FSPL();

    // Depots creation
    for (int i = 0; i < num_depots; i++) {
        double x = length_rand(re);
        double y = width_rand(re);

        depots.emplace_back(x, y);
    }

    // Sensors creation
    for (int i = 0; i < num_sensors; i++) {
        double x = length_rand(re);
        double y = width_rand(re);
        double data = data_rand(re);

        // pre-processing phase
        if (par.algorithm == 1){
            point dep = depots[0];
            depots.clear();
            depots.emplace_back(dep);
        }

        // Suitably fix this according to Degree of Irregularity (DOI)
        // https://www.sciencedirect.com/science/article/pii/S1574119218305406 (Section 4)
        // DOI=0 sphere, DOI>0 irregularities
        weibull_distribution<double> weibull_dist(shape, scale);

        vector<double> values(360);
        values[0] = -10; // dummy value for the next "while", so it enters for sure
        while (abs(values.back() - values.front()) > doi) {
            values[0] = 0; // fixed to 0
            for (int j = 1; j < 360; ++j) {
                values[j] = weibull_dist(re);
                uniform_int_distribution<int> sgn_dist(0, 1);  // Range [0, 1]
                int sgn = (sgn_dist(re) == 0) ? -1 : 1;
                values[j] = sgn * values[j] * doi;
            }
        }

        vector<double> r_doi(360);
        for (int angle = 0; angle < 360; angle++) {
            r_doi[angle] = sensor_radius * pow(10, (values[angle] * FSPL / 20.0));
        }

        sensors.emplace_back(x, y, data, r_doi);
    }

    // Sort sensors by x-coordinate
    sort(sensors.begin(), sensors.end(), [](const sensor &a, const sensor &b) {
        return a.get_pos_x() < b.get_pos_x();
    });

    for (int i = 0; i < num_sensors; i++) {
        sensors[i].set_id(i);
    }
}

int deployment::get_num_sensors() const {
    return num_sensors;
}

int deployment::get_area_length() const {
    return area_length;
}

int deployment::get_area_width() const {
    return area_width;
}

double deployment::get_sensor_radius() const {
    return sensor_radius;
}

vector<sensor> deployment::get_sensors() {
    return sensors;
}

double deployment::get_sensor_radius_doi() const {
    // Reduced by a percentage
    return sensor_radius_doi_percentage * sensor_radius;
}

vector<point> deployment::get_depots() {
    return depots;
}

double deployment::get_epsilon() const {
    return epsilon;
}

double deployment::get_doi() const {
    return doi;
}

ostream &operator<<(ostream &os, const deployment &d) {
    for (int i = 0; i < d.sensors.size(); i++) {
        sensor s = d.sensors[i];
        cout << i << ": " << s << endl;
    }

    return os;
}

double deployment::get_FSPL() const {
    double min_P_Rx = -100.0;
    double min_P_Tx = min_P_Rx + 20 * log10(4 * M_PI * f_c * sensor_radius / c);
    double FSPL = min_P_Tx - min_P_Rx;

    return FSPL;
}

double deployment::get_DTR(double distance) const {
    if (distance < 1) {
        distance = 1;
    }
    if (distance > sensor_radius) {
        distance = sensor_radius;
    }
    double P_Rx = P_Tx + 20 * log10(c / (4 * M_PI * f_c * distance));
    double P_Rx_W = pow(10, P_Rx / 10.0);  // Convert dB to watts
    double C = B * log2(1 + P_Rx_W / N);
    double C_MBps = C / 8e6;  // Convert bps to MB/s

    return C_MBps;
}

double deployment::get_energy_budget() const {
    return energy_budget;
}

double deployment::get_energy_cons_fly() const {
    return energy_cons_fly;
}

double deployment::get_energy_cons_hover() const {
    return energy_cons_hover;
}

bool deployment::check_feasibility() const {
    // flying
    double delta_x = area_length;
    double delta_y = area_width;
    double max_distance = sqrt(delta_x * delta_x + delta_y * delta_y) + sensor_radius;
    double energy_flying = max_distance * energy_cons_fly;

    // hovering
    double min_dtr = get_DTR(sensor_radius);
    double max_required_time = max_data / min_dtr;
    double energy_hovering = max_required_time * energy_cons_hover;

    double required_energy = 2*energy_flying + energy_hovering;
    if (required_energy > energy_budget) {
        cout << "Infeasible: required_energy=" << required_energy << " while budget=" << energy_budget << endl;
        return false;
    }
    return true;
}
