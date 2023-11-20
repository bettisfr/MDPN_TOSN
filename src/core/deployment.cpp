#include <random>
#include "deployment.h"

deployment::deployment(const input &par) {
    num_sensors = par.num_sensors;
    num_depots = par.num_depots;
    area_length = par.area_length;
    area_width = par.area_width;
    max_data = par.max_data;
    sensor_radius = par.sensor_radius;
    doi = par.doi;

    mt19937 re(par.seed);
    uniform_real_distribution<double> length_rand(0, area_length);
    uniform_real_distribution<double> width_rand(0, area_width);
    uniform_int_distribution<int> data_rand(0, max_data);

    double FSPL = get_FSPL();

    // Sensors creation
    for (int i = 0; i < num_sensors; i++) {
        double x = length_rand(re);
        double y = width_rand(re);
        double data = data_rand(re);

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

    // Depots creation
    for (int i = 0; i < num_depots; i++) {
        double x = length_rand(re);
        double y = width_rand(re);

        depots.emplace_back(x, y);
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

vector<point> deployment::get_depots() {
    return depots;
}

ostream &operator<<(ostream &os, const deployment &d) {
    for (int i = 0; i < d.sensors.size(); i++) {
        sensor s = d.sensors[i];
        cout << i << ": " << s << endl;
    }

    return os;
}

double deployment::get_FSPL() {
    double min_P_Rx = -100.0;
    double min_P_Tx = min_P_Rx + 20 * log10(4 * M_PI * f_c * sensor_radius / c);
    double FSPL = min_P_Tx - min_P_Rx;
    return FSPL;
}

double deployment::get_DTR(double distance) {
    double P_Rx = P_Tx + 20 * log10(c / (4 * M_PI * f_c * distance));
    double P_Rx_W = pow(10, P_Rx / 10.0);  // Convert dB to watts
    double C = B * log2(1 + P_Rx_W / N);
    double C_MBps = C / 8e6;  // Convert bps to MB/s
    return C_MBps;
}





