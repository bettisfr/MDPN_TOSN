#ifndef TOSN_DEFINITIONS_H
#define TOSN_DEFINITIONS_H

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

typedef tuple<double, double> point;

struct solution {
    // General info to be saved
    vector<vector<tuple<point, int>>> tours;
    int tours_number;
    vector<double> tours_costs;

    // For DOI
    int total_sensors = -1;
    int uncovered_sensors = -1;

    // For DTR
    double total_data = -1;
    double lost_data = -1;

    double running_time = -1;

    friend ostream& operator<<(ostream& os, const solution& out) {
        os << "Tours Number: " << out.tours_number << endl;
        os << "Tours Costs: ";
        for (const auto& cost : out.tours_costs) {
            os << fixed << setprecision(2) << cost << " ";
        }
        os << defaultfloat; // Reset to default precision mode
        os << endl;

        os << "Total Sensors: " << out.total_sensors << endl;
        os << "Uncovered Sensors: " << out.uncovered_sensors << endl;

        os << "Total Data: " << out.total_data << endl;
        os << "Lost Data: " << out.lost_data << endl;
        os << "Running Time: " << out.running_time << endl;

        return os;
    }
};

#endif // TOSN_DEFINITIONS_H
