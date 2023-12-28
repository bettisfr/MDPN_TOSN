#ifndef TOSN_DEFINITIONS_H
#define TOSN_DEFINITIONS_H

#include <iostream>

using namespace std;

typedef tuple<double, double> point;

//struct solution {
//    vector<vector<tuple<point, int>>> tours;
//    vector<double> tours_costs;
//
//    friend ostream& operator<<(ostream& os, const solution& sol) {
//        os << "Number of tours: " << sol.tours.size() << endl;
//        for (int i = 0; i < sol.tours.size(); i++) {
////            for (auto j: sol.tours[i]) {
////                os << "(" << get<0>(get<0>(j)) << ", " << get<1>(get<0>(j)) << ") - ID=" << get<1>(j) << endl;
////            }
//            os << "Tour cost[" << (i+1) << "]: " << sol.tours_costs[i] << endl;
////            os << "-------" << endl;
//        }
//        return os;
//    }
//};

struct output {
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

    friend ostream& operator<<(ostream& os, const output& out) {
        os << "Tours Number: " << out.tours_number << "\n";
        os << "Tours Costs: ";
        for (const auto& cost : out.tours_costs) {
            os << cost << " ";
        }
        os << "\n";

        os << "Total Sensors: " << out.total_sensors << "\n";
        os << "Uncovered Sensors: " << out.uncovered_sensors << "\n";

        os << "Total Data: " << out.total_data << "\n";
        os << "Lost Data: " << out.lost_data << "\n";

        return os;
    }
};

#endif // TOSN_DEFINITIONS_H
