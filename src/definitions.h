#ifndef TOSN_DEFINITIONS_H
#define TOSN_DEFINITIONS_H

#include <iostream>

using namespace std;

typedef tuple<double, double> point;

struct solution {
    vector<vector<tuple<point, int>>> tours;
    vector<double> costs;

    friend ostream& operator<<(ostream& os, const solution& sol) {
        os << "Number of tours: " << sol.tours.size() << endl;
        for (int i = 0; i < sol.tours.size(); i++) {
//            for (auto j: sol.tours[i]) {
//                os << "(" << get<0>(get<0>(j)) << ", " << get<1>(get<0>(j)) << ") - ID=" << get<1>(j) << endl;
//            }
            os << "Tour cost[" << (i+1) << "]: " << sol.costs[i] << endl;
//            os << "-------" << endl;
        }
        return os;
    }
};

#endif // TOSN_DEFINITIONS_H
