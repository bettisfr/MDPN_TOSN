#ifndef TOSN_OUTPUT_H
#define TOSN_OUTPUT_H

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>
#include <fstream>
#include <cmath>
#include <vector>

#include "input.h"
#include "definitions.h"

using namespace std;

double calculate_average(const vector<double>&);
double calculate_stddev(const vector<double>&, double);
void save_output(const input &par, vector<solution> &);

#endif //TOSN_OUTPUT_H
