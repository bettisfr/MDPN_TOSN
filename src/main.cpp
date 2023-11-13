#include <iostream>

#include "core/sensor.h"
#include "core/deployment.h"

using namespace std;

int main(int argc, char **argv) {

    int n = 100;
    int length = 1000;
    int width = 1000;
    double radius = 150;
    double doi = 1;

    deployment dep(n, length, width, radius, doi);

    return 0;
}
