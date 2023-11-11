#include <iostream>

#include "core/sensor.h"

using namespace std;

int main(int argc, char **argv) {

    sensor s(10, 20, 100, 50);
    cout << s << endl;

    return 0;
}
