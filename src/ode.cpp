//#include <utility>

//
// Created by gregor on 22.01.19.
//

#include "ode.h"

using namespace std;

double ODE::get_param(size_t n) { return _params.at(n); }

unsigned int ODE::size() { return _size; }

void FitzHughNagumo::operator()(const vector<double> &input, vector<double> &output, const double) {
    double x1(input[0]);
    double y1(input[1]);
    double x2(input[2]);
    double y2(input[3]);

    //valarray<double> output( 0.0, _size );
    output[0] = x1 - _third * (x1 * x1 * x1) - y1 + _gamma1 * x2;
    output[1] = _eps1 * (x1 + _a1);
    output[2] = x2 - _third * (x2 * x2 * x2) - y2 + _gamma2 * x1;
    output[3] = _eps2 * (x2 + _a2);

}
