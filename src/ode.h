//
// Created by gregor on 22.01.19.
//

#ifndef FHN_ODE_H
#define FHN_ODE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/scoped_ptr.hpp>
#include <vector>
#include <valarray>

using namespace std;

typedef vector<double> state_type;

struct pusher {
    std::vector<state_type> &m_states;
    std::vector<double> &m_times;

    pusher(std::vector<state_type> &states, std::vector<double> &times)
            : m_states(states), m_times(times) {}

    void operator()(const state_type &x, double t) {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

/////////////////////////////////
// class ODE
/////////////////////////////////

class ODE {
public:
    ODE() : _params(std::vector<double>()) {}

    ODE(const std::vector<double> &p) : _params(p) {}

    virtual ~ODE() {};

    unsigned int size();

protected:

    double get_param(size_t n);

    unsigned int _size;

private:

    const std::vector<double> _params;

};

typedef ODE *ode_ptr;

class FitzHughNagumo : public ODE {
public:

    FitzHughNagumo() : ODE() {}

    FitzHughNagumo(const std::vector<double> &p) : ODE(p),
                                                   _eps1(this->get_param(0)),
                                                   _gamma1(this->get_param(1)),
                                                   _a1(this->get_param(2)),
                                                   _eps2(this->get_param(3)),
                                                   _gamma2(this->get_param(4)),
                                                   _a2(this->get_param(5)),
                                                   _third(1.0 / 3.0) { this->_size = 4; }

    ~FitzHughNagumo() {};

    void operator()(const vector<double> &, vector<double> &, const double);

private:

    double _eps1;
    double _gamma1;
    double _a1;
    double _eps2;
    double _gamma2;
    double _a2;

    double _third;

};

#endif //FHN_ODE_H
