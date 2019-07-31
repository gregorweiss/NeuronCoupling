
#include <iostream>
#include <valarray>
#include <boost/numeric/odeint.hpp>

#include <fftw3.h>

#include "src/ArgParse.h"
#include "src/ode.h"
#include "src/npy.h"

double FHNfrequency(double eps1, double gamma1, double a1, double gamma2) {
    double tmp1(-1.0 * gamma1 * gamma2);
    double tmp2(2 * eps1);
    double tmp3(-1.0 * (1 - a1 * a1) * (1 - a1 * a1));
    double ret(tmp1 + tmp2 + tmp3);
    double sqr(tmp1 + tmp3);
    sqr *= (tmp1 + 2 * tmp2 + tmp3);
    ret += std::sqrt(sqr);
    ret /= 2.0;

    return ret;
}

using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;

int main(int argc, char *argv[]) {

    ArgParse &input(*new ArgParse(argc, argv));

    auto equilibrate = input.flag<int>("-eq", 1);
    auto eq_t = input.flag<double>("-eqt", 1000.0);
    auto t(input.flag<double>("-t", 500.0));
    auto dt(input.flag<double>("-dt", 0.01));
    auto nts(input.flag<size_t>("-nts", 1));
    auto adaptive(input.flag<int>("-adapt", 1));
    auto abserr(input.flag<double>("-abs", 1.0e-10));
    auto relerr(input.flag<double>("-rel", 1.0e-6));
    auto eps1(input.flag<double>("-eps1", 0.1));
    auto gamma1(input.flag<double>("-gamma1", 2.0));
    auto a1(input.flag<double>("-a1", 1.3));
    auto eps2(input.flag<double>("-eps2", 0.1));
    auto gamma2(input.flag<double>("-gamma2", -1.5));
    auto a2(input.flag<double>("-a2", 0.55));
    auto perturb(input.flag<double>("-e", 0.1));

    // Set the parameters for the FitzHugh-Nagumo model
    // Either read them from the parameter.dat
    // or simply hard code them here.
    vector<double> params(7, 0.0);
    params[0] = eps1;   // eps1
    params[1] = gamma1; // gamma1
    params[2] = a1;     // a1
    params[3] = eps2;   // eps2
    params[4] = gamma2; // gamma2
    params[5] = a2;     // a2

    cout << "# Frequency: " << FHNfrequency(eps1, gamma1, a1, gamma2) << endl;

    // Set inital values
    state_type init(4, 0.0);
    init[0] = -1.0 * params[2] + perturb;
    init[1] = params[2] * params[2] * params[2] / 3.0 - params[2] - params[1] * params[5];
    init[2] = -1.0 * params[5];
    init[3] = params[5] * params[5] * params[5] / 3.0 - params[5] - params[4] * params[2];

    vector<state_type> x_vec;
    vector<double> times;

    FitzHughNagumo fhn_obj(params);

    /* Various stepper tests */
    size_t steps(0);
    if (adaptive) {

        typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
        controlled_runge_kutta<runge_kutta_fehlberg78<state_type> > controlled_stepper;

        /* Let's first equilibrate */
        if (equilibrate) {
            vector<state_type> eq_x_vec;
            vector<double> eq_times;

            steps = integrate_adaptive(make_controlled<error_stepper_type>(abserr, relerr),
                                       fhn_obj, init,
                                       0.0, eq_t, dt,
                                       pusher(eq_x_vec, eq_times));

            init[0] = eq_x_vec[steps][0];
            init[1] = eq_x_vec[steps][1];
            init[2] = eq_x_vec[steps][2];
            init[3] = eq_x_vec[steps][3];
        }

        /* Run it */
        steps = integrate_adaptive(make_controlled<error_stepper_type>(abserr, relerr),
                                   fhn_obj, init,
                                   0.0, t, dt,
                                   pusher(x_vec, times));

    } else {

        std::cout << "Intergrating with non-adaptive stepper." << std::endl;
        runge_kutta4<state_type> stepper;

        /* Let's first equilibrate */
        if (equilibrate) {
            vector<state_type> eq_x_vec;
            vector<double> eq_times;

            steps = integrate_const(stepper, fhn_obj, init,
                                    0.0, eq_t, dt,
                                    pusher(eq_x_vec, eq_times));

            init = eq_x_vec[steps];
        }

        /* Run it */
        steps = integrate_const(stepper, fhn_obj, init,
                                0.0, t, dt,
                                pusher(x_vec, times));
    }

    /*
    Here we actually output in to python numpy output because the data is being post-processed in python.
    */
    vector<double> accum;
    for (auto &sub : x_vec)
        accum.insert(accum.end(), sub.begin(), sub.end());

    unsigned int n_dims;
    n_dims = 2;
    const long unsigned shape[] = {steps, init.size()};
    npy::SaveArrayAsNumpy<double>("output.npy", false, n_dims, shape, accum);

    return 0;
}
