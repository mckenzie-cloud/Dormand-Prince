
#include <iostream>
#include <array>              // std::array<>
#include <cmath>              // fmin(), fmax(), pow()
#include <chrono>
#include <functional>         // for function pointer.

constexpr size_t dimensions = 3;    // 3-dimensions (x, y, and z)

/* Dorman-Prince 4(5) butcher tableau*/
const double a0 = 1.0 / 5.0, a1 = 1.0 / 5.0;
const double b0 = 3.0 / 10.0, b1 = 3.0 / 40.0, b2 = 9.0 / 40.0;
const double c0 = 4.0 / 5.0, c1 = 44.0 / 45.0, c2 = -56.0 / 15.0, c3 = 32.0 / 9.0;
const double d0 = 8.0 / 9.0, d1 = 19372.0 / 6561.0, d2 = -25360.0 / 2187.0, d3 = 64448.0 / 6561.0, d4 = -212.0 / 729.0;
const double e0 = 1.0, e1 = 9017.0 / 3168.0, e2 = -355.0 / 33.0, e3 = 46732.0 / 5247.0, e4 = 49.0 / 176.0, e5 = -5103.0 / 18656.0;
const double f0 = 1.0, f1 = 35.0 / 384.0, f3 = 500.0 / 1113.0, f4 = 125.0 / 192.0, f5 = -2187.0 / 6784.0, f6 = 11.0 / 84.0;
const double g1 = 5179.0 / 57600.0, g3 = 7571.0 / 16695.0, g4 = 393.0 / 640.0, g5 = -92097.0 / 339200.0, g6 = 187.0 / 2100.0, g7 = 1.0 / 40.0;

static std::array<double, dimensions> Lorentz(double t, std::array<double, dimensions> &x)
{
	/* Parameters for Lorentz system. */
	double sigma = 10.0, beta = 8.0/3.0, rho = 24.7;
	std::array<double, dimensions> derivative{};
	derivative[0] = sigma * (x[1] - x[0]);              // dx
	derivative[1] = x[0] * (rho  - x[2]) - x[1];        // dy
	derivative[2] = x[0] * x[1] - beta * x[2];          // dz
	return derivative;
}

static double Error(std::array<double, dimensions>& yk_1, std::array<double, dimensions>& zk_1, 
    std::array<double, dimensions>& yk, double abs_tol, double rel_tol)
{
    /* We use L2-norm for error estimation. */
    double err = 0.0;
    for (size_t i = 0; i < dimensions; i++)
    {
        double diff = fabs(zk_1[i] - yk_1[i]);
        double max = fmax(fabs(yk[i]), fabs(yk_1[i]));
        double sci = abs_tol + max * rel_tol;
        err += (diff / sci) * (diff / sci);
    }
    err = sqrt(err / (double)dimensions);
    return err;
}

static void Integrate(double t0, double h, int T, double abs_tol, 
    double rel_tol, std::array<double, dimensions> &y0, 
    std::function<std::array<double, dimensions>(double, std::array<double, dimensions> &)> Func)
{ 

    /* Safety factors */
    double fac = 0.9, fac_min = 0.2, fac_max = 5.0;

	std::array<double, dimensions> yk{};

    // set initial value
    for (size_t i = 0; i < y0.size(); i++)
    {
        yk[i] = y0[i];
    }

    double dt = t0;
    int accepted = 0, rejected = 0;

    double total_steps = h;
    while (dt < T) {

        int is_rejected = 1;
        while (is_rejected)
        {
            std::array<double, dimensions> k{};
            // ------------------------------First-stage-----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i];
            }
            std::array<double, dimensions> n1 = Func(dt, k);

            // ------------------------------Second-stage----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (a1 * n1[i]);
            }
            std::array<double, dimensions> n2 = Func(dt, k);

            // ------------------------------Third-stage-----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (b1 * n1[i] + b2 * n2[i]);
            }
            std::array<double, dimensions> n3 = Func(dt, k);

            // ------------------------------Fourth-stage----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (c1 * n1[i] + c2 * n2[i] + c3 * n3[i]);
            }
            std::array<double, dimensions> n4 = Func(dt, k);

            // ------------------------------Fifth-stage-----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (d1 * n1[i] + d2 * n2[i] + d3 * n3[i] + d4 * n4[i]);
            }
            std::array<double, dimensions> n5 = Func(dt, k);

            // ------------------------------Sixth-stage--------------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (e1 * n1[i] + e2 * n2[i] + e3 * n3[i] + e4 * n4[i] + e5 * n5[i]);
            }
            std::array<double, dimensions> n6 = Func(dt, k);

            // ------------------------------Seventh-stage-----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (f1 * n1[i] + f3 * n3[i] + f4 * n4[i] + f5 * n5[i] + f6 * n6[i]);
            }
            std::array<double, dimensions> n7 = Func(dt, k);

            // ------------------------------------------------------------------------------
            std::array<double, dimensions> yk_1{};           /* Gives the fifth-order accurate solution */
            for (size_t i = 0; i < dimensions; i++)
            {
                yk_1[i] = k[i];
            }

            // ------------------------------------------------------------------------------
            std::array<double, dimensions> zk_1{};           /* Gives the fourth-order accurate solution */
            for (size_t i = 0; i < dimensions; i++)
            {
                zk_1[i] = yk[i] + h * (g1 * n1[i] + g3 * n3[i] + g4 * n4[i] + g5 * n5[i] + g6 * n6[i] + g7 * n7[i]);
            }

            double err = Error(yk_1, zk_1, yk, abs_tol, rel_tol);
            double h_new = 0.0;
            if (err <= 1.0) {
                for (size_t i = 0; i < yk.size(); i++)
                {
                    yk[i] = yk_1[i];
                }
                h_new = h * fmin(fac_max, fmax(fac_min, fac * pow(1.0 / err, 1.0 / 5.0)));
                h = h_new;
                dt += h;
                is_rejected = 0;
                accepted++;
            }
            else {
                h_new = h * fmin(1.0, fmax(fac_min, fac * pow(1.0 / err, 1.0 / 5.0)));
                h = h_new;
                rejected++;
            }
        }
    }
    std::cout << "Accepted: " << accepted << " Rejected: " << rejected << '\n';
}

/*
* In some applications dense output is required and in such a case RK methods must frequently shorten the step-size 
* and are therefore inefficient. Thus the new RK solvers should have the possibility of producing, if necessary, 
* reliable approximations to the solution at any point of the integration interval without step-size adjustment 
* and with little additional computational cost. This has been the main reason 
* for developing the so-called continuous or interpolatory RK methods. 
* 
* This implementation is an example of A continuous formula of order 4.
* See: Hairer, Norsett, Wanner: Solving Ordinary Differential Equations, Nonstiff Problems. I, pp.191-192
*/
void CalculateContinuousExtension(double t, double t_old, double h, std::array<double, dimensions> &derived_old,
    std::array<double, dimensions> &k3, std::array<double, dimensions> &k4, std::array<double, dimensions> &k5,
    std::array<double, dimensions> &k6, std::array<double, dimensions> &derived_new, std::array<double, dimensions> &y0,
    std::array<double, dimensions> &out)
{

    // k1 = &derive_old;
    // k3 = &k3;
    // k4 = &k4;
    // k5 = &k5;
    // k6 = &k6;
    // k7 = &derived_new;

    /* Evaluate interpolating polynomial for y[i] at location x, where t_old <= x <= t_old + h. */
    double theta = (t - t_old) / h;
    /* Continuous extension parameters. */
    const double X1 = 5.0 * (2558722523.0 - 31403016.0 * theta) / 11282082432.0;
    const double X3 = 100.0 * (882725551.0 - 15701508.0 * theta) / 32700410799.0;
    const double X4 = 25.0 * (443332067.0 - 31403016.0 * theta) / 1880347072.0;
    const double X5 = 32805.0 * (23143187.0 - 3489224.0 * theta) / 199316789632.0;
    const double X6 = 55.0 * (29972135.0 - 7076736.0 * theta) / 822651844.0;
    const double x7 = 10.0 * (7414447.0 - 829305.0 * theta) / 29380423.0;

    double theta_sqr = theta * theta;
    double two_theta = 2.0 * theta;
    double common_term1 = theta_sqr * (3.0 - two_theta);
    double common_term2 = theta_sqr * (theta_sqr - two_theta + 1.0);

    double b1_theta = common_term1 * f1 + theta * (theta_sqr - two_theta + 1.0) - common_term2 * X1;
    double b3_theta = common_term1 * f3 + common_term2 * X3;
    double b4_theta = common_term1 * f4 - common_term2 * X4;
    double b5_theta = common_term1 * f5 + common_term2 * X5;
    double b6_theta = common_term1 * f6 - common_term2 * X6;
    double b7_theta = theta_sqr * (theta - 1.0) + common_term2 * x7;

    for (size_t i = 0; i < dimensions; i++)
    {
        out[i] = y0[i] + h * (b1_theta * derived_old[i] + b3_theta * k3[i] + b4_theta * k4[i] + b5_theta * k5[i] + b6_theta * k6[i] + b7_theta * derived_new[i]);
    }
}

int main(void)
{
    /* Set initial value. */
    std::array<double, dimensions> y0 = { -8.0, 8.0, 27.0 };
    double t0 = 0.0, h0 = 1e-3, abs_tol = 1e-6, rel_tol = 1e-3;
    int T = 20000;
    Integrate(t0, h0, T, abs_tol, rel_tol, y0, Lorentz);
	return 0;
}