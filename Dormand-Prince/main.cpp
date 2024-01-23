
#include <iostream>
#include <array>              // std::array<>
#include <cmath>              // fmin(), fmax(), pow()
#include <chrono>

constexpr size_t dimensions = 3;    // 3-dimensions (x, y, and z)

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

static void Integrate(double t0, double h, int T, std::array<double, dimensions> &y0, double abs_tol, double rel_tol)
{
    /* Dormand-Prince butcher tableau*/
    double a0 = 1.0 / 5.0,  a1 = 1.0 / 5.0;
    double b0 = 3.0 / 10.0, b1 = 3.0 / 40.0,       b2 =  9.0 / 40.0;
    double c0 = 4.0 / 5.0,  c1 = 44.0 / 45.0,      c2 = -56.0 / 15.0,      c3 = 32.0 / 9.0;
    double d0 = 8.0 / 9.0,  d1 = 19372.0 / 6561.0, d2 = -25360.0 / 2187.0, d3 = 64448.0 / 6561.0, d4 = -212.0 / 729.0;
    double e0 = 1.0,        e1 = 9017.0 / 3168.0,  e2 = -355.0 / 33.0,     e3 = 46732.0 / 5247.0, e4 = 49.0 / 176.0,  e5 = -5103.0 / 18656.0;
    double f0 = 1.0,        f1 = 35.0 / 384.0,     f3 = 500.0 / 1113.0,    f4 = 125.0 / 192.0, f5 = -2187.0 / 6784.0, f6 = 11.0 / 84.0;
    double g1 = 5179.0 / 57600.0, g3 = 7571.0 / 16695.0, g4 = 393.0 / 640.0, g5 = -92097.0 / 339200.0, g6 = 187.0 / 2100.0, g7 = 1.0 / 40.0;
    

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
            std::array<double, dimensions> n1 = Lorentz(dt, k);

            // ------------------------------Second-stage----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (a1 * n1[i]);
            }
            std::array<double, dimensions> n2 = Lorentz(dt, k);

            // ------------------------------Third-stage-----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (b1 * n1[i] + b2 * n2[i]);
            }
            std::array<double, dimensions> n3 = Lorentz(dt, k);

            // ------------------------------Fourth-stage----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (c1 * n1[i] + c2 * n2[i] + c3 * n3[i]);
            }
            std::array<double, dimensions> n4 = Lorentz(dt, k);

            // ------------------------------Fifth-stage-----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (d1 * n1[i] + d2 * n2[i] + d3 * n3[i] + d4 * n4[i]);
            }
            std::array<double, dimensions> n5 = Lorentz(dt, k);

            // ------------------------------Sixth-stage--------------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (e1 * n1[i] + e2 * n2[i] + e3 * n3[i] + e4 * n4[i] + e5 * n5[i]);
            }
            std::array<double, dimensions> n6 = Lorentz(dt, k);

            // ------------------------------Seventh-stage-----------------------
            for (size_t i = 0; i < dimensions; i++) {
                k[i] = yk[i] + h * (f1 * n1[i] + f3 * n3[i] + f4 * n4[i] + f5 * n5[i] + f6 * n6[i]);
            }
            std::array<double, dimensions> n7 = Lorentz(dt, k);

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

int main(void)
{
    /* Set initial value. */
    std::array<double, dimensions> y0 = { -8.0, 8.0, 27.0 };
    double t0 = 0.0, h0 = 1e-3, abs_tol = 1e-6, rel_tol = 1e-3;
    int T = 20000;
    Integrate(t0, h0, T, y0, abs_tol, rel_tol);
	return 0;
}