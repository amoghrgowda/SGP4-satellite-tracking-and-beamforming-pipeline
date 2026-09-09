#include "SGP4Propagator.hpp"
#include <stdexcept>
#include <algorithm> // for std::clamp, which was introduced in cpp17 and snaps a value to particular range.
                     // useful for constraining the eccentricity (degree of deviation of our orbit from perfect circle) to a minimum>0, and max<1.
                     // <=0 means 'divide by zero' error. >1 means a parabolic or hyperbolic trajectory, and our newton-raphson formula goes into infinite loop.
#include <cmath>
#include <numbers>

namespace AstroStuff{
SGP4Propagator::SGP4Propagator(const TLE& tle, const GravitationalConstants& constants) : tle_(tle), consts_(constants){
    initialize();
}
// next, we initialize the required fields (as per the specification in STR-3).
void SGP4Propagator::initialize(){
    constexpr double ck2 = 0.5 * 0.001082616; // half of J2
    constexpr double ck4 = 0.375 * 0.00000165597; // -3/8 * J4
    constexpr double twoThirds = 2.0/3.0;

    double a1 = std::pow(consts_.xke / tle_.meanMotion, twoThirds);
    double cosI = std::cos(tle_.inclination);
    double theta2 = cosI * cosI;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double beta02 = 1.0 - tle_.eccentricity * tle_.eccentricity;
    double beta0 = std::sqrt(beta02);

    double del1 = 1.5 * (ck2 / (a1 * a1)) * (x3thm1 / (beta0 * beta02));
    double a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + 134.0 / 81.0 * del1)));
    double del0 = 1.5 * (ck2 / (a0 * a0)) * (x3thm1 / (beta0 * beta02));

    n0_double_prime_ = tle_.meanMotion / (1.0 + del0);
    a0_double_prime_ = a0 / (1.0 - del0);

    perigee_ = (a0_double_prime_ * (1.0 - tle_.eccentricity) - 1.0) * consts_.radiusEarth;

    // Atmospheric density scaling parameters
    double s = 78.0 / consts_.radiusEarth + 1.0;
    double qoms2t = 120.0 / consts_.radiusEarth;

    //depending on the perigee range, we set the imaginary reference point 's'. 
    // To prevent division by 0 if the height of satellite falls below 98 km, we simply fix s to 20
    if (perigee_ < 156.0) {
        s = perigee_ - 78.0;
        if (perigee_ <= 98.0) s = 20.0;
        s = s / consts_.radiusEarth + 1.0;
    }
    s4_ = s;
    qoms24_ = std::pow((qoms2t - s4_ + 1.0), 4.0);

    double pinvsq = 1.0 / (a0_double_prime_ * a0_double_prime_ * beta02 * beta02);
    tsi_ = 1.0 / (a0_double_prime_ - s4_);
    eta_ = a0_double_prime_ * tle_.eccentricity * tsi_;
    double eeta = tle_.eccentricity * eta_;
    double psisq = std::abs(1.0 - eta_ * eta_);
    double coef = qoms24_ * std::pow(tsi_, 4.0);
    double coef1 = coef / std::pow(psisq, 3.5);

    c2_ = coef1 * n0_double_prime_ * (a0_double_prime_ * (1.0 + 1.5 * eta_ * eta_ + eeta * (4.0 + eta_ * eta_))
          + 0.75 * ck2 * tsi_ / psisq * x3thm1 * (8.0 + 3.0 * eta_ * eta_ * (8.0 + eta_ * eta_)));
    c1_ = tle_.bstar * c2_;
    c3_ = (tle_.eccentricity > 1e-4) ? (coef * tsi_ * consts_.j3oj2 * n0_double_prime_ * std::sin(tle_.inclination) / tle_.eccentricity) : 0.0;
    c4_ = 2.0 * n0_double_prime_ * coef1 * a0_double_prime_ * beta02 * (eta_ * (2.0 + 0.5 * eta_ * eta_) + tle_.eccentricity * (0.5 + 2.0 * eta_ * eta_)
          - 2.0 * ck2 * tsi_ / (a0_double_prime_ * psisq) * (-3.0 * x3thm1 * (1.0 - 2.0 * eeta + eta_ * eta_ * (1.5 - 0.5 * eeta))
          + 0.75 * (1.0 - theta2) * (2.0 * eta_ * eta_ - eeta * (1.0 + eta_ * eta_)) * std::cos(2.0 * tle_.argPerigee)));
    c5_ = 2.0 * coef1 * a0_double_prime_ * beta02 * (1.0 + 2.75 * (eta_ * eta_ + eeta) + eeta * eta_ * eta_);

    // Secular rates
    double temp1 = 3.0 * ck2 * pinvsq * n0_double_prime_;
    double temp2 = temp1 * ck2 * pinvsq;
    double temp3 = 1.25 * ck4 * pinvsq * pinvsq * n0_double_prime_; // J4 term

    mDot_ = n0_double_prime_ + 0.5 * temp1 * beta0 * x3thm1 + 0.0625 * temp2 * beta0 * (13.0 - 78.0 * theta2 + 137.0 * theta2 * theta2);
    omegaDot_ = -0.5 * temp1 * (1.0 - 5.0 * theta2) + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta2 * theta2) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta2 * theta2);
    omegamDot_ = -temp1 * cosI + 0.0625 * temp2 * (4.0 * cosI - 19.0 * cosI * theta2) + 2.0 * temp3 * cosI * (3.0 - 7.0 * theta2);

    d2_ = 4.0 * a0_double_prime_ * tsi_ * c1_ * c1_;
    d3_ = d2_ * tsi_ * c1_ / 3.0;
    d4_ = d3_ * tsi_ * c1_ * 0.25;
}

double SGP4Propagator::solveKepler(double M, double e, double tol, int maxIter) {
    // Normalise M to [-pi, pi]
    M = std::fmod(M, 2.0 * std::numbers::pi);
    if (M < -std::numbers::pi) M += 2.0 * std::numbers::pi;
    if (M >  std::numbers::pi) M -= 2.0 * std::numbers::pi;

    double E = M;
    for (int i = 0; i < maxIter; ++i) {
        double f = E - e * std::sin(E) - M;
        if (std::abs(f) < tol) return E;
        double fPrime = 1.0 - e * std::cos(E);
        E -= f / fPrime;
    }
    return E;
}

StateVector SGP4Propagator::propagate(double tsince) const {
    constexpr double ck2 = 0.5 * 0.001082616;

    // Update secular terms
    double xmdf = tle_.meanAnomaly + mDot_ * tsince;
    double omgadf = tle_.argPerigee + omegaDot_ * tsince;
    double xnoddf = tle_.raan + omegamDot_ * tsince;

    double tsq = tsince * tsince;
    double xnode = xnoddf;

    double tempa = 1.0 - (c1_ * tsince + d2_ * tsq + d3_ * tsq * tsince + d4_ * tsq * tsq);
    double tempe = tle_.bstar * (c4_ * tsince + c5_ * (std::sin(xmdf) - std::sin(tle_.meanAnomaly)));
    double templ = 1.5 * c1_ * tsq;

    double a = a0_double_prime_ * tempa * tempa;
    double e = tle_.eccentricity - tempe;
    e = std::clamp(e, 1e-6, 0.999999);

    double beta2 = 1.0 - e * e;

    //Long period periodic terms (Lyddane)
    double sinI = std::sin(tle_.inclination);
    double cosI = std::cos(tle_.inclination);

    double axn = e * std::cos(omgadf);
    double temp_lp = 0.5 * consts_.j3oj2 * sinI / (a * beta2);
    double ayn = e * std::sin(omgadf) - temp_lp;
    double xl = xmdf + omgadf + xnode + n0_double_prime_ * templ 
                - 0.25 * consts_.j3oj2 * sinI / (a * beta2) * axn * (3.0 + 5.0 * cosI) / (1.0 + cosI);

    // Solve Kepler equation for Mean Anomaly
    // Mean anomaly recovered from mean longitude:
    double u_mean = std::fmod(xl - xnode, 2.0 * std::numbers::pi);
    if (u_mean < 0.0) u_mean += 2.0 * std::numbers::pi;

    double epw = u_mean;
    for (int i = 0; i < 15; ++i) {
        double sinEpw = std::sin(epw);
        double cosEpw = std::cos(epw);
        double f = epw - axn * sinEpw + ayn * cosEpw - u_mean;
        double fPrime = 1.0 - axn * cosEpw - ayn * sinEpw;
        double delta = f / fPrime;
        delta = std::clamp(delta, -0.95, 0.95);
        epw -= delta;
        if (std::abs(delta) < 1e-12) break;
    }

    double sinEpw = std::sin(epw);
    double cosEpw = std::cos(epw);

    double ecose = axn * cosEpw + ayn * sinEpw;
    double esine = axn * sinEpw - ayn * cosEpw;
    double el2 = axn * axn + ayn * ayn;
    double pl = a * (1.0 - el2);

    double rl = a * (1.0 - ecose);
    double betal = std::sqrt(1.0 - el2);

    // True anomaly & Argument of Latitude
    double temp_sin = esine / (1.0 + betal);
    double sinu = (a / rl) * (sinEpw - ayn - axn * temp_sin);
    double cosu = (a / rl) * (cosEpw - axn + ayn * temp_sin);
    double u = std::atan2(sinu, cosu); // True argument of latitude

    // Short period perturbations (J2)
    double sin2u = 2.0 * sinu * cosu;
    double cos2u = 1.0 - 2.0 * sinu * sinu;

    double temp_p = 1.0 / pl;
    double temp1 = ck2 * temp_p;
    double temp2 = temp1 * temp_p;

    double x3thm1_val = 3.0 * cosI * cosI - 1.0;
    double x7thm1_val = 7.0 * cosI * cosI - 1.0;
    double x1mth2_val = 1.0 - cosI * cosI;

    double rk = rl * (1.0 - 1.5 * temp2 * betal * x3thm1_val) + 0.5 * temp1 * x1mth2_val * cos2u;
    double uk = u - 0.25 * temp2 * x7thm1_val * sin2u;
    double xnodek = xnode + 1.5 * temp2 * cosI * sin2u;
    double xinck = tle_.inclination + 1.5 * temp2 * cosI * sinI * cos2u;

    //Unit vectors in TEME orbital plane
    double sinUk = std::sin(uk);
    double cosUk = std::cos(uk);
    double sinNodek = std::sin(xnodek);
    double cosNodek = std::cos(xnodek);
    double sinInck = std::sin(xinck);
    double cosInck = std::cos(xinck);

    // Satellite position unit vector
    Vector3D U{
        -sinNodek * cosInck * sinUk + cosNodek * cosUk,
         cosNodek * cosInck * sinUk + sinNodek * cosUk,
         sinInck * sinUk
    };

    // Satellite in-plane transverse unit vector
    Vector3D V{
        -sinNodek * cosInck * cosUk - cosNodek * sinUk,
         cosNodek * cosInck * cosUk - sinNodek * sinUk,
         sinInck * cosUk
    };

    // Output state vector
    StateVector state;
    state.position = {
        rk * U.x * consts_.radiusEarth,
        rk * U.y * consts_.radiusEarth,
        rk * U.z * consts_.radiusEarth
    };

    double rdotl = consts_.xke * std::sqrt(a) * esine / rl;
    double rvdotl = consts_.xke * std::sqrt(pl) / rl;

    double rdot = rdotl - n0_double_prime_ * temp1 * x1mth2_val * sin2u;
    double rvdot = rvdotl + n0_double_prime_ * temp1 * (x1mth2_val * cos2u + 1.5 * x3thm1_val);
    double vFactor = consts_.radiusEarth / 60.0;

    state.velocity = {
        (rdot * U.x + rvdot * V.x) * vFactor,
        (rdot * U.y + rvdot * V.y) * vFactor,
        (rdot * U.z + rvdot * V.z) * vFactor
    };

    return state;
}

} // End AstroStuff namespace