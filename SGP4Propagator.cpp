#include "SGP4Propagator.hpp"
#include <stdexcept>
#include <algorithm> // for std::clamp, which was introduced in cpp17 and snaps a value to particular range.
                     // useful for constraining the eccentricity (degree of deviation of our orbit from perfect circle) to a minimum>0, and max<1.
                     // <=0 means 'divide by zero' error. >1 means a parabolic or hyperbolic trajectory, and our newton-raphson formula goes into infinite loop.

namespace AstroStuff{
SGP4Propagator::SGP4Propagator(const TLE& tle, const GravitationalConstants& constants) : tle_(tle), consts_(constants){
    initialize();
}
// next, we initialize the required fields (as per the specification in STR-3).
void SGP4Propagator::initialize(){
    constexpr double ck2 = 0.5 * 0.001082616; // half of J2
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
    c3_ = (tle_.eccentricity > 1e-4) ? (coef * tsi_ * consts_.j3oj2 * n0_double_prime_ * consts_.radiusEarth * std::sin(tle_.inclination) / tle_.eccentricity) : 0.0;
    c4_ = 2.0 * n0_double_prime_ * coef1 * a0_double_prime_ * beta02 * (eta_ * (2.0 + 0.5 * eta_ * eta_) + tle_.eccentricity * (0.5 + 2.0 * eta_ * eta_)
          - 2.0 * ck2 * tsi_ / (a0_double_prime_ * psisq) * (-3.0 * x3thm1 * (1.0 - 2.0 * eeta + eta_ * eta_ * (1.5 - 0.5 * eeta))
          + 0.75 * (1.0 - theta2) * (2.0 * eta_ * eta_ - eeta * (1.0 + eta_ * eta_)) * std::cos(2.0 * tle_.argPerigee)));
    c5_ = 2.0 * coef1 * a0_double_prime_ * beta02 * (1.0 + 2.75 * (eta_ * eta_ + eeta) + eeta * eta_ * eta_);

    // Secular rates
    double temp1 = 3.0 * ck2 * pinvsq * n0_double_prime_;
    double temp2 = temp1 * ck2 * pinvsq;
    double temp3 = 1.25 * 0.00000165597 * pinvsq * pinvsq * n0_double_prime_; // J4 term

    mDot_ = n0_double_prime_ + 0.5 * temp1 * beta0 * x3thm1 + 0.0625 * temp2 * beta0 * (13.0 - 78.0 * theta2 + 137.0 * theta2 * theta2);
    omegaDot_ = -0.5 * temp1 * (1.0 - 5.0 * theta2) + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta2 * theta2) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta2 * theta2);
    omegamDot_ = -temp1 * cosI + 0.0625 * temp2 * (4.0 * cosI - 19.0 * cosI * theta2) + 2.0 * temp3 * cosI * (3.0 - 7.0 * theta2);

    d2_ = 4.0 * a0_double_prime_ * tsi_ * c1_ * c1_;
    d3_ = d2_ * tsi_ * c1_ / 3.0;
    d4_ = d3_ * tsi_ * c1_ * 0.25;
}

double SGP4Propagator::solveKepler(double M, double e, double tol, int maxIter) {
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
    // Secular updates
    double xmdf = tle_.meanAnomaly + mDot_ * tsince;
    double omgadf = tle_.argPerigee + omegaDot_ * tsince;
    double xnoddf = tle_.raan + omegamDot_ * tsince;

    double tsq = tsince * tsince;
    double xnode = xnoddf + d2_ * tsq + d3_ * tsq * tsince + d4_ * tsq * tsq;
    double tempa = 1.0 - (c1_ * tsince + d2_ * tsq + d3_ * tsq * tsince);
    double tempe = tle_.bstar * (c4_ * tsince + c5_ * (std::sin(xmdf) - std::sin(tle_.meanAnomaly)));
    double templ = 1.5 * c1_ * tsq; // secular drag L correction

    double a = a0_double_prime_ * tempa * tempa;
    double e = tle_.eccentricity - tempe;
    e = std::clamp(e, 1e-6, 0.999999);

    double xl = xmdf + omgadf + xnode + n0_double_prime_ * templ;
    double beta = std::sqrt(1.0 - e * e);
    double n = consts_.xke / std::pow(a, 1.5);

    // long period periodic corrections (Lyddane)
    double axn = e * std::cos(omgadf);
    double ayn = e * std::sin(omgadf) - 0.5 * consts_.j3oj2 * std::sin(tle_.inclination) / (a * beta * beta);
    double xl_long = xl - 0.25 * consts_.j3oj2 * std::sin(tle_.inclination) * axn * (3.0 + 5.0 * std::cos(tle_.inclination)) / (1.0 + std::cos(tle_.inclination));

    // Solving Kepler for eccentric anomaly
    double u = std::fmod(xl_long - xnode, 2.0 * std::numbers::pi);
    double E = solveKepler(u, e);

    // Short period perturbations & State Vector projection
    double sinE = std::sin(E);
    double cosE = std::cos(E);
    double ecosE = e * cosE;
    double esinE = e * sinE;

    double r = a * (1.0 - ecosE);
    double rDot = consts_.xke * std::sqrt(a) / r * esinE;
    double rfDot = consts_.xke * std::sqrt(a * (1.0 - e * e)) / r;

    // True anomaly & argument of latitude
    double sinv = (std::sqrt(1.0 - e * e) * sinE) / (1.0 - ecosE);
    double cosv = (cosE - e) / (1.0 - ecosE);
    double v = std::atan2(sinv, cosv);
    double u_lat = v + omgadf;

    // Unit vectors in TEME orbital plane
    double sin2u = std::sin(2.0 * u_lat);
    double cos2u = std::cos(2.0 * u_lat);
    double mr = r * (1.0 - 1.5 * 0.000541308 * (1.0 / (a * beta * beta)) * (3.0 * std::cos(tle_.inclination) * std::cos(tle_.inclination) - 1.0))
              + 0.5 * 0.000541308 * (1.0 / (a * beta * beta)) * std::sin(tle_.inclination) * std::sin(tle_.inclination) * cos2u;

    double sinu = std::sin(u_lat);
    double cosu = std::cos(u_lat);
    double sinI = std::sin(tle_.inclination);
    double cosI = std::cos(tle_.inclination);
    double sinNode = std::sin(xnode);
    double cosNode = std::cos(xnode);

    // TEME unit orientation vectors
    Vector3D P{
        cosNode * cosu - sinNode * sinu * cosI,
        sinNode * cosu + cosNode * sinu * cosI,
        sinu * sinI
    };

    Vector3D Q{
        -cosNode * sinu - sinNode * cosu * cosI,
        -sinNode * sinu + cosNode * cosu * cosI,
        cosu * sinI
    };

    // State vectors converted to km and km/s
    StateVector state;
    state.position = {
        mr * P.x * consts_.radiusEarth,
        mr * P.y * consts_.radiusEarth,
        mr * P.z * consts_.radiusEarth
    };

    double vFactor = (consts_.radiusEarth * consts_.xke) / 60.0;
    state.velocity = {
        (rDot * P.x + rfDot * Q.x) * vFactor,
        (rDot * P.y + rfDot * Q.y) * vFactor,
        (rDot * P.z + rfDot * Q.z) * vFactor
    };

    return state;
}

} // End AstroStuff namespace