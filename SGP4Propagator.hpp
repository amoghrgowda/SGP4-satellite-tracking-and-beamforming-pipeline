#pragma once
#include "SGP4Types.hpp"
#include <numbers>

// Intended to be a core calculation engine of SGP4 algo. 
// Input is some delta-time and output is a StateVector (3D pos & velocity).

namespace AstroStuff{

class SGP4Propagator{
public:
    explicit SGP4Propagator(const TLE& tle, const GravitationalConstants& constants = WGS72);

    StateVector propagate(double deltaMinutes) const;

private:
    TLE tle_;
    GravitationalConstants consts_;

    //We need brouwer elements/coefficients to be initialized 
    // Elements considered from STR-3 report, for Low Earth Orbit (LEO) ONLY!!! - keeping it simple.

    // Thanks to Daniel Warner on github (dnwrnr) for these coefficient lists in C++ 

    //Brouwer constants:
    double a0_double_prime_{0.0};
    double n0_double_prime_{0.0};

    // Atmospheric layers
    double s4_{0.0};
    double qoms24_{0.0};
    double perigee_{0.0};

    //Drag geometry
    double tsi_{0.0};
    double eta_{0.0};

    //Drag cross-terms
    double c1_{0.0};
    double c2_{0.0};
    double c3_{0.0};
    double c4_{0.0};
    double c5_{0.0};

    //Decay polynomials
    double d2_{0.0};
    double d3_{0.0};
    double d4_{0.0};

    //secular rates
    double mDot_{0.0};
    double omegaDot_{0.0};
    double omegamDot_{0.0};

    // skipped the inclination powers and periodic helpers here, which are stated in Daniel Warner's port.
    // Will calculate them as and when required.

    void initialize();

    // The satellite does NOT move in a perfect circular orbit. Therefore, the speed varies at different points in the ellipse. 
    // We need kepler's equation to find the position of the satellite at any given time.
    static double solveKepler(double M, double e, double tol = 1e-12, int maxIter = 50);

};
}   // end namespace AstroStuff