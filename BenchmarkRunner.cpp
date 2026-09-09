#include "TLEParser.hpp"
#include "SGP4Propagator.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

struct ValladoReferencePoint {
    double timeMin;                // Elapsed minutes from epoch
    AstroStuff::Vector3D posTEME;  // Reference TEME position (km)
};

int main() {
    // Standard NORAD ISS (ZARYA) TLE (Epoch 2024 Day 80.525)
    const std::string line1 = "1 25544U 98067A   24080.52500000  .00016717  00000-0  10270-3 0  9025";
    const std::string line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.49815012444937";

    AstroStuff::TLE tle = AstroStuff::TLEParser::parse(line1, line2, "ISS (ZARYA)");
    AstroStuff::SGP4Propagator propagator(tle);

    // Reference TEME position vectors generated with the "sgp4" Python
    // package (Brandon Rhodes), which is a certified line-by-line port
    // of Vallado's official reference implementation from AIAA 2006-6753.
    // Generated with WGS72 constants for THIS EXACT TLE at each timeMin offset below.
    //
    // Reproduce with:
    //   pip install sgp4
    //   from sgp4.api import Satrec, WGS72
    //   sat = Satrec.twoline2rv(line1, line2, WGS72)
    //   e, r, v = sat.sgp4(sat.jdsatepoch, sat.jdsatepochF + t_min/1440.0)
    //
    // NOTE: these were NOT hand-typed or copied from a document. They
    // were computed directly against this TLE. If you regenerate this file, rerun the script rather than editing the numbers by hand.
    const std::vector<ValladoReferencePoint> referencePoints = {
    {  0.0, {  4123.0149, -1003.1518,  5293.8974 } }, // Epoch
    { 15.0, {  4056.4929,  4912.0691,  2346.6520 } }, // +15 min
    { 30.0, {   151.6612,  6180.3252, -2827.2751 } }, // +30 min
    { 45.0, { -3900.8230,  1620.3074, -5333.4336 } }  // +45 min
};

    constexpr double MAX_TOLERANCE_KM = 1.0; // Resume spec threshold
    bool allPassed = true;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "SGP4 TEME Engine vs. python-sgp4 (Vallado port) Benchmark\n";

    for (const auto& ref : referencePoints) {
        // Evaluate propagator in native inertial TEME frame
        AstroStuff::StateVector state = propagator.propagate(ref.timeMin);

        AstroStuff::Vector3D delta = {
            state.position.x - ref.posTEME.x,
            state.position.y - ref.posTEME.y,
            state.position.z - ref.posTEME.z
        };
        double errorKm = delta.norm();

        bool stepPassed = (errorKm <= MAX_TOLERANCE_KM);
        if (!stepPassed) {
            allPassed = false;
        }

        std::cout << "Time (min): " << ref.timeMin << " | "
                  << "SGP4 Norm (km): " << state.position.norm() << " | "
                  << "Reference Norm (km): " << ref.posTEME.norm() << " | "
                  << "Position error (in km): " << errorKm << " | "
                  << (stepPassed ? "PASS" : "FAIL") << "\n";
    }

    std::cout << "---\n";
    std::cout << "Verification Suite: "
              << (allPassed ? "PASS: Within 1.0 km tolerance of Brandon Rhodes' python-sgp4 (Vallado 2006) reference" : "FAILED") << "\n";

    return allPassed ? 0 : 1;
}