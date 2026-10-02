#include "TLEParser.hpp"
#include "SGP4Propagator.hpp"
#include "CoordinateTransforms.hpp"
#include "PhasedArrayBeamFormer.hpp"
#include <iostream>
#include <iomanip>

int main() {
    constexpr double RAD_TO_DEG = 180.0 / std::numbers::pi;

    // Ground Station: Sydney (-33.8688 S, 151.2093 E, 50m)
    AstroStuff::GeodeticLocation gs{
        -33.8688 * (std::numbers::pi / 180.0),
        151.2093 * (std::numbers::pi / 180.0),
        0.050
    };

    // S-Band 16-Element Phased Array Setup
    PAB::PhasedArrayConfig arrayConfig{
        .numElements = 16,
        .carriedFrequency = 2.45e9, // 2.45 GHz S-Band
        .elementSpacingMeters = 0.0612, // lambda / 2 spacing
        .arrayAzimuthRad = 0.0 // Aligned North-South
    };
    PAB::PhasedArrayBeamFormer beamformer(arrayConfig);

    std::string l1 = "1 25544U 98067A   24080.52500000  .00016717  00000-0  10270-3 0  9025";
    std::string l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.49815012444937";

    AstroStuff::TLE tle = AstroStuff::TLEParser::parse(l1, l2, "ISS (ZARYA)");
    AstroStuff::SGP4Propagator propagator(tle);

    std::cout << "Real-Time Phased Array Beamsteering & Phase Delay Engine:\n";
    std::cout << "Carrier: " << arrayConfig.carriedFrequency / 1e9 << " GHz | Elements: " << arrayConfig.numElements
              << " | Spacing: " << arrayConfig.elementSpacingMeters * 100.0 << " cm (lambda/2)\n";
    std::cout << "----\n";
    std::cout << "Time(m) | Az(deg) | El(deg) | ScanAngle | HPBW(deg) | Elem 0 Ph | Elem 1 Ph | Elem 4 Ph | Elem 8 Ph | Elem 15 Ph\n";
    std::cout << "----\n";

    for (double t = 454.0; t <= 461.0; t += 1.0) {
        double currentJD = tle.epochJulianDate + (t / 1440.0);
        AstroStuff::StateVector stateTeme = propagator.propagate(t);
        
        AstroStuff::AERCoordinates aer = AstroStuff::CoordinateTransforms::computeLookAngles(stateTeme, currentJD, gs);
        PAB::BeamWeights beam = beamformer.computeSteeringWeights(aer.azimuth, aer.elevation);
        
        double scanDeg = beam.scanAngleRad * RAD_TO_DEG;
        double hpbwDeg = beamformer.estimateHPBW(beam.scanAngleRad) * RAD_TO_DEG;

        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(7) << t << " | "
                  << std::setw(7) << aer.azimuth * RAD_TO_DEG << " | "
                  << std::setw(7) << aer.elevation * RAD_TO_DEG << " | "
                  << std::setw(8) << scanDeg << "deg | "
                  << std::setw(8) << hpbwDeg << "deg | "
                  << std::setw(8) << beam.phaseDelays[0] * RAD_TO_DEG << "deg | "
                  << std::setw(8) << beam.phaseDelays[1] * RAD_TO_DEG << "deg | "
                  << std::setw(8) << beam.phaseDelays[4] * RAD_TO_DEG << "deg | "
                  << std::setw(8) << beam.phaseDelays[8] * RAD_TO_DEG << "deg | "
                  << std::setw(8) << beam.phaseDelays[15] * RAD_TO_DEG << "deg\n";
    }
    std::cout << "EoF\n";

    return 0;
}