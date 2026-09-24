#pragma once
#include <vector>
#include <complex>
#include "CoordinateTransforms.hpp"

namespace PAB{

struct PhasedArrayConfig{
    int numElements{16};    //number of antenna elements
    double carriedFrequency{2.4e9}; // S-band 2.4 GHz
    double elementSpacingMeters{0.0};
    double arrayAzimuthRad{0.0}; // Baseline orientation relative to North ( 0 means North-South orientation)
};

struct BeamWeights{
    double scanAngleRad{0.0};   // current direction that the antenna is pointing to
    std::vector<double> phaseDelays; // in radians, and per-element
    std::vector<std::complex<double>> weights; // phase shifting is simple to do in complex numbers than to use trigonometry.
};

struct ArrayPatternSample{  // represents a discrete data point on a radiation pattern plot. 
// When testing an antenna array, we sweep an angle across the sky to see how much power the beam radiates in every angle.
    double angleDeg{0.0}; // Observation scan angle -90 degree to +90 degree (left to right of the broadside, aka, perpendicular to array)
    double gainLinear{0.0}; // If all 16 arrays line up constructively, gainLinear = 16. 0 if they are destructive line up.
    double gainDbi{0.0}; // antenna pattern in decibels (dB) instead of linear to make changes more visible on graph plot.
};

class PhasedArrayBeamFormer{
public:
    static constexpr double SPEED_OF_LIGHT = 299792458.0; // in m/s 
    explicit PhasedArrayBeamFormer(const PhasedArrayConfig& config);
    double computeScanAngle(double azimuthRad, double elevationRad) const; //project topocentric to AzEl look angles into ULA scan angles.
    BeamWeights computeSteeringWeights(double azimuthRad, double elevationRad) const; // computer per-element phase shift
    BeamWeights computeWeightsFromScanAngle(double scanAngleRad) const; // Computes steering weights directly from 1D scan angle
    
    //Synthesizes the far-field Array Factor radiation pattern
    std::vector<ArrayPatternSample> synthesizePattern(const BeamWeights& beam, int numPoints = 361) const;
    
    // Computes key RF beam metrics (HPBW, First Sidelobe Level)
    double estimateHPBW(double scanAngleRad) const;

    double getWavelength() const { 
        return wavelength_; 
    }

    double getElementSpacing() const { 
        return spacing_; 
    }

private:
    PhasedArrayConfig config_;
    double wavelength_;
    double spacing_;
    double waveNumber_;

};

} // end namespace PAB