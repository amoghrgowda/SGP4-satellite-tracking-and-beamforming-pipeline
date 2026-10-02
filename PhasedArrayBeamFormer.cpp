#include "PhasedArrayBeamFormer.hpp"
#include <cmath>
#include <algorithm>

namespace PAB{

PhasedArrayBeamFormer::PhasedArrayBeamFormer(const PhasedArrayConfig& config) : config_(config) {
    wavelength_ = SPEED_OF_LIGHT / config_.carriedFrequency;
    spacing_ = (config_.elementSpacingMeters > 0.0) ? config_.elementSpacingMeters : (0.5 * wavelength_);
    waveNumber_ = (2.0 * std::numbers::pi) / wavelength_;
}

double PhasedArrayBeamFormer::computeScanAngle(double azRad, double elRad) const {
    // Project line-of-sight vector onto array baseline. 
    // This is done because the delay can only be applied to antenna X relative to antenna Y, in the same line (1D).
    // Antenna does not understand 3D line-of-sight vector of the satellite.
    double deltaAz = azRad - config_.arrayAzimuthRad;
    double sinTheta = std::cos(elRad) * std::cos(deltaAz);
    
    // Clamp to valid arcsin domain [-1.0, 1.0]
    sinTheta = std::clamp(sinTheta, -1.0, 1.0);
    return std::asin(sinTheta);
}

BeamWeights PhasedArrayBeamFormer::computeWeightsFromScanAngle(double scanAngleRad) const {
    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    
    BeamWeights result;
    result.scanAngleRad = scanAngleRad;
    result.phaseDelays.resize(config_.numElements);
    result.weights.resize(config_.numElements);

    double normFactor = 1.0 / std::sqrt(static_cast<double>(config_.numElements));
    double spatialFreq = waveNumber_ * spacing_ * std::sin(scanAngleRad);

    for (int n = 0; n < config_.numElements; ++n) {
        // Continuous progressive phase for steering
        double phase = -static_cast<double>(n) * spatialFreq;
        
        // Wrap commanded phase delay to [0, 2pi)
        double wrappedPhase = std::fmod(phase, TWO_PI);
        if (wrappedPhase < 0.0) wrappedPhase += TWO_PI;
        
        result.phaseDelays[n] = wrappedPhase;
        
        result.weights[n] = std::complex<double>(
            normFactor * std::cos(phase),
            normFactor * std::sin(phase)
        );
    }

    return result;
}

BeamWeights PhasedArrayBeamFormer::computeSteeringWeights(double azRad, double elRad)const {
    double scanAngle = computeScanAngle(azRad, elRad);
    return computeWeightsFromScanAngle(scanAngle);
}

// Now we simulate/mock the far field radiation pattern of the antenna array
// 
std::vector<ArrayPatternSample> PhasedArrayBeamFormer::synthesizePattern(const BeamWeights& beam, int numPoints) const{

    constexpr double DEG_TO_RAD = std::numbers::pi / 180.0;
    std::vector<ArrayPatternSample> pattern(numPoints);

    double maxGain = 0.0;
    std::vector<double> linearGains(numPoints);

    // sweeping the observation angle to find the total combined signal strength in every observable direction
    // at the end of this for loop, we will have 'max gain'/point at which beam peaks, and datapoints to plot on the graph if, we want 
    for (int i = 0; i < numPoints; ++i) {
        double angleDeg = -90.0 + (180.0 * i) / (numPoints - 1);
        double thetaRad = angleDeg * DEG_TO_RAD;

        std::complex<double> af{0.0, 0.0};
        double spatialFreq = waveNumber_ * spacing_ * std::sin(thetaRad);

        for (int n = 0; n < config_.numElements; ++n) {
            double elemPhase = static_cast<double>(n) * spatialFreq;
            std::complex<double> steeringElem(std::cos(elemPhase), std::sin(elemPhase));
            
            af += std::conj(beam.weights[n]) * steeringElem;
        }

        double gain = std::abs(af);
        linearGains[i] = gain;
        if (gain > maxGain) {
            maxGain = gain;
        }
        pattern[i].angleDeg = angleDeg;
        pattern[i].gainLinear = gain;
    }

    // Normalize to peak directivity (dBi)-Decibels relative to Isotropic: 
    // Compares antenna gain to a hypothetical ideal antenna that radiates equally in all 3D directions (0 dBi). 
    // A 16-element array typically has a peak directivity of around +12 dBi.
    for (int i = 0; i < numPoints; ++i) {
        if (maxGain > 1e-12) {
            double normPower = linearGains[i] / maxGain;
            pattern[i].gainDbi = 20.0 * std::log10(std::max(normPower, 1e-6));
        } else {
            pattern[i].gainDbi = -120.0;
        }
    }

    return pattern;
}

double PhasedArrayBeamFormer::estimateHPBW(double scanAngleRad) const {
    double cosTheta = std::cos(scanAngleRad);
    if (std::abs(cosTheta) < 1e-4) return std::numbers::pi;
    return (0.886 * wavelength_) / (config_.numElements * spacing_ * cosTheta);
}

}
