#include "PhaseArrayBeamFormer"


namespace PAB{

PhasedArrayBeamformer::PhasedArrayBeamformer(const PhasedArrayConfig& config) : config_(config) {
    wavelength_ = SPEED_OF_LIGHT / config_.carrierFrequencyHz;
    spacing_ = (config_.elementSpacingMeters > 0.0) ? config_.elementSpacingMeters : (0.5 * wavelength_);
    waveNumber_ = (2.0 * std::numbers::pi) / wavelength_;
}

double PhasedArrayBeamformer::computeScanAngle(double azRad, double elRad) const {
    // Project line-of-sight vector onto array baseline. 
    // This is done because the delay can only be applied to antenna X relative to antenna Y, in the same line (1D).
    // Antenna does not understand 3D line-of-sight vector of the satellite.
    double deltaAz = azRad - config_.arrayAzimuthRad;
    double sinTheta = std::cos(elRad) * std::cos(deltaAz);
    
    // Clamp to valid arcsin domain [-1.0, 1.0]
    sinTheta = std::clamp(sinTheta, -1.0, 1.0);
    return std::asin(sinTheta);
}

BeamWeights PhasedArrayBeamformer::computeWeightsFromScanAngle(double scanAngleRad) const {
    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    
    BeamWeights result;
    result.scanAngleRad = scanAngleRad;
    result.phaseDelaysRad.resize(config_.numElements);
    result.weights.resize(config_.numElements);

    double normFactor = 1.0 / std::sqrt(static_cast<double>(config_.numElements));
    double spatialFreq = waveNumber_ * spacing_ * std::sin(scanAngleRad);

    for (int n = 0; n < config_.numElements; ++n) {
        // Continuous progressive phase for steering
        double phase = -static_cast<double>(n) * spatialFreq;
        
        // Wrap commanded phase delay to [0, 2pi)
        double wrappedPhase = std::fmod(phase, TWO_PI);
        if (wrappedPhase < 0.0) wrappedPhase += TWO_PI;
        
        result.phaseDelaysRad[n] = wrappedPhase;
        
        result.weights[n] = std::complex<double>(
            normFactor * std::cos(phase),
            normFactor * std::sin(phase)
        );
    }

    return result;
}
}
