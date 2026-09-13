#include "CoordinateTransforms.hpp"
#include <numbers>

namespace AstroStuff{

double CoordinateTransforms::computeGMST(double jd) {
    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    constexpr double SECONDS_PER_DAY = 86400.0;

    double t = (jd - 2451545.0) / 36525.0;
    
    // GMST in seconds
    double gmstSec = 24110.54841 + 8640184.812866 * t + 0.093104 * t * t - 6.2e-6 * t * t * t;
    
    // Convert to radians (cuz C++ math libraries only accept angle in radians) 
    // and  we normalize to [0, 2pi) to preserve precision bits and reduce the storage taken up. Basically to reduce the magnitude.
    double gmstRad = std::fmod(gmstSec * (TWO_PI / SECONDS_PER_DAY), TWO_PI);
    if (gmstRad < 0.0) {
        gmstRad += TWO_PI;
    }
    return gmstRad;
}

// TEME (Which is a non-rotating frame) is converted to ECEF (rotating frame) position, by rotating about the shared Z-axis by GMST.
// (Both TEME and ECEF have earth-centered origin, and share the same Z axis)
StateVector CoordinateTransforms::temeToEcef(const StateVector& teme, double gmst) {
    double cosG = std::cos(gmst);
    double sinG = std::sin(gmst);

    StateVector ecef;
    // rotation on position
    ecef.position.x =  cosG * teme.position.x + sinG * teme.position.y;
    ecef.position.y = -sinG * teme.position.x + cosG * teme.position.y;
    ecef.position.z = teme.position.z;

    // Rotation on velocity with Coriolis term (omega_E x r_ecef)
    double vxRot =  cosG * teme.velocity.x + sinG * teme.velocity.y;
    double vyRot = -sinG * teme.velocity.x + cosG * teme.velocity.y;
    double vzRot = teme.velocity.z;

    ecef.velocity.x = vxRot + EARTH_OMEGA * ecef.position.y;
    ecef.velocity.y = vyRot - EARTH_OMEGA * ecef.position.x;
    ecef.velocity.z = vzRot;

    return ecef;
}

Vector3D CoordinateTransforms::geodeticToEcef(const GeodeticLocation& loc) {
    double sinLat = std::sin(loc.latitude);
    double cosLat = std::cos(loc.latitude);
    double sinLon = std::sin(loc.longitude);
    double cosLon = std::cos(loc.longitude);

    double n = WGS84_A / std::sqrt(1.0 - WGS84_E2 * sinLat * sinLat);

    return {
        (n + loc.altitude) * cosLat * cosLon,
        (n + loc.altitude) * cosLat * sinLon,
        (n * (1.0 - WGS84_E2) + loc.altitude) * sinLat
    };
}

Vector3D CoordinateTransforms::ecefToEnu(const Vector3D& dEcef, double lat, double lon) {
    double sinLat = std::sin(lat);
    double cosLat = std::cos(lat);
    double sinLon = std::sin(lon);
    double cosLon = std::cos(lon);

    return {
        -sinLon * dEcef.x + cosLon * dEcef.y,
        -sinLat * cosLon * dEcef.x - sinLat * sinLon * dEcef.y + cosLat * dEcef.z,
         cosLat * cosLon * dEcef.x + cosLat * sinLon * dEcef.y + sinLat * dEcef.z
    };
}

Vector3D CoordinateTransforms::ecefVelToEnu(const Vector3D& dVel, double lat, double lon) {
    // we need the exact same matrix multiplication here as well, so lets just re-use ecefToEnu
    // cuz DRY - DO NOT REPEAT YOURSELF
    return ecefToEnu(dVel, lat, lon);
}

AERCoordinates CoordinateTransforms::computeLookAngles(
    const StateVector& satTeme,
    double jd,
    const GeodeticLocation& gs) 
{
    constexpr double TWO_PI = 2.0 * std::numbers::pi;

    double gmst = computeGMST(jd);
    StateVector satEcef = temeToEcef(satTeme, gmst);
    Vector3D gsEcef = geodeticToEcef(gs);

    Vector3D dPosEcef{
        satEcef.position.x - gsEcef.x,
        satEcef.position.y - gsEcef.y,
        satEcef.position.z - gsEcef.z
    };

    Vector3D enuPos = ecefToEnu(dPosEcef, gs.latitude, gs.longitude);
    Vector3D enuVel = ecefVelToEnu(satEcef.velocity, gs.latitude, gs.longitude);

    double range = enuPos.norm();
    double elevation = std::asin(enuPos.z / range);
    double azimuth = std::atan2(enuPos.x, enuPos.y);
    if (azimuth < 0.0) {
        azimuth += TWO_PI;
    }

    double rangeRate = (enuPos.x * enuVel.x + enuPos.y * enuVel.y + enuPos.z * enuVel.z) / range;

    return { azimuth, elevation, range, rangeRate };
}

}