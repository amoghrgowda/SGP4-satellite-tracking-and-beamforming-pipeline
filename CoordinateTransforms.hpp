#pragma once
#include "SGP4Types.hpp"
#include <numbers>
#include <cmath>

namespace AstroStuff{
    //Geodetic location is NOT the same as astronomical coords (meaning, it does NOT wobble and is fixed for convenience).
    
// Geodetic latitude = The angle between perpendicular of a location (reference ellipsoid WGS-84, and NOT perpendicular of the terrain or gravity or something else) 
// -and the fixed (non-wobbling) Equatorial plane. 
struct GeodeticLocation {
    double latitude;
    double longitude;
    double altitude;
};

struct AERCoordinates {
    double azimuth;
    double elevation;
    double range;     // Kilometers (slant range)
    double rangeRate; // representing how quickly the straight-line distance (slant range) is changing per second.
};

class CoordinateTransforms {
public:
    //WGS-84 constants
    static constexpr double WGS84_A = 6378.137;
    static constexpr double WGS84_F = 1.0 / 298.257223563;
    static constexpr double WGS84_E2 = 2.0 * WGS84_F - WGS84_F * WGS84_F; // Square of Eccentricity (E)
    static constexpr double EARTH_OMEGA = 7.292115146706979e-5; // rad/s

    // Computes Greenwich Mean Sidereal Time (in radians) for a given Julian Date
    static double computeGMST(double julianDate);

    // Transforms TEME State Vector to ECEF State Vector
    static StateVector temeToEcef(const StateVector& temeState, double gmstRad);

    // Transforms Geodetic coordinates (of our ground station) to ECEF position vector
    static Vector3D geodeticToEcef(const GeodeticLocation& loc);

    // Transforms relative ECEF (= ECEF state - position) vector to local Topocentric/ENU vector
    static Vector3D ecefToEnu(const Vector3D& deltaEcef, double latRad, double lonRad);

    // Transforms relative ECEF velocity to ENU velocity vector (for driving antenna motors)
    static Vector3D ecefVelToEnu(const Vector3D& deltaVelEcef, double latRad, double lonRad);

    // Computes AER -> Azimuth, Elevation, and Range
    static AERCoordinates computeLookAngles(
        const StateVector& satTemeState,
        double currentJulianDate,
        const GeodeticLocation& groundStation
    );
};
}