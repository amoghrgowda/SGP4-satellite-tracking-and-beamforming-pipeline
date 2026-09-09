#pragma once
#include <string>
#include <cmath>

namespace AstroStuff{

struct Vector3D{
    double x {0.0};
    double y {0.0};
    double z {0.0};

    double norm() const{
       return std::sqrt(x*x + y*y + z*z); // we need this euclidean length (scalar length or wtv). 
    }                               // This magnitude is helpful to calculate altitude of satellite and the orbital speed   
};

struct StateVector{
    // both in TEME (True Equator, Mean Equinox)
    Vector3D position;
    Vector3D velocity;
};

struct TLE{
    // Two Line Element components (provided by the overlords at US space force):
    std::string name;
    int satNumber;
    char classification;
    int idLaunchYear;
    int idLaunchNumber;
    std::string idLaunchPiece;
    
    // Metadata
    int epochYear;
    double epochDay;
    double epochJulianDate;
    
    //decay terms
    double bstar;          // Drag term (1/Earth radii)
    double meanMotionDot;   // First time derivative of mean motion
    double meanMotionDDot;  // Second time derivative of mean motion
    
    // orbital parameters
    double inclination;    
    double raan;           
    double eccentricity; 
    double argPerigee;     
    double meanAnomaly;   
    double meanMotion;     
    
    //Mission counter
    int revNumberAtEpoch;
};

struct GravitationalConstants{
    // Will follow the WGS-72 standard for SGP4
    double mu;         // standard gravitational parameter. 

    // Don't ask what the rest do. Let the Math geeks handle it.
    double radiusEarth;
    double xke;        
    double j2;         
    double j3;         
    double j4;         
    double j3oj2;
};

inline GravitationalConstants const WGS72 = {
    398600.8,                   
    6378.135,                    
    0.0743669161331734132,     
    0.001082616,               
    -0.00000253881,             
    -0.00000165597,               
    -0.00000253881 / 0.001082616 
};

} // end namespace 