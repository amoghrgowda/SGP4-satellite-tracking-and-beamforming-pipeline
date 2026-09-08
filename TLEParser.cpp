#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <numbers>
#include "SGP4Types.hpp"

namespace AstroStuff{

class TLEParser{
private:
    static constexpr double DEG_TO_RAD = std::numbers::pi / 180.0;
    static constexpr double TWO_PI = std::numbers:pi * 2.0;
    static constexpr double MINUTES_PER_DAY = 1440.0;

    // Parses TLE data format into SGP4 expected decimal float format.
    // TLE format (10270-3) is to be converted to SGP4 expectated format (0.10270 * 10^-3)
    static double parseDecimalAssumedScientific(const std::string& str) {
        std::string cleaned;
        for (char c : str) {
            if (c != ' ')
                cleaned += c;
        }
        if (cleaned.empty() || cleaned == "0") return 0.0;

        char sign = '+';
        size_t startIdx = 0;
        if (cleaned[0] == '-' || cleaned[0] == '+') {
            sign = cleaned[0];
            startIdx = 1;
        }

        size_t expPos = cleaned.find_last_of("+-", cleaned.length() - 1);
        if (expPos == std::string::npos || expPos < startIdx) {
            return std::stod(cleaned);
        }

        std::string mantissa = cleaned.substr(startIdx, expPos - startIdx);
        std::string exponent = cleaned.substr(expPos);

        double val = std::stod("0." + mantissa) * std::pow(10.0, std::stod(exponent));
        return (sign == '-') ? -val : val;
    }
    // Julian date here is fixed (not a sliding window) for 100 years (1957 to 2056) for simplicity
    // Julian date starts in 4713 BC.
    static double epochToJulianDate(int year, double dayFraction) {
        int fullYear = (year < 57) ? (2000 + year) : (1900 + year);
        
        // Computing Julian date for Jan 0.0 of fullYear (using the Jean Meeus's algorithm)
        int y = fullYear - 1;   // shifting year by 1 backwards -
        // -because Jan and Feb belong to prev year for simple calculations of date -
        // -and not worrying about leap or non leap feb days count ruining the cumulative number of days in all the other months.
        int a = y / 100;    // calc century
        int b = 2 - a + (a / 4);    // Gregorian calendar correction term
        double jdJan0 = std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * 14) + 0.0 + b - 1524.5;
        
        return jdJan0 + dayFraction;
    }

public:
    static TLE parse(const std::string& line1, const std::string& line2, const std::string& name = "SAT") {
        if (line1.length() < 68 || line2.length() < 68) {
            throw std::runtime_error("Invalid TLE line length.");
        }

        TLE tle;
        tle.name = name;
        // parsing line 1
        tle.satNumber = std::stoi(line1.substr(2, 5));
        tle.classification = line1[7];
        tle.idLaunchYear = std::stoi(line1.substr(9, 2));
        tle.idLaunchNumber = std::stoi(line1.substr(11, 3));
        tle.idLaunchPiece = line1.substr(14, 3);
        
        tle.epochYear = std::stoi(line1.substr(18, 2));
        tle.epochDay = std::stod(line1.substr(20, 12));
        tle.epochJulianDate = epochToJulianDate(tle.epochYear, tle.epochDay);
        
        tle.meanMotionDot = std::stod(line1.substr(33, 10)) * (TWO_PI / (MINUTES_PER_DAY * MINUTES_PER_DAY));
        tle.meanMotionDDot = parseDecimalAssumedScientific(line1.substr(44, 8)) * (TWO_PI / (MINUTES_PER_DAY * MINUTES_PER_DAY * MINUTES_PER_DAY));
        tle.bstar = parseDecimalAssumedScientific(line1.substr(53, 8));

        // Line 2 parsing
        tle.inclination = std::stod(line2.substr(8, 8)) * DEG_TO_RAD;
        tle.raan = std::stod(line2.substr(17, 8)) * DEG_TO_RAD;
        tle.eccentricity = std::stod("0." + line2.substr(26, 7));
        tle.argPerigee = std::stod(line2.substr(34, 8)) * DEG_TO_RAD;
        tle.meanAnomaly = std::stod(line2.substr(43, 8)) * DEG_TO_RAD;
        
        // Convert all units to suitable SGP4 internal units
        double revsPerDay = std::stod(line2.substr(52, 11));
        tle.meanMotion = revsPerDay * (TWO_PI / MINUTES_PER_DAY);
        
        tle.revNumberAtEpoch = std::stoi(line2.substr(63, 5));

        return tle;
    }
};
}