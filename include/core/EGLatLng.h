// Copyright 2005 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

// Author: ericv@google.com (Eric Veach)

#ifndef S2_EGLatLng_H_
#define S2_EGLatLng_H_

#include <cmath>
#include <iosfwd>
#include <ostream>
#include <string>

#include "base/integral_types.h"
#include "_fp_contract_off.h"
#include "core/r2.h"
#include "core/EG1DAngle.h"
#include "util/math/vector.h"

// This class represents a point on the unit sphere as a pair
// of latitude-longitude coordinates.  Like the rest of the "geometry"
// package, the intent is to represent spherical geometry as a mathematical
// abstraction, so functions that are specifically related to the Earth's
// geometry (e.g. easting/northing conversions) should be put elsewhere.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class EGLatLng {
public:
    // Constructor.  The latitude and longitude are allowed to be outside
    // the is_valid() range.  However, note that most methods that accept
    // EGLatLngs expect them to be normalized (see Normalized() below).
    EGLatLng(EG1DAngle lat, EG1DAngle lng);

    // The default constructor sets the latitude and longitude to zero.  This is
    // mainly useful when declaring arrays, STL containers, etc.
    EGLatLng();

    // Convert a direction vector (not necessarily unit length) to an EGLatLng.
    explicit EGLatLng(const EGPoint& p);

    // Returns an EGLatLng for which is_valid() will return false.
    static EGLatLng Invalid();

    // Convenience functions -- shorter than calling EG1DAngle::Radians(), etc.
    static EGLatLng FromRadians(double lat_radians, double lng_radians);
    static EGLatLng FromDegrees(double lat_degrees, double lng_degrees);
    static EGLatLng FromE5(int32 lat_e5, int32 lng_e5);
    static EGLatLng FromE6(int32 lat_e6, int32 lng_e6);
    static EGLatLng FromE7(int32 lat_e7, int32 lng_e7);

    // Convenience functions -- to use when args have been fixed32s in protos.
    //
    // The arguments are static_cast into int32, so very large unsigned values
    // are treated as negative numbers.
    static EGLatLng FromUnsignedE6(uint32 lat_e6, uint32 lng_e6);
    static EGLatLng FromUnsignedE7(uint32 lat_e7, uint32 lng_e7);

    // Methods to compute the latitude and longitude of a point separately.
    static EG1DAngle Latitude(const EGPoint& p);
    static EG1DAngle Longitude(const EGPoint& p);

    // Accessor methods.
    EG1DAngle lat() const { return EG1DAngle::Radians(coords_[0]); }
    EG1DAngle lng() const { return EG1DAngle::Radians(coords_[1]); }
    const R2Point& coords() const { return coords_; }

    // Return true if the latitude is between -90 and 90 degrees inclusive
    // and the longitude is between -180 and 180 degrees inclusive.
    bool is_valid() const;

    // Clamps the latitude to the range [-90, 90] degrees, and adds or subtracts
    // a multiple of 360 degrees to the longitude if necessary to reduce it to
    // the range [-180, 180].
    EGLatLng Normalized() const;

    // Converts a normalized EGLatLng to the equivalent unit-length vector.
    // The maximum error in the result is 1.5 * DBL_EPSILON.  (This does not
    // include the error of converting degrees, E5, E6, or E7 to radians.)
    //
    // Can be used just like an EGPoint constructor.  For example:
    //   S2Cap cap;
    //   cap.AddPoint(EGPoint(latlng));
    explicit operator EGPoint() const;

    // Converts to an EGPoint (equivalent to the operator above).
    EGPoint ToPoint() const;

    // Returns the distance (measured along the surface of the sphere) to the
    // given EGLatLng, implemented using the Haversine formula.  This is
    // equivalent to
    //
    //   EG1DAngle(ToPoint(), o.ToPoint())
    //
    // except that this function is slightly faster, and is also somewhat less
    // accurate for distances approaching 180 degrees (see EG1DAngle.h for
    // details).  Both EGLatLngs must be normalized.
    EG1DAngle GetDistance(const EGLatLng& o) const;

    // Simple arithmetic operations for manipulating latitude-longitude pairs.
    // The results are not normalized (see Normalized()).
    friend EGLatLng operator+(const EGLatLng& a, const EGLatLng& b);
    friend EGLatLng operator-(const EGLatLng& a, const EGLatLng& b);
    friend EGLatLng operator*(double m, const EGLatLng& a);
    friend EGLatLng operator*(const EGLatLng& a, double m);

    bool operator==(const EGLatLng& o) const { return coords_ == o.coords_; }
    bool operator!=(const EGLatLng& o) const { return coords_ != o.coords_; }
    bool operator<(const EGLatLng& o) const { return coords_ < o.coords_; }
    bool operator>(const EGLatLng& o) const { return coords_ > o.coords_; }
    bool operator<=(const EGLatLng& o) const { return coords_ <= o.coords_; }
    bool operator>=(const EGLatLng& o) const { return coords_ >= o.coords_; }

    bool ApproxEquals(const EGLatLng& o,
                      EG1DAngle max_error = EG1DAngle::Radians(1e-15)) const {
        return coords_.aequal(o.coords_, max_error.radians());
    }

    // Exports the latitude and longitude in degrees, separated by a comma.
    // e.g. "94.518000,150.300000"
    string ToStringInDegrees() const;
    void ToStringInDegrees(string* s) const;

private:
    // Internal constructor.
    explicit EGLatLng(const R2Point& coords) : coords_(coords) {}

    // This is internal to avoid ambiguity about which units are expected.
    EGLatLng(double lat_radians, double lng_radians)
            : coords_(lat_radians, lng_radians) {}

    R2Point coords_;
};

// Hasher for EGLatLng.
// Example use: std::unordered_map<EGLatLng, int, EGLatLngHash>.
struct EGLatLngHash {
    size_t operator()(const EGLatLng& lat_lng) const {
        return GoodFastHash<R2Point>()(lat_lng.coords());
    }
};

//////////////////   Implementation details follow   ////////////////////


inline EGLatLng::EGLatLng(EG1DAngle lat, EG1DAngle lng)
        : coords_(lat.radians(), lng.radians()) {}

inline EGLatLng::EGLatLng() : coords_(0, 0) {}

inline EGLatLng EGLatLng::FromRadians(double lat_radians, double lng_radians) {
    return EGLatLng(lat_radians, lng_radians);
}

inline EGLatLng EGLatLng::FromDegrees(double lat_degrees, double lng_degrees) {
    return EGLatLng(EG1DAngle::Degrees(lat_degrees), EG1DAngle::Degrees(lng_degrees));
}

inline EGLatLng EGLatLng::FromE5(int32 lat_e5, int32 lng_e5) {
    return EGLatLng(EG1DAngle::E5(lat_e5), EG1DAngle::E5(lng_e5));
}

inline EGLatLng EGLatLng::FromE6(int32 lat_e6, int32 lng_e6) {
    return EGLatLng(EG1DAngle::E6(lat_e6), EG1DAngle::E6(lng_e6));
}

inline EGLatLng EGLatLng::FromE7(int32 lat_e7, int32 lng_e7) {
    return EGLatLng(EG1DAngle::E7(lat_e7), EG1DAngle::E7(lng_e7));
}

inline EGLatLng EGLatLng::FromUnsignedE6(uint32 lat_e6, uint32 lng_e6) {
    return EGLatLng(EG1DAngle::UnsignedE6(lat_e6), EG1DAngle::UnsignedE6(lng_e6));
}

inline EGLatLng EGLatLng::FromUnsignedE7(uint32 lat_e7, uint32 lng_e7) {
    return EGLatLng(EG1DAngle::UnsignedE7(lat_e7), EG1DAngle::UnsignedE7(lng_e7));
}

inline EGLatLng EGLatLng::Invalid() {
    // These coordinates are outside the bounds allowed by is_valid().
    return EGLatLng(M_PI, 2 * M_PI);
}

inline EG1DAngle EGLatLng::Latitude(const EGPoint& p) {
    // We use atan2 rather than asin because the input vector is not necessarily
    // unit length, and atan2 is much more accurate than asin near the poles.
    return EG1DAngle::Radians(atan2(p[2], sqrt(p[0]*p[0] + p[1]*p[1])));
}

inline EG1DAngle EGLatLng::Longitude(const EGPoint& p) {
    // Note that atan2(0, 0) is defined to be zero.
    return EG1DAngle::Radians(atan2(p[1], p[0]));
}

inline bool EGLatLng::is_valid() const {
    return (std::fabs(lat().radians()) <= M_PI_2 &&
            std::fabs(lng().radians()) <= M_PI);
}

inline EGLatLng::operator EGPoint() const {
    return ToPoint();
}

inline EGLatLng operator+(const EGLatLng& a, const EGLatLng& b) {
    return EGLatLng(a.coords_ + b.coords_);
}

inline EGLatLng operator-(const EGLatLng& a, const EGLatLng& b) {
    return EGLatLng(a.coords_ - b.coords_);
}

inline EGLatLng operator*(double m, const EGLatLng& a) {
    return EGLatLng(m * a.coords_);
}

inline EGLatLng operator*(const EGLatLng& a, double m) {
    return EGLatLng(m * a.coords_);
}

std::ostream& operator<<(std::ostream& os, const EGLatLng& ll);

#endif  // S2_EGLatLng_H_
