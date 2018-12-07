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

#ifndef EG_EG1DAngle_H_
#define EG_EG1DAngle_H_

#include <cmath>
#include <limits>
#include <ostream>
#include <type_traits>

#include "base/integral_types.h"
#include "_fp_contract_off.h"
#include "EGPoint.h"
#include "util/math/mathutil.h"
#include "util/math/vector.h"
#include "base/def.h"

class EGLatLng;

#ifndef SWIG
#define IFNDEF_SWIG(x) x
#else
#define IFNDEF_SWIG(x)
#endif

// This class represents a one-dimensional angle (as opposed to a
// two-dimensional solid angle).  It has methods for converting angles to
// or from radians, degrees, and the E5/E6/E7 representations (i.e. degrees
// multiplied by 1e5/1e6/1e7 and rounded to the nearest integer).
//
// The internal representation is a double-precision value in radians, so
// conversion to and from radians is exact.  Conversions between E5, E6, E7,
// and Degrees are not always exact; for example, Degrees(3.1) is different
// from E6(3100000) or E7(310000000).  However, the following properties are
// guaranteed for any integer "n", provided that "n" is in the input range of
// both functions:
//
//     Degrees(n) == E6(1000000 * n)
//     Degrees(n) == E7(10000000 * n)
//          E6(n) == E7(10 * n)
//
// The corresponding properties are *not* true for E5, so if you use E5 then
// don't test for exact equality when comparing to other formats such as
// Degrees or E7.
//
// The following conversions between degrees and radians are exact:
//
//          Degrees(180) == Radians(M_PI)
//       Degrees(45 * k) == Radians(k * M_PI / 4)  for k == 0..8
//
// These identities also hold when the arguments are scaled up or down by any
// power of 2.  Some similar identities are also true, for example,
// Degrees(60) == Radians(M_PI / 3), but be aware that this type of identity
// does not hold in general.  For example, Degrees(3) != Radians(M_PI / 60).
//
// Similarly, the conversion to radians means that Angle::Degrees(x).degrees()
// does not always equal "x".  For example,
//
//         EG1DAngle::Degrees(45 * k).degrees() == 45 * k      for k == 0..8
//   but       EG1DAngle::Degrees(60).degrees() != 60.
//
// This means that when testing for equality, you should allow for numerical
// errors (EXPECT_DOUBLE_EQ) or convert to discrete E5/E6/E7 values first.
//
// CAVEAT: All of the above properties depend on "double" being the usual
// 64-bit IEEE 754 type (which is true on almost all modern platforms).
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class EG1DAngle {
public:
    // These methods construct EG1DAngle objects from their measure in radians
    // or degrees.
    static constexpr EG1DAngle Radians(double radians);
    static constexpr EG1DAngle Degrees(double degrees);
    static constexpr EG1DAngle E5(int32 e5);
    static constexpr EG1DAngle E6(int32 e6);
    static constexpr EG1DAngle E7(int32 e7);

    // Convenience functions -- to use when args have been fixed32s in protos.
    //
    // The arguments are static_cast into int32, so very large unsigned values
    // are treated as negative numbers.
    static constexpr EG1DAngle UnsignedE6(uint32 e6);
    static constexpr EG1DAngle UnsignedE7(uint32 e7);

    // The default constructor yields a zero angle.  This is useful for STL
    // containers and class methods with output arguments.
    IFNDEF_SWIG(constexpr) EG1DAngle() : radians_(0) {}

    // Return an angle larger than any finite angle.
    static constexpr EG1DAngle Infinity();

    // A explicit shorthand for the default constructor.
    static constexpr EG1DAngle Zero();

    // Return the angle between two points, which is also equal to the distance
    // between these points on the unit sphere.  The points do not need to be
    // normalized.  This function has a maximum error of 3.25 * DBL_EPSILON (or
    // 2.5 * DBL_EPSILON for angles up to 1 radian).
    EG1DAngle(const EGPoint& x, const EGPoint& y);

    // Like the constructor above, but return the angle (i.e., distance) between
    // two EGLatLng points.  This function has about 15 digits of accuracy for
    // small distances but only about 8 digits of accuracy as the distance
    // approaches 180 degrees (i.e., nearly-antipodal points).
    EG1DAngle(const EGLatLng& x, const EGLatLng& y);

    constexpr double radians() const;
    constexpr double degrees() const;

    int32 e5() const;
    int32 e6() const;
    int32 e7() const;

    // Return the absolute value of an angle.
    EG1DAngle abs() const;

    // Comparison operators.
    friend bool operator==(EG1DAngle x, EG1DAngle y);
    friend bool operator!=(EG1DAngle x, EG1DAngle y);
    friend bool operator<(EG1DAngle x, EG1DAngle y);
    friend bool operator>(EG1DAngle x, EG1DAngle y);
    friend bool operator<=(EG1DAngle x, EG1DAngle y);
    friend bool operator>=(EG1DAngle x, EG1DAngle y);

    // Simple arithmetic operators for manipulating EG1DAngles.
    friend EG1DAngle operator-(EG1DAngle a);
    friend EG1DAngle operator+(EG1DAngle a, EG1DAngle b);
    friend EG1DAngle operator-(EG1DAngle a, EG1DAngle b);
    friend EG1DAngle operator*(double m, EG1DAngle a);
    friend EG1DAngle operator*(EG1DAngle a, double m);
    friend EG1DAngle operator/(EG1DAngle a, double m);
    friend double operator/(EG1DAngle a, EG1DAngle b);
    EG1DAngle& operator+=(EG1DAngle a);
    EG1DAngle& operator-=(EG1DAngle a);
    EG1DAngle& operator*=(double m);
    EG1DAngle& operator/=(double m);

    // Trigonmetric functions (not necessary but slightly more convenient).
    friend double sin(EG1DAngle a);
    friend double cos(EG1DAngle a);
    friend double tan(EG1DAngle a);

    // Return the angle normalized to the range (-180, 180] degrees.
    EG1DAngle Normalized() const;

    // Normalize this angle to the range (-180, 180] degrees.
    void Normalize();

    // When EG1DAngle is used as a key in one of the btree container types
    // (util/btree), indicate that linear rather than binary search should be
    // used.  This is much faster when the comparison function is cheap.
    typedef std::true_type goog_btree_prefer_linear_node_search;

private:
    explicit IFNDEF_SWIG(constexpr) EG1DAngle(double radians) : radians_(radians) {}
    double radians_;
};


//////////////////   Implementation details follow   ////////////////////


inline constexpr EG1DAngle EG1DAngle::Infinity() {
  return EG1DAngle(std::numeric_limits<double>::infinity());
}

inline constexpr EG1DAngle EG1DAngle::Zero() {
  return EG1DAngle(0);
}

inline constexpr double EG1DAngle::radians() const {
  return radians_;
}

inline constexpr double EG1DAngle::degrees() const {
  return (180 / M_PI) * radians_;
}

// Note that the E5, E6, and E7 conversion involve two multiplications rather
// than one.  This is mainly for backwards compatibility (changing this would
// break many tests), but it does have the nice side effect that conversions
// between Degrees, E6, and E7 are exact when the arguments are integers.

inline int32 EG1DAngle::e5() const {
  return MathUtil::FastIntRound(1e5 * degrees());
}

inline int32 EG1DAngle::e6() const {
  return MathUtil::FastIntRound(1e6 * degrees());
}

inline int32 EG1DAngle::e7() const {
  return MathUtil::FastIntRound(1e7 * degrees());
}

inline EG1DAngle EG1DAngle::abs() const {
  return EG1DAngle(std::fabs(radians_));
}

inline bool operator==(EG1DAngle x, EG1DAngle y) {
  return x.radians() == y.radians();
}

inline bool operator!=(EG1DAngle x, EG1DAngle y) {
  return x.radians() != y.radians();
}

inline bool operator<(EG1DAngle x, EG1DAngle y) {
  return x.radians() < y.radians();
}

inline bool operator>(EG1DAngle x, EG1DAngle y) {
  return x.radians() > y.radians();
}

inline bool operator<=(EG1DAngle x, EG1DAngle y) {
  return x.radians() <= y.radians();
}

inline bool operator>=(EG1DAngle x, EG1DAngle y) {
  return x.radians() >= y.radians();
}

inline EG1DAngle operator-(EG1DAngle a) {
  return EG1DAngle::Radians(-a.radians());
}

inline EG1DAngle operator+(EG1DAngle a, EG1DAngle b) {
  return EG1DAngle::Radians(a.radians() + b.radians());
}

inline EG1DAngle operator-(EG1DAngle a, EG1DAngle b) {
  return EG1DAngle::Radians(a.radians() - b.radians());
}

inline EG1DAngle operator*(double m, EG1DAngle a) {
  return EG1DAngle::Radians(m * a.radians());
}

inline EG1DAngle operator*(EG1DAngle a, double m) {
  return EG1DAngle::Radians(m * a.radians());
}

inline EG1DAngle operator/(EG1DAngle a, double m) {
  return EG1DAngle::Radians(a.radians() / m);
}

inline double operator/(EG1DAngle a, EG1DAngle b) {
  return a.radians() / b.radians();
}

inline EG1DAngle& EG1DAngle::operator+=(EG1DAngle a) {
  radians_ += a.radians();
  return *this;
}

inline EG1DAngle& EG1DAngle::operator-=(EG1DAngle a) {
  radians_ -= a.radians();
  return *this;
}

inline EG1DAngle& EG1DAngle::operator*=(double m) {
  radians_ *= m;
  return *this;
}

inline EG1DAngle& EG1DAngle::operator/=(double m) {
  radians_ /= m;
  return *this;
}

inline double sin(EG1DAngle a) {
  return sin(a.radians());
}

inline double cos(EG1DAngle a) {
  return cos(a.radians());
}

inline double tan(EG1DAngle a) {
  return tan(a.radians());
}

inline constexpr EG1DAngle EG1DAngle::Radians(double radians) {
  return EG1DAngle(radians);
}

inline constexpr EG1DAngle EG1DAngle::Degrees(double degrees) {
  return EG1DAngle((M_PI / 180) * degrees);
}

inline constexpr EG1DAngle EG1DAngle::E5(int32 e5) {
  return Degrees(1e-5 * e5);
}

inline constexpr EG1DAngle EG1DAngle::E6(int32 e6) {
  return Degrees(1e-6 * e6);
}

inline constexpr EG1DAngle EG1DAngle::E7(int32 e7) {
  return Degrees(1e-7 * e7);
}

inline constexpr EG1DAngle EG1DAngle::UnsignedE6(uint32 e6) {
  return E6(static_cast<int32>(e6));
}

inline constexpr EG1DAngle EG1DAngle::UnsignedE7(uint32 e7) {
  return E7(static_cast<int32>(e7));
}

// Writes the angle in degrees with 7 digits of precision after the
// decimal point, e.g. "17.3745904".
std::ostream& operator<<(std::ostream& os, EG1DAngle a);

#undef IFNDEF_SWIG

#endif  // EG_EG1DAngle_H_
