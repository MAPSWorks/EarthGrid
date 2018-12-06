// Copyright 2013 Google Inc. All Rights Reserved.
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

#ifndef EG_S1CHORD_ANGLE_H_
#define EG_S1CHORD_ANGLE_H_

#include <cmath>
#include <limits>
#include <ostream>
#include <type_traits>

#include "_fp_contract_off.h"
#include "core/EG1DAngle.h"
#include "core/EGPointUtil.h"

// EG1DChordAngle represents the angle subtended by a chord (i.e., the straight
// line segment connecting two points on the sphere).  Its representation
// makes it very efficient for computing and comparing distances, but unlike
// EG1DAngle it is only capable of representing angles between 0 and Pi radians.
// Generally, EG1DChordAngle should only be used in loops where many angles need
// to be calculated and compared.  Otherwise it is simpler to use EG1DAngle.
//
// EG1DChordAngle also loses some accuracy as the angle approaches Pi radians.
// Specifically, the representation of (Pi - x) radians has an error of about
// (1e-15 / x), with a maximum error of about 2e-8 radians (about 13cm on the
// Earth's surface).  For comparison, for angles up to 90 degrees (10000km)
// the worst-case representation error is about 2e-16 radians (1 nanometer),
// which is about the same as EG1DAngle.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class EG1DChordAngle {
public:
    // The default constructor yields a zero angle.  This is useful for STL
    // containers and class methods with output arguments.
    EG1DChordAngle() : length2_(0) {}

    // Construct the EG1DChordAngle corresponding to the distance between the two
    // given points.  The points must be unit length.
    EG1DChordAngle(const EGPoint& x, const EGPoint& y);

    // Return the zero chord angle.
    static EG1DChordAngle Zero();

    // Return a chord angle of 90 degrees (a "right angle").
    static EG1DChordAngle Right();

    // Return a chord angle of 180 degrees (a "straight angle").  This is the
    // maximum finite chord angle.
    static EG1DChordAngle Straight();

    // Return a chord angle larger than any finite chord angle.  The only valid
    // operations on Infinity() are comparisons, EG1DAngle conversions, and
    // Successor() / Predecessor().
    static EG1DChordAngle Infinity();

    // Return a chord angle smaller than Zero().  The only valid operations on
    // Negative() are comparisons, EG1DAngle conversions, and Successor() /
    // Predecessor().
    static EG1DChordAngle Negative();

    // Conversion from an EG1DAngle.  Angles outside the range [0, Pi] are handled
    // as follows: Infinity() is mapped to Infinity(), negative angles are
    // mapped to Negative(), and finite angles larger than Pi are mapped to
    // Straight().
    //
    // Note that this operation is relatively expensive and should be avoided.
    // To use EG1DChordAngle effectively, you should structure your code so that
    // input arguments are converted to EG1DChordAngles at the beginning of your
    // algorithm, and results are converted back to EG1DAngles only at the end.
    explicit EG1DChordAngle(EG1DAngle angle);

    // Convenience methods implemented by converting from an EG1DAngle.
    static EG1DChordAngle Radians(double radians);
    static EG1DChordAngle Degrees(double degrees);
    static EG1DChordAngle E5(int32 e5);
    static EG1DChordAngle E6(int32 e6);
    static EG1DChordAngle E7(int32 e7);

    // Construct an EG1DChordAngle that is an upper bound on the given EG1DAngle,
    // i.e. such that FastUpperBoundFrom(x).ToAngle() >= x.  Unlike the EG1DAngle
    // constructor above, this method is very fast, and the bound is accurate to
    // within 1% for distances up to about 3100km on the Earth's surface.
    static EG1DChordAngle FastUpperBoundFrom(EG1DAngle angle);

    // Construct an EG1DChordAngle from the squared chord length.  Note that the
    // argument is automatically clamped to a maximum of 4.0 to handle possible
    // roundoff errors.  The argument must be non-negative.
    static EG1DChordAngle FromLength2(double length2);

    // Converts to an EG1DAngle.  Can be used just like an EG1DAngle constructor:
    //
    //   EG1DChordAngle x = ...;
    //   return EG1DAngle(x);
    //
    // Infinity() is converted to EG1DAngle::Infinity(), and Negative() is
    // converted to an unspecified negative EG1DAngle.
    //
    // Note that the conversion uses trigonometric functions and therefore
    // should be avoided in inner loops.
    explicit operator EG1DAngle() const;

    // Converts to an EG1DAngle (equivalent to the operator above).
    EG1DAngle ToAngle() const;

    // Convenience methods implemented by calling ToAngle() first.  Note that
    // because of the EG1DAngle conversion these methods are relatively expensive
    // (despite their lowercase names), so the results should be cached if they
    // are needed inside loops.
    double radians() const;
    double degrees() const;
    int32 e5() const;
    int32 e6() const;
    int32 e7() const;

    // All operators and functions are declared here so that we can put them all
    // in one place.  (The compound assignment operators must be put here.)

    // Comparison operators.
    friend bool operator==(EG1DChordAngle x, EG1DChordAngle y);
    friend bool operator!=(EG1DChordAngle x, EG1DChordAngle y);
    friend bool operator<(EG1DChordAngle x, EG1DChordAngle y);
    friend bool operator>(EG1DChordAngle x, EG1DChordAngle y);
    friend bool operator<=(EG1DChordAngle x, EG1DChordAngle y);
    friend bool operator>=(EG1DChordAngle x, EG1DChordAngle y);

    // Comparison predicates.
    bool is_zero() const;
    bool is_negative() const;
    bool is_infinity() const;
    bool is_special() const;  // Negative or infinity.

    // Only addition and subtraction of EG1DChordAngles is supported.  These
    // methods add or subtract the corresponding EG1DAngles, and clamp the result
    // to the range [0, Pi].  Both arguments must be non-negative and
    // non-infinite.
    //
    // REQUIRES: !a.is_special() && !b.is_special()
    friend EG1DChordAngle operator+(EG1DChordAngle a, EG1DChordAngle b);
    friend EG1DChordAngle operator-(EG1DChordAngle a, EG1DChordAngle b);
    EG1DChordAngle& operator+=(EG1DChordAngle a);
    EG1DChordAngle& operator-=(EG1DChordAngle a);

    // Trigonmetric functions.  It is more accurate and efficient to call these
    // rather than first converting to an EG1DAngle.
    friend double sin(EG1DChordAngle a);
    friend double cos(EG1DChordAngle a);
    friend double tan(EG1DChordAngle a);

    // Returns sin(a)**2, but computed more efficiently.
    friend double sin2(EG1DChordAngle a);

    // The squared length of the chord.  (Most clients will not need this.)
    double length2() const { return length2_; }

    // Returns the smallest representable EG1DChordAngle larger than this object.
    // This can be used to convert a "<" comparison to a "<=" comparison.  For
    // example:
    //
    //   EGClosestEdgeQuery query(...);
    //   EG1DChordAngle limit = ...;
    //   if (query.IsDistanceLess(target, limit.Successor())) {
    //     // Distance to "target" is less than or equal to "limit".
    //   }
    //
    // Note the following special cases:
    //   Negative().Successor() == Zero()
    //   Straight().Successor() == Infinity()
    //   Infinity().Successor() == Infinity()
    EG1DChordAngle Successor() const;

    // Like Successor(), but returns the largest representable EG1DChordAngle less
    // than this object.
    //
    // Note the following special cases:
    //   Infinity().Predecessor() == Straight()
    //   Zero().Predecessor() == Negative()
    //   Negative().Predecessor() == Negative()
    EG1DChordAngle Predecessor() const;

    // Returns a new EG1DChordAngle that has been adjusted by the given error
    // bound (which can be positive or negative).  "error" should be the value
    // returned by one of the error bound methods below.  For example:
    //    EG1DChordAngle a(x, y);
    //    EG1DChordAngle a1 = a.PlusError(a.GetEGPointConstructorMaxError());
    EG1DChordAngle PlusError(double error) const;

    // Return the maximum error in length2() for the EG1DChordAngle(x, y)
    // constructor, assuming that "x" and "y" are normalized to within the
    // bounds guaranteed by EGPoint::Normalize().  (The error is defined with
    // respect to the true distance after the points are projected to lie
    // exactly on the sphere.)
    double GetEGPointConstructorMaxError() const;

    // Return the maximum error in length2() for the EG1DAngle constructor.
    double GetEG1DAngleConstructorMaxError() const;

    // Return true if the internal representation is valid.  Negative() and
    // Infinity() are both considered valid.
    bool is_valid() const;

    // When EG1DChordAngle is used as a key in one of the btree container types
    // (util/btree), indicate that linear rather than binary search should be
    // used.  This is much faster when the comparison function is cheap.
    typedef std::true_type goog_btree_prefer_linear_node_search;

private:
    // EG1DChordAngles are represented by the squared chord length, which can
    // range from 0 to 4.  Infinity() uses an infinite squared length.
    explicit EG1DChordAngle(double length2) : length2_(length2) {
      EG_DCHECK(is_valid());
    }
    double length2_;
};


//////////////////   Implementation details follow   ////////////////////


inline EG1DChordAngle::EG1DChordAngle(const EGPoint& x, const EGPoint& y) {
  EG_DCHECK(EG::IsUnitLength(x));
  EG_DCHECK(EG::IsUnitLength(y));
  // The squared distance may slightly exceed 4.0 due to roundoff errors.
  // The maximum error in the result is 2 * DBL_EPSILON * length2_.
  length2_ = std::min(4.0, (x - y).Norm2());
  EG_DCHECK(is_valid());
}

inline EG1DChordAngle EG1DChordAngle::FromLength2(double length2) {
  return EG1DChordAngle(std::min(4.0, length2));
}

inline EG1DChordAngle EG1DChordAngle::Zero() {
  return EG1DChordAngle(0);
}

inline EG1DChordAngle EG1DChordAngle::Right() {
  return EG1DChordAngle(2);
}

inline EG1DChordAngle EG1DChordAngle::Straight() {
  return EG1DChordAngle(4);
}

inline EG1DChordAngle EG1DChordAngle::Infinity() {
  return EG1DChordAngle(std::numeric_limits<double>::infinity());
}

inline EG1DChordAngle EG1DChordAngle::Negative() {
  return EG1DChordAngle(-1);
}

inline EG1DChordAngle EG1DChordAngle::Radians(double radians) {
  return EG1DChordAngle(EG1DAngle::Radians(radians));
}

inline EG1DChordAngle EG1DChordAngle::Degrees(double degrees) {
  return EG1DChordAngle(EG1DAngle::Degrees(degrees));
}

inline EG1DChordAngle EG1DChordAngle::E5(int32 e5) {
  return EG1DChordAngle(EG1DAngle::E5(e5));
}

inline EG1DChordAngle EG1DChordAngle::E6(int32 e6) {
  return EG1DChordAngle(EG1DAngle::E6(e6));
}

inline EG1DChordAngle EG1DChordAngle::E7(int32 e7) {
  return EG1DChordAngle(EG1DAngle::E7(e7));
}

inline EG1DChordAngle EG1DChordAngle::FastUpperBoundFrom(EG1DAngle angle) {
  // This method uses the distance along the surface of the sphere as an upper
  // bound on the distance through the sphere's interior.
  return EG1DChordAngle::FromLength2(angle.radians() * angle.radians());
}

inline EG1DChordAngle::operator EG1DAngle() const {
  return ToAngle();
}

inline double EG1DChordAngle::radians() const {
  return ToAngle().radians();
}

inline double EG1DChordAngle::degrees() const {
  return ToAngle().degrees();
}

inline int32 EG1DChordAngle::e5() const {
  return ToAngle().e5();
}

inline int32 EG1DChordAngle::e6() const {
  return ToAngle().e6();
}

inline int32 EG1DChordAngle::e7() const {
  return ToAngle().e7();
}

inline bool EG1DChordAngle::is_zero() const {
  return length2_ == 0;
}

inline bool EG1DChordAngle::is_negative() const {
  // TODO(ericv): Consider stricter check here -- only allow Negative().
  return length2_ < 0;
}

inline bool EG1DChordAngle::is_infinity() const {
  return length2_ == std::numeric_limits<double>::infinity();
}

inline bool EG1DChordAngle::is_special() const {
  return is_negative() || is_infinity();
}

inline bool operator==(EG1DChordAngle x, EG1DChordAngle y) {
  return x.length2() == y.length2();
}

inline bool operator!=(EG1DChordAngle x, EG1DChordAngle y) {
  return x.length2() != y.length2();
}

inline bool operator<(EG1DChordAngle x, EG1DChordAngle y) {
  return x.length2() < y.length2();
}

inline bool operator>(EG1DChordAngle x, EG1DChordAngle y) {
  return x.length2() > y.length2();
}

inline bool operator<=(EG1DChordAngle x, EG1DChordAngle y) {
  return x.length2() <= y.length2();
}

inline bool operator>=(EG1DChordAngle x, EG1DChordAngle y) {
  return x.length2() >= y.length2();
}

inline EG1DChordAngle& EG1DChordAngle::operator+=(EG1DChordAngle a) {
  return (*this = *this + a);
}

inline EG1DChordAngle& EG1DChordAngle::operator-=(EG1DChordAngle a) {
  return (*this = *this - a);
}

// Outputs the chord angle as the equivalent EG1DAngle.
std::ostream& operator<<(std::ostream& os, EG1DChordAngle a);

#endif  // EG_S1CHORD_ANGLE_H_
