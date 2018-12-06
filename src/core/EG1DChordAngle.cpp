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

#include "core/EG1DChordAngle.h"

#include <cfloat>
#include <cmath>

#include "core/EG1DAngle.h"
#include "core/EGPointUtil.h"

using std::max;
using std::min;
// Android with gnustl has ::nextafter but not std::nextafter.
// https://github.com/android-ndk/ndk/issues/82
// Check for gnustl with _GLIBCXX_CMATH, which is its cmath include
// guard.
#if !defined(__ANDROID__) || !defined(_GLIBCXX_CMATH)
using std::nextafter;
#endif

static constexpr double kMaxLength2 = 4.0;

EG1DChordAngle::EG1DChordAngle(EG1DAngle angle) {
    if (angle.radians() < 0) {
        *this = Negative();
    } else if (angle == EG1DAngle::Infinity()) {
        *this = Infinity();
    } else {
        // The chord length is 2 * sin(angle / 2).
        double length = 2 * sin(0.5 * min(M_PI, angle.radians()));
        length2_ = length * length;
    }
    EG1DAngle(is_valid());
}

EG1DAngle EG1DChordAngle::ToAngle() const {
    if (is_negative()) return EG1DAngle::Radians(-1);
    if (is_infinity()) return EG1DAngle::Infinity();
    return EG1DAngle::Radians(2 * asin(0.5 * sqrt(length2_)));
}

bool EG1DChordAngle::is_valid() const {
    return (length2_ >= 0 && length2_ <= kMaxLength2) || is_special();
}

EG1DChordAngle EG1DChordAngle::Successor() const {
    if (length2_ >= kMaxLength2) return Infinity();
    if (length2_ < 0.0) return Zero();
    return EG1DChordAngle(nextafter(length2_, 10.0));
}

EG1DChordAngle EG1DChordAngle::Predecessor() const {
    if (length2_ <= 0.0) return Negative();
    if (length2_ > kMaxLength2) return Straight();
    return EG1DChordAngle(nextafter(length2_, -10.0));
}

EG1DChordAngle EG1DChordAngle::PlusError(double error) const {
    // If angle is Negative() or Infinity(), don't change it.
    // Otherwise clamp it to the valid range.
    if (is_special()) return *this;
    return EG1DChordAngle(max(0.0, min(kMaxLength2, length2_ + error)));
}

double EG1DChordAngle::GetEGPointConstructorMaxError() const {
    // There is a relative error of 2.5 * DBL_EPSILON when computing the squared
    // distance, plus a relative error of 2 * DBL_EPSILON and an absolute error
    // of (16 * DBL_EPSILON**2) because the lengths of the input points may
    // differ from 1 by up to (2 * DBL_EPSILON) each.  (This is the maximum
    // length error in S2Point::Normalize.)
    return 4.5 * DBL_EPSILON * length2_ + 16 * DBL_EPSILON * DBL_EPSILON;
}

double EG1DChordAngle::GetEG1DAngleConstructorMaxError() const {
    // Assuming that an accurate math library is being used, the sin() call and
    // the multiply each have a relative error of 0.5 * DBL_EPSILON.
    return DBL_EPSILON * length2_;
}

EG1DChordAngle operator+(EG1DChordAngle a, EG1DChordAngle b) {
    // Note that this method is much more efficient than converting the chord
    // angles to EG1DAngles and adding those.  It requires only one square root
    // plus a few additions and multiplications.
    EG_DCHECK(!a.is_special());
    EG_DCHECK(!b.is_special());

    // Optimization for the common case where "b" is an error tolerance
    // parameter that happens to be set to zero.
    double a2 = a.length2(), b2 = b.length2();
    if (b2 == 0) return a;

    // Clamp the angle sum to at most 180 degrees.
    if (a2 + b2 >= kMaxLength2) return EG1DChordAngle::Straight();

    // Let "a" and "b" be the (non-squared) chord lengths, and let c = a+b.
    // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
    // Then the formula below can be derived from c = 2 * sin(A+B) and the
    // relationships   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
    //                 cos(X) = sqrt(1 - sin^2(X)) .

    double x = a2 * (1 - 0.25 * b2);  // is_valid() => non-negative
    double y = b2 * (1 - 0.25 * a2);  // is_valid() => non-negative
    return EG1DChordAngle(min(kMaxLength2, x + y + 2 * sqrt(x * y)));
}

EG1DChordAngle operator-(EG1DChordAngle a, EG1DChordAngle b) {
    // See comments in operator+().
    EG_DCHECK(!a.is_special());
    EG_DCHECK(!b.is_special());
    double a2 = a.length2(), b2 = b.length2();
    if (b2 == 0) return a;
    if (a2 <= b2) return EG1DChordAngle::Zero();
    double x = a2 * (1 - 0.25 * b2);
    double y = b2 * (1 - 0.25 * a2);
    return EG1DChordAngle(max(0.0, x + y - 2 * sqrt(x * y)));
}

double sin2(EG1DChordAngle a) {
    EG_DCHECK(!a.is_special());
    // Let "a" be the (non-squared) chord length, and let A be the corresponding
    // half-angle (a = 2*sin(A)).  The formula below can be derived from:
    //   sin(2*A) = 2 * sin(A) * cos(A)
    //   cos^2(A) = 1 - sin^2(A)
    // This is much faster than converting to an angle and computing its sine.
    return a.length2() * (1 - 0.25 * a.length2());
}

double sin(EG1DChordAngle a) {
    return sqrt(sin2(a));
}

double cos(EG1DChordAngle a) {
    // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
    EG_DCHECK(!a.is_special());
    return 1 - 0.5 * a.length2();
}

double tan(EG1DChordAngle a) {
    return sin(a) / cos(a);
}

std::ostream& operator<<(std::ostream& os, EG1DChordAngle a) {
    return os << a.ToAngle();
}
