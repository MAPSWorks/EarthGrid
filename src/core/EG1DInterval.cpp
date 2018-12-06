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

#include "core/EG1DInterval.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

//#include "s2/base/logging.h"

using std::fabs;
using std::max;

EG1DInterval EG1DInterval::FromPoint(double p) {
  if (p == -M_PI) p = M_PI;
  return EG1DInterval(p, p, ARGS_CHECKED);
}

double EG1DInterval::GetCenter() const {
  double center = 0.5 * (lo() + hi());
  if (!is_inverted()) return center;
  // Return the center in the range (-Pi, Pi].
  return (center <= 0) ? (center + M_PI) : (center - M_PI);
}

double EG1DInterval::GetLength() const {
  double length = hi() - lo();
  if (length >= 0) return length;
  length += 2 * M_PI;
  // Empty intervals have a negative length.
  return (length > 0) ? length : -1;
}

EG1DInterval EG1DInterval::Complement() const {
  if (lo() == hi()) return Full();   // Singleton.
  return EG1DInterval(hi(), lo(), ARGS_CHECKED);  // Handles empty and full.
}

double EG1DInterval::GetComplementCenter() const {
  if (lo() != hi()) {
    return Complement().GetCenter();
  } else {  // Singleton.
    return (hi() <= 0) ? (hi() + M_PI) : (hi() - M_PI);
  }
}

bool EG1DInterval::FastContains(double p) const {
  if (is_inverted()) {
    return (p >= lo() || p <= hi()) && !is_empty();
  } else {
    return p >= lo() && p <= hi();
  }
}

bool EG1DInterval::Contains(double p) const {
  // Works for empty, full, and singleton intervals.
  EG_DCHECK_LE(fabs(p), M_PI);
  if (p == -M_PI) p = M_PI;
  return FastContains(p);
}

bool EG1DInterval::InteriorContains(double p) const {
  // Works for empty, full, and singleton intervals.
  EG_DCHECK_LE(fabs(p), M_PI);
  if (p == -M_PI) p = M_PI;

  if (is_inverted()) {
    return p > lo() || p < hi();
  } else {
    return (p > lo() && p < hi()) || is_full();
  }
}

bool EG1DInterval::Contains(const EG1DInterval& y) const {
  // It might be helpful to compare the structure of these tests to
  // the simpler Contains(double) method above.

  if (is_inverted()) {
    if (y.is_inverted()) return y.lo() >= lo() && y.hi() <= hi();
    return (y.lo() >= lo() || y.hi() <= hi()) && !is_empty();
  } else {
    if (y.is_inverted()) return is_full() || y.is_empty();
    return y.lo() >= lo() && y.hi() <= hi();
  }
}

bool EG1DInterval::InteriorContains(const EG1DInterval& y) const {
  if (is_inverted()) {
    if (!y.is_inverted()) return y.lo() > lo() || y.hi() < hi();
    return (y.lo() > lo() && y.hi() < hi()) || y.is_empty();
  } else {
    if (y.is_inverted()) return is_full() || y.is_empty();
    return (y.lo() > lo() && y.hi() < hi()) || is_full();
  }
}

bool EG1DInterval::Intersects(const EG1DInterval& y) const {
  if (is_empty() || y.is_empty()) return false;
  if (is_inverted()) {
    // Every non-empty inverted interval contains Pi.
    return y.is_inverted() || y.lo() <= hi() || y.hi() >= lo();
  } else {
    if (y.is_inverted()) return y.lo() <= hi() || y.hi() >= lo();
    return y.lo() <= hi() && y.hi() >= lo();
  }
}

bool EG1DInterval::InteriorIntersects(const EG1DInterval& y) const {
  if (is_empty() || y.is_empty() || lo() == hi()) return false;
  if (is_inverted()) {
    return y.is_inverted() || y.lo() < hi() || y.hi() > lo();
  } else {
    if (y.is_inverted()) return y.lo() < hi() || y.hi() > lo();
    return (y.lo() < hi() && y.hi() > lo()) || is_full();
  }
}

inline static double PositiveDistance(double a, double b) {
  // Compute the distance from "a" to "b" in the range [0, 2*Pi).
  // This is equivalent to (remainder(b - a - M_PI, 2 * M_PI) + M_PI),
  // except that it is more numerically stable (it does not lose
  // precision for very small positive distances).
  double d = b - a;
  if (d >= 0) return d;
  // We want to ensure that if b == Pi and a == (-Pi + eps),
  // the return result is approximately 2*Pi and not zero.
  return (b + M_PI) - (a - M_PI);
}

double EG1DInterval::GetDirectedHausdorffDistance(const EG1DInterval& y) const {
  if (y.Contains(*this)) return 0.0;  // this includes the case *this is empty
  if (y.is_empty()) return M_PI;  // maximum possible distance on S1

  double y_complement_center = y.GetComplementCenter();
  if (Contains(y_complement_center)) {
    return PositiveDistance(y.hi(), y_complement_center);
  } else {
    // The Hausdorff distance is realized by either two hi() endpoints or two
    // lo() endpoints, whichever is farther apart.
    double hi_hi = EG1DInterval(y.hi(), y_complement_center).Contains(hi()) ?
        PositiveDistance(y.hi(), hi()) : 0;
    double lo_lo = EG1DInterval(y_complement_center, y.lo()).Contains(lo()) ?
        PositiveDistance(lo(), y.lo()) : 0;
    EG_DCHECK(hi_hi > 0 || lo_lo > 0);
    return max(hi_hi, lo_lo);
  }
}

void EG1DInterval::AddPoint(double p) {
  EG_DCHECK_LE(fabs(p), M_PI);
  if (p == -M_PI) p = M_PI;

  if (FastContains(p)) return;
  if (is_empty()) {
    set_hi(p);
    set_lo(p);
  } else {
    // Compute distance from p to each endpoint.
    double dlo = PositiveDistance(p, lo());
    double dhi = PositiveDistance(hi(), p);
    if (dlo < dhi) {
      set_lo(p);
    } else {
      set_hi(p);
    }
    // Adding a point can never turn a non-full interval into a full one.
  }
}

double EG1DInterval::Project(double p) const {
  EG_DCHECK(!is_empty());
  EG_DCHECK_LE(fabs(p), M_PI);
  if (p == -M_PI) p = M_PI;
  if (FastContains(p)) return p;
  // Compute distance from p to each endpoint.
  double dlo = PositiveDistance(p, lo());
  double dhi = PositiveDistance(hi(), p);
  return (dlo < dhi) ? lo() : hi();
}

EG1DInterval EG1DInterval::FromPointPair(double p1, double p2) {
  EG_DCHECK_LE(fabs(p1), M_PI);
  EG_DCHECK_LE(fabs(p2), M_PI);
  if (p1 == -M_PI) p1 = M_PI;
  if (p2 == -M_PI) p2 = M_PI;
  if (PositiveDistance(p1, p2) <= M_PI) {
    return EG1DInterval(p1, p2, ARGS_CHECKED);
  } else {
    return EG1DInterval(p2, p1, ARGS_CHECKED);
  }
}

EG1DInterval EG1DInterval::Expanded(double margin) const {
  if (margin >= 0) {
    if (is_empty()) return *this;
    // Check whether this interval will be full after expansion, allowing
    // for a 1-bit rounding error when computing each endpoint.
    if (GetLength() + 2 * margin + 2 * DBL_EPSILON >= 2 * M_PI) return Full();
  } else {
    if (is_full()) return *this;
    // Check whether this interval will be empty after expansion, allowing
    // for a 1-bit rounding error when computing each endpoint.
    if (GetLength() + 2 * margin - 2 * DBL_EPSILON <= 0) return Empty();
  }
  EG1DInterval result(remainder(lo() - margin, 2*M_PI),
                    remainder(hi() + margin, 2*M_PI));
  if (result.lo() <= -M_PI) result.set_lo(M_PI);
  return result;
}

EG1DInterval EG1DInterval::Union(const EG1DInterval& y) const {
  // The y.is_full() case is handled correctly in all cases by the code
  // below, but can follow three separate code paths depending on whether
  // this interval is inverted, is non-inverted but contains Pi, or neither.

  if (y.is_empty()) return *this;
  if (FastContains(y.lo())) {
    if (FastContains(y.hi())) {
      // Either this interval contains y, or the union of the two
      // intervals is the Full() interval.
      if (Contains(y)) return *this;  // is_full() code path
      return Full();
    }
    return EG1DInterval(lo(), y.hi(), ARGS_CHECKED);
  }
  if (FastContains(y.hi())) return EG1DInterval(y.lo(), hi(), ARGS_CHECKED);

  // This interval contains neither endpoint of y.  This means that either y
  // contains all of this interval, or the two intervals are disjoint.
  if (is_empty() || y.FastContains(lo())) return y;

  // Check which pair of endpoints are closer together.
  double dlo = PositiveDistance(y.hi(), lo());
  double dhi = PositiveDistance(hi(), y.lo());
  if (dlo < dhi) {
    return EG1DInterval(y.lo(), hi(), ARGS_CHECKED);
  } else {
    return EG1DInterval(lo(), y.hi(), ARGS_CHECKED);
  }
}

EG1DInterval EG1DInterval::Intersection(const EG1DInterval& y) const {
  // The y.is_full() case is handled correctly in all cases by the code
  // below, but can follow three separate code paths depending on whether
  // this interval is inverted, is non-inverted but contains Pi, or neither.

  if (y.is_empty()) return Empty();
  if (FastContains(y.lo())) {
    if (FastContains(y.hi())) {
      // Either this interval contains y, or the region of intersection
      // consists of two disjoint subintervals.  In either case, we want
      // to return the shorter of the two original intervals.
      if (y.GetLength() < GetLength()) return y;  // is_full() code path
      return *this;
    }
    return EG1DInterval(y.lo(), hi(), ARGS_CHECKED);
  }
  if (FastContains(y.hi())) return EG1DInterval(lo(), y.hi(), ARGS_CHECKED);

  // This interval contains neither endpoint of y.  This means that either y
  // contains all of this interval, or the two intervals are disjoint.

  if (y.FastContains(lo())) return *this;  // is_empty() okay here
  EG_DCHECK(!Intersects(y));
  return Empty();
}

bool EG1DInterval::ApproxEquals(const EG1DInterval& y, double max_error) const {
  // Full and empty intervals require special cases because the "endpoints"
  // are considered to be positioned arbitrarily.
  if (is_empty()) return y.GetLength() <= 2 * max_error;
  if (y.is_empty()) return GetLength() <= 2 * max_error;
  if (is_full()) return y.GetLength() >= 2 * (M_PI - max_error);
  if (y.is_full()) return GetLength() >= 2 * (M_PI - max_error);

  // The purpose of the last test below is to verify that moving the endpoints
  // does not invert the interval, e.g. [-1e20, 1e20] vs. [1e20, -1e20].
  return (fabs(remainder(y.lo() - lo(), 2 * M_PI)) <= max_error &&
          fabs(remainder(y.hi() - hi(), 2 * M_PI)) <= max_error &&
          fabs(GetLength() - y.GetLength()) <= 2 * max_error);
}
