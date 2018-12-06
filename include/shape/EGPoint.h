//
// Created by pgl on 2018/12/6.
//

#ifndef EARTHGRID_EGPOINT_H
#define EARTHGRID_EGPOINT_H


#include "../util/math/vector.h"  // IWYU pragma: export
#include "../util/math/vector3_hash.h"

// An EGPoint represents a point on the unit sphere as a 3D vector.  Usually
// points are normalized to be unit length, but some methods do not require
// this.  See util/math/vector.h for the methods available.  Among other
// things, there are overloaded operators that make it convenient to write
// arithmetic expressions (e.g. (1-x)*p1 + x*p2).
using EGPoint = Vector3_d;

// EGPointHash can be used with standard containers (e.g., unordered_set) or
// nonstandard extensions (e.g., hash_map).  It is defined such that if two
// EGPoints compare equal to each other, they have the same hash.  (This
// requires that positive and negative zero hash to the same value.)
using EGPointHash = GoodFastHash<EGPoint>;

#endif //EARTHGRID_EGPOINT_H
