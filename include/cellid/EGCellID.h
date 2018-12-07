//
// Created by pgl on 2018/12/6.
//

#ifndef EARTHGRID_EGCellID_H
#define EARTHGRID_EGCellID_H

#include <cstddef>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "base/integral_types.h"
#include "base/logging.h"
#include "base/port.h"
#include "_fp_contract_off.h"
#include "core/r2.h"
#include "core/r2rect.h"
#include "core/EG1DAngle.h"
#include "core/EGCoords.h"
#include "third_party/absl/strings/string_view.h"
#include "util/bits/bits.h"
#include "util/coding/coder.h"

class  EGLatLng;

#ifndef SWIG
#define IFNDEF_SWIG(x) x
#else
#define IFNDEF_SWIG(x)
#endif

// An EGCellID is a 64-bit unsigned integer that uniquely identifies a
// cell in the EG cell decomposition.  It has the following format:
//
//   id = [face][face_pos]
//
//   face:     a 3-bit number (range 0..5) encoding the cube face.
//
//   face_pos: a 61-bit number encoding the position of the center of this
//             cell along the Hilbert curve over this face (see the Wiki
//             pages for details).
//
// Sequentially increasing cell ids follow a continuous space-filling curve
// over the entire sphere.  They have the following properties:
//
//  - The id of a cell at level k consists of a 3-bit face number followed
//    by k bit pairs that recursively select one of the four children of
//    each cell.  The next bit is always 1, and all other bits are 0.
//    Therefore, the level of a cell is determined by the position of its
//    lowest-numbered bit that is turned on (for a cell at level k, this
//    position is 2 * (kMaxLevel - k).)
//
//  - The id of a parent cell is at the midpoint of the range of ids spanned
//    by its children (or by its descendants at any level).
//
// Leaf cells are often used to represent points on the unit sphere, and
// this class provides methods for converting directly between these two
// representations.  For cells that represent 2D regions rather than
// discrete point, it is better to use the EGCell class.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class EGCellID {
public:
    // The extra position bit (61 rather than 60) let us encode each cell as its
    // Hilbert curve position at the cell center (which is halfway along the
    // portion of the Hilbert curve that fills that cell).
    static const int kFaceBits = 3;
    static const int kNumFaces = 6;
    static const int kMaxLevel = EG::kMaxCellLevel;  // Valid levels: 0..kMaxLevel
    static const int kPosBits = 2 * kMaxLevel + 1;
    static const int kMaxSize = 1 << kMaxLevel;

    explicit IFNDEF_SWIG(constexpr) EGCellID(uint64 id) : id_(id) {}

    // Construct a leaf cell containing the given point "p".  Usually there is
    // is exactly one such cell, but for points along the edge of a cell, any
    // adjacent cell may be (deterministically) chosen.  This is because
    // EGCellIDs are considered to be closed sets.  The returned cell will
    // always contain the given point, i.e.
    //
    //   EGCell(EGCellID(p)).Contains(p)
    //
    // is always true.  The point "p" does not need to be normalized.
    //
    // If instead you want every point to be contained by exactly one EGCell,
    // you will need to convert the EGCellIDs to EGLoops (which implement point
    // containment this way).
    explicit EGCellID(const EGPoint& p);

    // Construct a leaf cell containing the given normalized EGLatLng.
    explicit EGCellID(const EGLatLng& ll);

    // The default constructor returns an invalid cell id.
    IFNDEF_SWIG(constexpr) EGCellID() : id_(0) {}
    static constexpr EGCellID None() { return EGCellID(); }

    // Returns an invalid cell id guaranteed to be larger than any
    // valid cell id.  Useful for creating indexes.
    static constexpr EGCellID Sentinel() { return EGCellID(~uint64{0}); }

    // Return the cell corresponding to a given EG cube face.
    static EGCellID FromFace(int face);

    // Return a cell given its face (range 0..5), Hilbert curve position within
    // that face (an unsigned integer with EGCellID::kPosBits bits), and level
    // (range 0..kMaxLevel).  The given position will be modified to correspond
    // to the Hilbert curve position at the center of the returned cell.  This
    // is a static function rather than a constructor in order to indicate what
    // the arguments represent.
    static EGCellID FromFacePosLevel(int face, uint64 pos, int level);

    // Return the direction vector corresponding to the center of the given
    // cell.  The vector returned by ToPointRaw is not necessarily unit length.
    // This method returns the same result as EGCell::GetCenter().
    //
    // The maximum directional error in ToPoint() (compared to the exact
    // mathematical result) is 1.5 * DBL_EPSILON radians, and the maximum length
    // error is 2 * DBL_EPSILON (the same as Normalize).
    EGPoint ToPoint() const { return ToPointRaw().Normalize(); }
    EGPoint ToPointRaw() const;

    // Return the center of the cell in (s,t) coordinates (see EGcoords.h).
    R2Point GetCenterST() const;

    // Return the edge length of this cell in (s,t)-space.
    double GetSizeST() const;

    // Return the edge length in (s,t)-space of cells at the given level.
    static double GetSizeST(int level);

    // Return the bounds of this cell in (s,t)-space.
    R2Rect GetBoundST() const;

    // Return the center of the cell in (u,v) coordinates (see EGcoords.h).
    // Note that the center of the cell is defined as the point at which it is
    // recursively subdivided into four children; in general, it is not at the
    // midpoint of the (u,v) rectangle covered by the cell.
    R2Point GetCenterUV() const;

    // Return the bounds of this cell in (u,v)-space.
    R2Rect GetBoundUV() const;

    // Expand a rectangle in (u,v)-space so that it contains all points within
    // the given distance of the boundary, and return the smallest such
    // rectangle.  If the distance is negative, then instead shrink this
    // rectangle so that it excludes all points within the given absolute
    // distance of the boundary.
    //
    // Distances are measured *on the sphere*, not in (u,v)-space.  For example,
    // you can use this method to expand the (u,v)-bound of an EGCellID so that
    // it contains all points within 5km of the original cell.  You can then
    // test whether a point lies within the expanded bounds like this:
    //
    //   R2Point uv;
    //   if (EG::FaceXYZtoUV(face, point, &uv) && bound.Contains(uv)) { ... }
    //
    // Limitations:
    //
    //  - Because the rectangle is drawn on one of the six cube-face planes
    //    (i.e., {x,y,z} = +/-1), it can cover at most one hemisphere.  This
    //    limits the maximum amount that a rectangle can be expanded.  For
    //    example, EGCellID bounds can be expanded safely by at most 45 degrees
    //    (about 5000 km on the Earth's surface).
    //
    //  - The implementation is not exact for negative distances.  The resulting
    //    rectangle will exclude all points within the given distance of the
    //    boundary but may be slightly smaller than necessary.
    static R2Rect ExpandedByDistanceUV(const R2Rect& uv, EG1DAngle distance);

    // Return the (face, si, ti) coordinates of the center of the cell.  Note
    // that although (si,ti) coordinates span the range [0,2**31] in general,
    // the cell center coordinates are always in the range [1,2**31-1] and
    // therefore can be represented using a signed 32-bit integer.
    int GetCenterSiTi(int* psi, int* pti) const;

    // Return the EGLatLng corresponding to the center of the given cell.
    EGLatLng ToLatLng() const;

    // The 64-bit unique identifier for this cell.
    uint64 id() const { return id_; }

    // Return true if id() represents a valid cell.
    //
    // All methods require is_valid() to be true unless otherwise specified
    // (although not all methods enforce this).
    bool is_valid() const;

    // Which cube face this cell belongs to, in the range 0..5.
    int face() const;

    // The position of the cell center along the Hilbert curve over this face,
    // in the range 0..(2**kPosBits-1).
    uint64 pos() const;

    // Return the subdivision level of the cell (range 0..kMaxLevel).
    int level() const;

    // Return the edge length of this cell in (i,j)-space.
    int GetSizeIJ() const;

    // Like the above, but return the size of cells at the given level.
    static int GetSizeIJ(int level);

    // Return true if this is a leaf cell (more efficient than checking
    // whether level() == kMaxLevel).
    bool is_leaf() const;

    // Return true if this is a top-level face cell (more efficient than
    // checking whether level() == 0).
    bool is_face() const;

    // Return the child position (0..3) of this cell within its parent.
    // REQUIRES: level() >= 1.
    int child_position() const;

    // Return the child position (0..3) of this cell's ancestor at the given
    // level within its parent.  For example, child_position(1) returns the
    // position of this cell's level-1 ancestor within its top-level face cell.
    // REQUIRES: 1 <= level <= this->level().
    int child_position(int level) const;

    // These methods return the range of cell ids that are contained within this
    // cell (including itself).  The range is *inclusive* (i.e. test using >=
    // and <=) and the return values of both methods are valid leaf cell ids.
    // In other words, a.contains(b) if and only if
    //
    //     (b >= a.range_min() && b <= a.range_max())
    //
    // If you want to iterate through all the descendants of this cell at a
    // particular level, use child_begin(level) and child_end(level) instead.
    // Also see maximum_tile(), which can be used to iterate through a range of
    // cells using EGCellIDs at different levels that are as large as possible.
    //
    // If you need to convert the range to a semi-open interval [min, limit)
    // (e.g., in order to use a key-value store that only supports semi-open
    // range queries), do not attempt to define "limit" as range_max.next().
    // The problem is that leaf EGCellIDs are 2 units apart, so the semi-open
    // interval [min, limit) includes an additional value (range_max.id() + 1)
    // which is happens to be a valid EGCellID about one-third of the time and
    // is *never* contained by this cell.  (It always correpsonds to a cell that
    // is larger than this one.)  You can define "limit" as (range_max.id() + 1)
    // if necessary (which is not always a valid EGCellID but can still be used
    // with FromToken/ToToken), or you can convert range_max() to the key space
    // of your key-value store and define "limit" as Successor(key).
    //
    // Note that Sentinel().range_min() == Sentinel.range_max() == Sentinel().
    EGCellID range_min() const;
    EGCellID range_max() const;

    // Return true if the given cell is contained within this one.
    bool contains(EGCellID other) const;

    // Return true if the given cell intersects this one.
    bool intersects(EGCellID other) const;

    // Return the cell at the previous level or at the given level (which must
    // be less than or equal to the current level).
    EGCellID parent() const;
    EGCellID parent(int level) const;

    // Return the immediate child of this cell at the given traversal order
    // position (in the range 0 to 3).  This cell must not be a leaf cell.
    EGCellID child(int position) const;

    // Iterator-style methods for traversing the immediate children of a cell or
    // all of the children at a given level (greater than or equal to the current
    // level).  Note that the end value is exclusive, just like standard STL
    // iterators, and may not even be a valid cell id.  You should iterate using
    // code like this:
    //
    //   for(EGCellID c = id.child_begin(); c != id.child_end(); c = c.next())
    //     ...
    //
    // The convention for advancing the iterator is "c = c.next()" rather
    // than "++c" to avoid possible confusion with incrementing the
    // underlying 64-bit cell id.
    EGCellID child_begin() const;
    EGCellID child_begin(int level) const;
    EGCellID child_end() const;
    EGCellID child_end(int level) const;

    // Return the next/previous cell at the same level along the Hilbert curve.
    // Works correctly when advancing from one face to the next, but
    // does *not* wrap around from the last face to the first or vice versa.
    EGCellID next() const;
    EGCellID prev() const;

    // This method advances or retreats the indicated number of steps along the
    // Hilbert curve at the current level, and returns the new position.  The
    // position is never advanced past End() or before Begin().
    EGCellID advance(int64 steps) const;

    // Returns the number of steps that this cell is from Begin(level()). The
    // return value is always non-negative.
    int64 distance_from_begin() const;

    // Like next() and prev(), but these methods wrap around from the last face
    // to the first and vice versa.  They should *not* be used for iteration in
    // conjunction with child_begin(), child_end(), Begin(), or End().  The
    // input must be a valid cell id.
    EGCellID next_wrap() const;
    EGCellID prev_wrap() const;

    // This method advances or retreats the indicated number of steps along the
    // Hilbert curve at the current level, and returns the new position.  The
    // position wraps between the first and last faces as necessary.  The input
    // must be a valid cell id.
    EGCellID advance_wrap(int64 steps) const;

    // Return the largest cell with the same range_min() and such that
    // range_max() < limit.range_min().  Returns "limit" if no such cell exists.
    // This method can be used to generate a small set of EGCellIDs that covers
    // a given range (a "tiling").  This example shows how to generate a tiling
    // for a semi-open range of leaf cells [start, limit):
    //
    //   for (EGCellID id = start.maximum_tile(limit);
    //        id != limit; id = id.next().maximum_tile(limit)) { ... }
    //
    // Note that in general the cells in the tiling will be of different sizes;
    // they gradually get larger (near the middle of the range) and then
    // gradually get smaller (as "limit" is approached).
    EGCellID maximum_tile(EGCellID limit) const;

    // Returns the level of the lowest common ancestor of this cell and "other",
    // that is, the maximum level such that parent(level) == other.parent(level).
    // Returns -1 if the two cells do not have any common ancestor (i.e., they
    // are from different faces).
    int GetCommonAncestorLevel(EGCellID other) const;

    // Iterator-style methods for traversing all the cells along the Hilbert
    // curve at a given level (across all 6 faces of the cube).  Note that the
    // end value is exclusive (just like standard STL iterators), and is not a
    // valid cell id.
    static EGCellID Begin(int level);
    static EGCellID End(int level);

    // Methods to encode and decode cell ids to compact text strings suitable
    // for display or indexing.  Cells at lower levels (i.e. larger cells) are
    // encoded into fewer characters.  The maximum token length is 16.
    //
    // Tokens preserve ordering, i.e. ToToken(x) < ToToken(y) iff x < y.
    //
    // ToToken() returns a string by value for convenience; the compiler
    // does this without intermediate copying in most cases.
    //
    // These methods guarantee that FromToken(ToToken(x)) == x even when
    // "x" is an invalid cell id.  All tokens are alphanumeric strings.
    // FromToken() returns EGCellID::None() for malformed inputs.
    string ToToken() const;
    static EGCellID FromToken(const char* token, size_t length);
    static EGCellID FromToken(const string& token);

    // Use encoder to generate a serialized representation of this cell id.
    // Can also encode an invalid cell.
    void Encode(Encoder* const encoder) const;

    // Decodes an EGCellID encoded by Encode(). Returns true on success.
    bool Decode(Decoder* const decoder);

    // Creates a human readable debug string.  Used for << and available for
    // direct usage as well.  The format is "f/dd..d" where "f" is a digit in
    // the range [0-5] representing the EGCellID face, and "dd..d" is a string
    // of digits in the range [0-3] representing each child's position with
    // respect to its parent.  (Note that the latter string may be empty.)
    //
    // For example "4/" represents EGCellID::FromFace(4), and "3/02" represents
    // EGCellID::FromFace(3).child(0).child(2).
    string ToString() const;

    // Converts a string in the format returned by ToString() to an EGCellID.
    // Returns EGCellID::None() if the string could not be parsed.
    //
    // The method name includes "Debug" in order to avoid possible confusion
    // with FromToken() above.
    static EGCellID FromDebugString(absl::string_view str);

    // Return the four cells that are adjacent across the cell's four edges.
    // Neighbors are returned in the order defined by EGCell::GetEdge.  All
    // neighbors are guaranteed to be distinct.
    void GetEdgeNeighbors(EGCellID neighbors[4]) const;

    // Return the neighbors of closest vertex to this cell at the given level,
    // by appending them to "output".  Normally there are four neighbors, but
    // the closest vertex may only have three neighbors if it is one of the 8
    // cube vertices.
    //
    // Requires: level < this->level(), so that we can determine which vertex is
    // closest (in particular, level == kMaxLevel is not allowed).
    void AppendVertexNeighbors(int level, std::vector<EGCellID>* output) const;

    // Append all neighbors of this cell at the given level to "output".  Two
    // cells X and Y are neighbors if their boundaries intersect but their
    // interiors do not.  In particular, two cells that intersect at a single
    // point are neighbors.  Note that for cells adjacent to a face vertex, the
    // same neighbor may be appended more than once.
    //
    // REQUIRES: nbr_level >= this->level().
    void AppendAllNeighbors(int nbr_level, std::vector<EGCellID>* output) const;

    /////////////////////////////////////////////////////////////////////
    // Low-level methods.

    // Return a leaf cell given its cube face (range 0..5) and
    // i- and j-coordinates (see EGcoords.h).
    static EGCellID FromFaceIJ(int face, int i, int j);

    // Return the (face, i, j) coordinates for the leaf cell corresponding to
    // this cell id.  Since cells are represented by the Hilbert curve position
    // at the center of the cell, the returned (i,j) for non-leaf cells will be
    // a leaf cell adjacent to the cell center.  If "orientation" is non-nullptr,
    // also return the Hilbert curve orientation for the current cell.
    int ToFaceIJOrientation(int* pi, int* pj, int* orientation) const;

    // Return the lowest-numbered bit that is on for this cell id, which is
    // equal to (uint64{1} << (2 * (kMaxLevel - level))).  So for example,
    // a.lsb() <= b.lsb() if and only if a.level() >= b.level(), but the
    // first test is more efficient.
    uint64 lsb() const { return id_ & (~id_ + 1); }

    // Return the lowest-numbered bit that is on for cells at the given level.
    static uint64 lsb_for_level(int level) {
        return uint64{1} << (2 * (kMaxLevel - level));
    }

    // Return the bound in (u,v)-space for the cell at the given level containing
    // the leaf cell with the given (i,j)-coordinates.
    static R2Rect IJLevelToBoundUV(int ij[2], int level);

    // When EGCellID is used as a key in one of the btree container types
    // (util/btree), indicate that linear rather than binary search should be
    // used.  This is much faster when the comparison function is cheap.
    typedef std::true_type goog_btree_prefer_linear_node_search;

private:
    // This is the offset required to wrap around from the beginning of the
    // Hilbert curve to the end or vice versa; see next_wrap() and prev_wrap().
    // SWIG doesn't understand uint64{}, so use static_cast.
    static const uint64 kWrapOffset = static_cast<uint64>(kNumFaces) << kPosBits;

    // Given a face and a point (i,j) where either i or j is outside the valid
    // range [0..kMaxSize-1], this function first determines which neighboring
    // face "contains" (i,j), and then returns the leaf cell on that face which
    // is adjacent to the given face and whose distance from (i,j) is minimal.
    static EGCellID FromFaceIJWrap(int face, int i, int j);

    // Inline helper function that calls FromFaceIJ if "same_face" is true,
    // or FromFaceIJWrap if "same_face" is false.
    static EGCellID FromFaceIJSame(int face, int i, int j, bool same_face);

    uint64 id_;
} ABSL_ATTRIBUTE_PACKED;  // Necessary so that structures containing EGCellID's
// can be ABSL_ATTRIBUTE_PACKED.

inline bool operator==(EGCellID x, EGCellID y) {
    return x.id() == y.id();
}

inline bool operator!=(EGCellID x, EGCellID y) {
    return x.id() != y.id();
}

inline bool operator<(EGCellID x, EGCellID y) {
    return x.id() < y.id();
}

inline bool operator>(EGCellID x, EGCellID y) {
    return x.id() > y.id();
}

inline bool operator<=(EGCellID x, EGCellID y) {
    return x.id() <= y.id();
}

inline bool operator>=(EGCellID x, EGCellID y) {
    return x.id() >= y.id();
}

inline EGCellID EGCellID::FromFace(int face) {
    return EGCellID((static_cast<uint64>(face) << kPosBits) + lsb_for_level(0));
}

inline EGCellID EGCellID::FromFacePosLevel(int face, uint64 pos, int level) {
    EGCellID cell((static_cast<uint64>(face) << kPosBits) + (pos | 1));
    return cell.parent(level);
}

inline int EGCellID::GetCenterSiTi(int* psi, int* pti) const {
    // First we compute the discrete (i,j) coordinates of a leaf cell contained
    // within the given cell.  Given that cells are represented by the Hilbert
    // curve position corresponding at their center, it turns out that the cell
    // returned by ToFaceIJOrientation is always one of two leaf cells closest
    // to the center of the cell (unless the given cell is a leaf cell itself,
    // in which case there is only one possibility).
    //
    // Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
    // jmin) be the coordinates of its lower left-hand corner, the leaf cell
    // returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
    // (imin + s/2 - 1, jmin + s/2 - 1).  The first case is the one we want.
    // We can distinguish these two cases by looking at the low bit of "i" or
    // "j".  In the second case the low bit is one, unless s == 2 (i.e. the
    // level just above leaf cells) in which case the low bit is zero.
    //
    // In the code below, the expression ((i ^ (int(id_) >> 2)) & 1) is true
    // if we are in the second case described above.
    int i, j;
    int face = ToFaceIJOrientation(&i, &j, nullptr);
    int delta = is_leaf() ? 1 : ((i ^ (static_cast<int>(id_) >> 2)) & 1) ? 2 : 0;

    // Note that (2 * {i,j} + delta) will never overflow a 32-bit integer.
    *psi = 2 * i + delta;
    *pti = 2 * j + delta;
    return face;
}

inline bool EGCellID::is_valid() const {
    return (face() < kNumFaces && (lsb() & 0x1555555555555555ULL));
}

inline int EGCellID::face() const {
    return id_ >> kPosBits;
}

inline uint64 EGCellID::pos() const {
    return id_ & (~uint64{0} >> kFaceBits);
}

inline int EGCellID::level() const {
    // We can't just EG_DCHECK(is_valid()) because we want level() to be
    // defined for end-iterators, i.e. EGCellID::End(kLevel).  However there is
    // no good way to define EGCellID::None().level(), so we do prohibit that.
    EG_DCHECK(id_ != 0);

    // A special case for leaf cells is not worthwhile.
    return kMaxLevel - (Bits::FindLSBSetNonZero64(id_) >> 1);
}

inline int EGCellID::GetSizeIJ() const {
    return GetSizeIJ(level());
}

inline double EGCellID::GetSizeST() const {
    return GetSizeST(level());
}

inline int EGCellID::GetSizeIJ(int level) {
    return 1 << (kMaxLevel - level);
}

inline double EGCellID::GetSizeST(int level) {
    return EG::IJtoSTMin(GetSizeIJ(level));
}

inline bool EGCellID::is_leaf() const {
    return int(id_) & 1;
}

inline bool EGCellID::is_face() const {
    return (id_ & (lsb_for_level(0) - 1)) == 0;
}

inline int EGCellID::child_position() const {
    // No need for a special implementation; the compiler optimizes this well.
    return child_position(level());
}

inline int EGCellID::child_position(int level) const {
    EG_DCHECK(is_valid());
    EG_DCHECK_GE(level, 1);
    EG_DCHECK_LE(level, this->level());
    return static_cast<int>(id_ >> (2 * (kMaxLevel - level) + 1)) & 3;
}

inline EGCellID EGCellID::range_min() const {
    return EGCellID(id_ - (lsb() - 1));
}

inline EGCellID EGCellID::range_max() const {
    return EGCellID(id_ + (lsb() - 1));
}

inline bool EGCellID::contains(EGCellID other) const {
    EG_DCHECK(is_valid());
    EG_DCHECK(other.is_valid());
    return other >= range_min() && other <= range_max();
}

inline bool EGCellID::intersects(EGCellID other) const {
    EG_DCHECK(is_valid());
    EG_DCHECK(other.is_valid());
    return other.range_min() <= range_max() && other.range_max() >= range_min();
}

inline EGCellID EGCellID::parent(int level) const {
    EG_DCHECK(is_valid());
    EG_DCHECK_GE(level, 0);
    EG_DCHECK_LE(level, this->level());
    uint64 new_lsb = lsb_for_level(level);
    return EGCellID((id_ & (~new_lsb + 1)) | new_lsb);
}

inline EGCellID EGCellID::parent() const {
    EG_DCHECK(is_valid());
    EG_DCHECK(!is_face());
    uint64 new_lsb = lsb() << 2;
    return EGCellID((id_ & (~new_lsb + 1)) | new_lsb);
}

inline EGCellID EGCellID::child(int position) const {
    EG_DCHECK(is_valid());
    EG_DCHECK(!is_leaf());
    // To change the level, we need to move the least-significant bit two
    // positions downward.  We do this by subtracting (4 * new_lsb) and adding
    // new_lsb.  Then to advance to the given child cell, we add
    // (2 * position * new_lsb).
    uint64 new_lsb = lsb() >> 2;
    return EGCellID(id_ + (2 * position + 1 - 4) * new_lsb);
}

inline EGCellID EGCellID::child_begin() const {
    EG_DCHECK(is_valid());
    EG_DCHECK(!is_leaf());
    uint64 old_lsb = lsb();
    return EGCellID(id_ - old_lsb + (old_lsb >> 2));
}

inline EGCellID EGCellID::child_begin(int level) const {
    EG_DCHECK(is_valid());
    EG_DCHECK_GE(level, this->level());
    EG_DCHECK_LE(level, kMaxLevel);
    return EGCellID(id_ - lsb() + lsb_for_level(level));
}

inline EGCellID EGCellID::child_end() const {
    EG_DCHECK(is_valid());
    EG_DCHECK(!is_leaf());
    uint64 old_lsb = lsb();
    return EGCellID(id_ + old_lsb + (old_lsb >> 2));
}

inline EGCellID EGCellID::child_end(int level) const {
    EG_DCHECK(is_valid());
    EG_DCHECK_GE(level, this->level());
    EG_DCHECK_LE(level, kMaxLevel);
    return EGCellID(id_ + lsb() + lsb_for_level(level));
}

inline EGCellID EGCellID::next() const {
    return EGCellID(id_ + (lsb() << 1));
}

inline EGCellID EGCellID::prev() const {
    return EGCellID(id_ - (lsb() << 1));
}

inline EGCellID EGCellID::next_wrap() const {
    EG_DCHECK(is_valid());
    EGCellID n = next();
    if (n.id_ < kWrapOffset) return n;
    return EGCellID(n.id_ - kWrapOffset);
}

inline EGCellID EGCellID::prev_wrap() const {
    EG_DCHECK(is_valid());
    EGCellID p = prev();
    if (p.id_ < kWrapOffset) return p;
    return EGCellID(p.id_ + kWrapOffset);
}

inline EGCellID EGCellID::Begin(int level) {
    return FromFace(0).child_begin(level);
}

inline EGCellID EGCellID::End(int level) {
    return FromFace(5).child_end(level);
}

std::ostream& operator<<(std::ostream& os, EGCellID id);

// Hasher for EGCellID.
// Example use: std::unordered_map<EGCellID, int, EGCellIDHash>.
struct EGCellIDHash {
    size_t operator()(EGCellID id) const {
        return std::hash<uint64>()(id.id());
    }
};


#undef IFNDEF_SWIG

#endif //EARTHGRID_EGCellID_H
