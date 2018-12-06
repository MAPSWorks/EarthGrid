//
// Created by pgl on 2018/12/6.
//

#ifndef EARTHGRID_EGREGION_H
#define EARTHGRID_EGREGION_H

#include <vector>

class EGSphereCap;
class EGCell;
class EGCellId;
class EGLatLongRect;

// An EGRegion represents a two-dimensional object over the unit sphere.
// It is an abstract interface with various concrete subtypes.
// EGRegion是一个用于表球面二维对象的抽象类，定义了各种子类的基本接口。
// 该接口的主要作用是将复杂的区域对象转换为更为简单的区域对象，因此仅封装了和估计有关的接口。
// 其他接口的实现则交由子类去实现。
// The main purpose of this interface is to allow complex regions to be
// approximated as simpler regions.  So rather than having a wide variety
// of virtual methods that are implemented by all subtypes, the interface
// is restricted to methods that are useful for computing approximations.

class EGRegion {
public:
    virtual ~EGRegion() {}

    // Returns a deep copy of the object.
    //
    // Note that each subtype of EGRegion returns a pointer to an object of its
    // own type (e.g., EGCap::Clone() returns an EGCap*).
    virtual EGRegion* Clone() const = 0;

    // Returns a bounding spherical cap that contains the object.  The bound may
    // not be tight.
    // 返回区域对象的球面外包球冠。
    virtual EGSphereCap GetCapBound() const = 0;

    // Returns a bounding latitude-longitude rectangle that contains the object.
    // The bound may not be tight.
    // 返回区域对象的经纬度外包矩形
    virtual EGLatLongRect GetRectBound() const = 0;

    // Returns a small collection of EGCellIds whose union covers the object.
    // The cells are not sorted, may have redundancies (such as cells that
    // contain other cells), and may cover much more area than necessary.
    //
    // This method is not intended for direct use by client code.  Clients
    // should typically use EGRegionCoverer::GetCovering, which has options to
    // control the size and accuracy of the covering.  Alternatively, if you
    // want a fast covering and don't care about accuracy, consider calling
    // EGRegionCoverer::GetFastCovering (which returns a cleaned-up version of
    // the covering computed by this method).
    //
    // GetCellUnionBound() implementations should attempt to return a small
    // covering (ideally 4 cells or fewer) that covers the object and can be
    // computed quickly.  The result is used by EGRegionCoverer as a starting
    // point for further refinement.
    //
    // TODO(ericv): Remove the default implementation.
    // 返回区域对象的最小外包格元列表（4个或者更少）
    virtual void GetCellUnionBound(std::vector<EGCellId> *cell_ids) const;

    // Returns true if the object completely contains the given cell, otherwise
    // returns false.
    virtual bool Contains(const EGCell& cell) const = 0;

    // If this method returns false, the object does not intersect the given
    // cell.  Otherwise, either object intersects the cell, or the intersection
    // relationship could not be determined.
    //
    // Note that there is currently exactly one implementation of this method
    // (EGLatLngRect::MayIntersect) that takes advantage of the semantics above
    // to be more efficient.  For all other EGRegion subtypes, this method
    // returns true if the object intersect the cell and false otherwise.
    virtual bool MayIntersect(const EGCell& cell) const = 0;

    // Returns true if and only if the given point is contained by the object.
    // The point 'p' is generally required to be unit length, although some
    // subtypes may relax this restriction.
    virtual bool Contains(const EGPoint& p) const = 0;

    //////////////////////////////////////////////////////////////////////////
    // Many EGRegion subtypes also define the following non-virtual methods.
    //////////////////////////////////////////////////////////////////////////

    // Appends a serialized representation of the object to "encoder".
    //
    // The representation chosen is left up to the sub-classes but it should
    // satisfy the following constraints:
    // - It should encode a version number.
    // - It should be deserializable using the corresponding Decode method.
    // - Performance, not space, should be the chief consideration. Encode() and
    //   Decode() should be implemented such that the combination is equivalent
    //   to calling Clone().
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    //
    // void Encode(Encoder* const encoder) const;

    // Decodes an EGRegion encoded with Encode().  Note that this method
    // requires that an EGRegion object of the appropriate concrete type has
    // already been constructed.  It is not possible to decode regions of
    // unknown type.
    //
    // Whenever the Decode method is changed to deal with new serialized
    // representations, it should be done so in a manner that allows for
    // older versions to be decoded i.e. the version number in the serialized
    // representation should be used to decide how to decode the data.
    //
    // Returns true on success.
    //
    // bool Decode(Decoder* const decoder);

    // Provides the same functionality as Decode, except that decoded regions
    // are allowed to point directly into the Decoder's memory buffer rather
    // than copying the data.  This method can be much faster for regions that
    // have a lot of data (such as polygons), but the decoded object is only
    // valid within the scope (lifetime) of the Decoder's memory buffer.
    //
    // bool DecodeWithinScope(Decoder* const decoder);
};


#endif //EARTHGRID_EGREGION_H
