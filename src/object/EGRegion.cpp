//
// Created by pgl on 2018/12/6.
//

#include "object/EGRegion.h"

#include <vector>
#include "cellid/EGCellID.h"
#include "object/EGSphereCap.h"

void EGRegion::GetCellUnionBound(std::vector<EGCellId> *cell_ids) const {
    return GetCapBound().GetCellUnionBound(cell_ids);
}
