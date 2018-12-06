//
// Created by pgl on 2018/12/6.
//

#include "EGRegion.h"

#include "cellid/EGCellID.h"

#include <vector>

#include "s2cap.h"

void S2Region::GetCellUnionBound(std::vector<EGCellId> *cell_ids) const {
    return GetCapBound().GetCellUnionBound(cell_ids);
}
