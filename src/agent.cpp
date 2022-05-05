#include "circuit.h"

#define FEATURE_NUM 5

using namespace std;

State::State() {}

State::State(vector<cell *> cells) { cell_list = cells; }

State::~State() {}

Agent::Agent() : moveFailCnt(0), targetGcell(-1), Pov(0) {}

Agent::~Agent() {}

vector<cell *> circuit::non_group_movable_cells(int gcellid) {
    vector<cell *> cell_list;
    field *theGcell = &Gcells[gcellid];

    for (int i = 0; i < theGcell->fieldCells.size(); i++) {
        cell *theCell = theGcell->fieldCells[i];
        if (theCell->isFixed || theCell->isPlaced || theCell->inGroup) continue;
        cell_list.push_back(theCell);
    }
    return cell_list;
}

void Agent::state_init(circuit *ck, vector<cell *> cells, int gcellid) {
    targetGcell = gcellid;
    state.cell_list = cells;

    for (int i = 0; i < (int) state.cell_list.size(); i++) {
        cell *tarCell = state.cell_list[i];
        tarCell->sindex = i;
    }
    cout << " Cell list is initialized!" << endl;

    feature_init(ck);
}

void Agent::feature_init(circuit *ck) {

    for (auto theCell: state.cell_list) {
        /* id, moveTry, width, height is initialized or unchanged */
        ck->get_ovcells(theCell); // overlapNum & ovcells init.
    }

    for (auto theField: ck->Gcells[targetGcell].subFields) {
        ck->get_field_cells_n_overlapnum(theField->id);
        ck->field_area_n_block_area(theField->id);
    }
    cout << " Feature of cell list is initalized!" << endl;

    // Update for state.features
    for (int i = 0; i < (int) state.cell_list.size(); i++) {
        cell *theCell = state.cell_list[i];

        vector<double> cellFeature;
        cellFeature.push_back(theCell->id);
        cellFeature.push_back(theCell->moveTry);
        cellFeature.push_back(theCell->width);
        cellFeature.push_back(theCell->height);
        cellFeature.push_back(theCell->min2SrchDist);
        cellFeature.push_back(theCell->overlapNum);
        cellFeature.push_back(ck->Fields[theCell->fieldId].cellArea);
        cellFeature.push_back(ck->Fields[theCell->fieldId].fieldArea);
        cellFeature.push_back(ck->Fields[theCell->fieldId].fieldOverlap);
        cellFeature.push_back(ck->Fields[theCell->fieldId].placedArea);

        state.features.push_back(cellFeature);
    }
    cout << "Feature 2D-vector is initialized!" << endl;
}

bool Agent::is_done() {
    if (state.cell_list.size() == 0)
        cout << "[WARNING] No cell list to place" << endl;

    for (auto theCell: state.cell_list) {
        if (!theCell->moveTry) return false;
    }
    return true;
}

int Agent::action(circuit *ck, int tarID) {
    int moveType = 1;
    ck->REWARD = 0.0;
    tarCell = &ck->cells[tarID];

    if (!tarCell->isPlaced) {    // shift_moved cells are only process feature_update
        if (ck->map_move(tarCell, "init_coord") == false) {
            moveType = 2;
            cout << " Map move fail! Trying shift move.." << endl;
            if (ck->shift_move(tarCell, "init_coord") == false) {
                moveType = 0;
                cout << "[WARNING] " << tarCell->name << " -> Move fail!" << endl;
            }
        }
        tarCell->moveTry = true;
    }
    if (moveType == 0) moveFailCnt++;

    return moveType;
}

void Agent::feature_update(circuit *ck, int tarID, int moveType) {
    if (moveType == 0) return;
    tarCell = &ck->cells[tarID];
    ck->Rtree_update(tarCell, moveType);

    // For overlapNum
    ck->get_ovcells(tarCell);
    for (auto cell: tarCell->ovcells) { ck->get_ovcells(cell); }

    // For fieldId
    int fId1 = tarCell->fieldId;
    int fId2 = ck->get_field_id(tarCell);
    if (fId1 == fId2) {
        // Update fieldOverlap only.
        ck->field_overlapnum_update(fId1);
    } else {
        // For both Fields...
        // Update fieldCells and fieldOverlap
        ck->get_field_cells_n_overlapnum(fId1);
        ck->get_field_cells_n_overlapnum(fId2);
        tarCell->fieldId = fId2;
    }

    // Update for state.features
    // return only effected cells.. cellid -> effected cell list -> update state.features -> returning effected features
    effCells.clear();
    vector<int>().swap(effCells);

    // 1. tarCell
    effCells.push_back(tarCell->sindex);
    // 2. Ovcells
    for (int i = 0; i < (int) tarCell->ovcells.size(); i++) {
        if (tarCell->ovcells[i]->sindex == -1) continue;
        if (tarCell->ovcells[i]->moveTry == true) continue;
        effCells.push_back(tarCell->ovcells[i]->sindex);
    }
    // 3. FieldCells
    for (int i = 0; i < (int) ck->Fields[fId1].fieldCells.size(); i++) {
        if (ck->Fields[fId1].fieldCells[i]->sindex == -1) continue;
        if (ck->Fields[fId1].fieldCells[i]->moveTry == true) continue;
        effCells.push_back(ck->Fields[fId1].fieldCells[i]->sindex);
    }
    if (fId1 != fId2) {
        for (int i = 0; i < (int) ck->Fields[fId2].fieldCells.size(); i++) {
            if (ck->Fields[fId2].fieldCells[i]->sindex == -1) continue;
            if (ck->Fields[fId2].fieldCells[i]->moveTry == true) continue;
            effCells.push_back(ck->Fields[fId2].fieldCells[i]->sindex);
        }
    }
    // 4. shift_move
    if (ck->overlap_region_cells.size()) {
        assert(moveType == 2);
        for (int i = 0; i < ck->overlap_region_cells.size(); i++)
            effCells.push_back(ck->overlap_region_cells[i]->id);
        ck->overlap_region_cells.clear();
        vector<cell *>().swap(ck->overlap_region_cells);
    }

    sort(effCells.begin(), effCells.end());
    effCells.erase(unique(effCells.begin(), effCells.end()), effCells.end());

    for (int i = 0; i < (int) effCells.size(); i++) {
        cell * updateCell = state.cell_list[effCells[i]];
        // 0: id
        state.features[effCells[i]][1] = (double) updateCell->moveTry;
        // 2: width, 3: height, 4: minsrchDist
        state.features[effCells[i]][5] = (double) updateCell->overlapNum;
        state.features[effCells[i]][6] = (double) ck->Fields[updateCell->fieldId].cellArea;
        state.features[effCells[i]][7] = (double) ck->Fields[updateCell->fieldId].fieldArea;
        state.features[effCells[i]][8] = (double) ck->Fields[updateCell->fieldId].fieldOverlap;
        state.features[effCells[i]][9] = (double) ck->Fields[updateCell->fieldId].placedArea;
    }
}
