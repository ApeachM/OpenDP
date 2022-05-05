#include "circuit.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef bg::model::polygon<bg::model::d2::point_xy<double>> polygon;

// RTrees
typedef bgi::rtree<std::pair<box, cell *>, bgi::quadratic<6> > cellboxRtree;
typedef bgi::rtree<std::pair<point, cell *>, bgi::quadratic<6> > cellpointRtree;
typedef bgi::rtree<std::pair<box, field *>, bgi::quadratic<6> > fieldboxRtree;

typedef std::pair<box, cell *> boxValue;
typedef std::pair<point, cell *> ptValue;

// Update Fields coordinate & area
void circuit::field_box_init() {
    double unit = field_unit;

    double fieldUnit = unit * rowHeight;
    int x_fieldNum = (int) ceil((rx - lx) / fieldUnit);
    int y_fieldNum = (int) ceil((ty - by) / fieldUnit);
    int numField = x_fieldNum * y_fieldNum;
    Fields.resize(numField);

    int id = 0;
    for (int y = 0; y < y_fieldNum; y++) {
        for (int x = 0; x < x_fieldNum; x++) {
            unsigned fieldId = y * x_fieldNum + x;
            Fields[fieldId].id = fieldId;
            Fields[fieldId].xLL = lx + x * fieldUnit;
            Fields[fieldId].yLL = by + y * fieldUnit;
            Fields[fieldId].xUR = min(Fields[fieldId].xLL + fieldUnit, rx);
            Fields[fieldId].yUR = min(Fields[fieldId].yLL + fieldUnit, ty);
            Fields[fieldId].fieldArea =
                    1.0 * (Fields[fieldId].xUR - Fields[fieldId].xLL) * (Fields[fieldId].yUR - Fields[fieldId].yLL);
            if (Fields[fieldId].fieldArea == 0)
                cout << "[WARNING] Zero init-Area of Field-" << fieldId << endl;
            //else if(Fields[fieldId].fieldArea != fieldUnit*fieldUnit)
            //    cout << "[WARNING] Non-fieldUnit^2 init-Area of Field-" << fieldId << endl;

        }
    }
}

bool SortUpUtil(const field &fa, const field &fb) {
    return (fa.fieldCells.size() > fb.fieldCells.size());
}

bool SortId(const cell *ca, const cell *cb) {
    return (ca->id < cb->id);
}

// Update Gcells coordinate & area
void circuit::gcell_box_init() {
    double unit = gcell_unit;
    double gcellUnit = unit * rowHeight; // default unit=50.0
    int x_gcellNum = (int) ceil((rx - lx) / gcellUnit);
    int y_gcellNum = (int) ceil((ty - by) / gcellUnit);
    int numGcells = x_gcellNum * y_gcellNum;
    Gcells.resize(numGcells);
    //GcAvgDisp.resize(numGcells);
    //GcAvgDisp.shrink_to_fit();

    for (int y = 0; y < y_gcellNum; y++) {
        for (int x = 0; x < x_gcellNum; x++) {
            unsigned gcellId = y * x_gcellNum + x;
            Gcells[gcellId].pid = gcellId;
            Gcells[gcellId].xLL = lx + x * gcellUnit;
            Gcells[gcellId].yLL = by + y * gcellUnit;
            Gcells[gcellId].xUR = min(Gcells[gcellId].xLL + gcellUnit, rx);
            Gcells[gcellId].yUR = min(Gcells[gcellId].yLL + gcellUnit, ty);
            Gcells[gcellId].fieldArea =
                    1.0 * (Gcells[gcellId].xUR - Gcells[gcellId].xLL) * (Gcells[gcellId].yUR - Gcells[gcellId].yLL);
            if (Gcells[gcellId].fieldArea == 0)
                cout << "[WARNING] Zero init-Area of Gcell-" << gcellId << endl;
            //else if(Gcells[gcellId].fieldArea != gcellUnit*gcellUnit)
            //    cout << "[WARNING] Non-gcellUnit^2 init-Area of Gcell-" << gcellId << endl;
        }
    }
}

void circuit::gcell_post_init() {
    // Get belonging cells for all Gcell
    for (int i = 0; i < (int) Gcells.size(); i++) {
        get_gcell_cells(i);
    }
    // Sort Gcells in order of utilization
    sort(Gcells.begin(), Gcells.end(), SortUpUtil);
    for (int i = 0; i < (int) Gcells.size(); i++) {
        Gcells[i].id = i;
    }
}

// Create Rtree cell-cell, field-gcell
void circuit::Rtree_init() {
    if (btw_cell_rtree != nullptr && cell_fg_rtree != nullptr) {
        cout << " Cell-Cell Rtree & Cell-Field/Gcell Rtree already exist!" << endl;
        return;
    }
    cout << " Generate Cell-Cell Rtree & Cell-Field/Gcell Rtree" << endl;

    // cell box -> find cell box
    btw_cell_rtree = (void *) (new cellboxRtree);
    auto *ccRtree = (cellboxRtree *) btw_cell_rtree;

    // field box -> find cell point
    cell_fg_rtree = (void *) (new cellpointRtree);
    auto *cfgRtree = (cellpointRtree *) cell_fg_rtree;

//    // gcell box -> find field box
//    field_gcell_rtree = (void*) (new fieldboxRtree);
//    fieldboxRtree* fgRtree = (fieldboxRtree*) field_gcell_rtree;
//
//    // cell point -> find field box
//    field_rtree = (void*) (new fieldboxRtree);
//    fieldboxRtree* fRtree = (fieldboxRtree*) field_rtree;

    // For all cells -> +-0.1 is trick for bgi::intersects
    for (int i = 0; i < (int) cells.size(); i++) {
        cell *theCell = &cells[i];
        box b(point((double)theCell->x_coord + 0.1, (double) theCell->y_coord + 0.1),
        point((double) (theCell->x_coord + theCell->width - 0.1),
              (double) (theCell->y_coord + theCell->height - 0.1)));
        point pt((double) theCell->x_coord + 0.1, (double) theCell->y_coord + 0.1);
        ccRtree->insert(make_pair(b, theCell));
        cfgRtree->insert(make_pair(pt, theCell));
    }
}


void circuit::Rtree_init_gcells_n_fields() {
    if (field_gcell_rtree != nullptr && field_rtree != nullptr) {
        cout << " Field Rtree & Field-Gcell Rtree already exist!" << endl;
        return;
    }

    cout << " Generate Field Rtree & Field-Gcell Rtree" << endl;

    // gcell box -> find field box
    field_gcell_rtree = (void *) (new fieldboxRtree);
    auto *fgRtree = (fieldboxRtree *) field_gcell_rtree;

    // cell point -> find field box
    field_rtree = (void *) (new fieldboxRtree);
    auto *fRtree = (fieldboxRtree *) field_rtree;

    for (auto & Field : Fields) {
        field *theField = &Field;
        box bfg(point(theField->xLL + 0.1, theField->yLL + 0.1),
                point(theField->xUR - 0.1, theField->yUR - 0.1));
        box bf(point(theField->xLL, theField->yLL), point(theField->xUR - 0.1, theField->yUR - 0.1));
        fgRtree->insert(make_pair(bfg, theField));
        fRtree->insert(make_pair(bf, theField));
    }
}

void circuit::get_fields_in_gcells() {
    auto *fgRtree = (fieldboxRtree *) field_gcell_rtree;

    for (auto & Gcell : Gcells) {
        field *theGcell = &Gcell;
        vector<pair<box, field *>> fields;
        /// Use fgRtree
        box queryBox(point(theGcell->xLL, theGcell->yLL),
                     point(theGcell->xUR, theGcell->yUR));
        fgRtree->query(bgi::intersects(queryBox), back_inserter(fields));
        //theGcell->subFields.reserve(fields.size());

        for (auto field: fields) {
            theGcell->subFields.push_back(field.second);
        }
    }
}

void circuit::Rtree_update(cell *tarCell, int moveType) {
    cellboxRtree *ccRtree = (cellboxRtree *) btw_cell_rtree;
    cellpointRtree *cfgRtree = (cellpointRtree *) cell_fg_rtree;

    double xLL = max(lx, (double) tarCell->init_x_coord) + 0.1;
    double yLL = max(by, (double) tarCell->init_y_coord) + 0.1;
    double xUR = tarCell->init_x_coord + tarCell->width - 0.1;
    double yUR = tarCell->init_y_coord + tarCell->height - 0.1;

    double xLL2 = tarCell->x_coord + 0.1;
    double yLL2 = tarCell->y_coord + 0.1;
    double xUR2 = tarCell->x_coord + tarCell->width - 0.1;
    double yUR2 = tarCell->y_coord + tarCell->height - 0.1;

    point init_pt(xLL, yLL);
    point pt(xLL2, yLL2);
    box init_b(point(xLL, yLL), point(xUR, yUR));
    box b(point(xLL2, yLL2), point(xUR2, yUR2));

    ptValue init_ptv = make_pair(init_pt, tarCell);
    ptValue ptv = make_pair(pt, tarCell);
    boxValue init_bv = make_pair(init_b, tarCell);
    boxValue bv = make_pair(b, tarCell);

    if (moveType == 1) { // map_move
        ccRtree->remove(init_bv);
        ccRtree->insert(bv);
        cfgRtree->remove(init_ptv);
        cfgRtree->insert(ptv);
    } else if (moveType == 2) { // shift_move -> get cell lists
        cout << "[WARNING] Shift move Rtree update is not completed!" << endl;
        ///////////////////////////// Need implementation.. ///
    }
}

void circuit::Rtree_clear() {
    if (btw_cell_rtree)
        delete (cellboxRtree *) btw_cell_rtree;
    if (cell_fg_rtree)
        delete (cellpointRtree *) cell_fg_rtree;
    //if(field_rtree)
    //    delete (fieldboxRtree*) field_rtree;
    //if(field_gcell_rtree)
    //    delete (fieldboxRtree*) field_gcell_rtree;
    btw_cell_rtree = nullptr;
    cell_fg_rtree = nullptr;
    //field_rtree       = nullptr;
    //field_gcell_rtree = nullptr;
}


// Get cell overlaps (cell->ovcells)
void circuit::get_ovcells(cell *theCell) {
    theCell->ovcells.clear();
    //cell* theCell = &cells[cellId];
    cellboxRtree *ccRtree = (cellboxRtree *) btw_cell_rtree;
    vector<pair<box, cell *>> ovcells;

//    box queryBox(point((double)theCell->x_coord, (double) theCell->y_coord),
//    point((double) (theCell->x_coord + theCell->width),
//          (double) (theCell->y_coord, theCell->y_coord + theCell->height)));
    box queryBox(point((double)theCell->x_coord, (double) theCell->y_coord),
    point((double) (theCell->x_coord + theCell->width),
          (double) (theCell->y_coord + theCell->height)));
    ccRtree->query(bgi::intersects(queryBox), back_inserter(ovcells));
    //theCell->ovcells.reserve(ovcells.size());

    for (auto cell: ovcells) {
        if (cell.second == theCell) continue;
        theCell->ovcells.push_back(cell.second);
    }
    theCell->overlapNum = (int) theCell->ovcells.size();
}

// Get fieldId for the cell
int circuit::get_field_id(cell *theCell) {
    fieldboxRtree *fRtree = (fieldboxRtree *) field_rtree;
    vector<pair<box, field *>> fields;

    point queryPt((double) theCell->x_coord, (double) theCell->y_coord);
    fRtree->query(bgi::intersects(queryPt), back_inserter(fields));

    if (fields.size() != 1) {
        cout << "[ERROR] the cell (id: " << theCell->id << ") belongs to more than one fields!" << endl;
        exit(0);
    }

    //if(theCell->fieldId != (fields[0].second)->id) {
    //    cout << " Cell-" << theCell->id << "(" << theCell->name << ")" << " move field ("
    //         << theCell->fieldId << "->" << (fields[0].second)->id << ")" << endl;
    //}
    return (fields[0].second)->id;
}


// Get field cells (field->fieldCells) & cell-overlap-num (fieldOverlap)
void circuit::get_field_cells_n_overlapnum(int fieldId) {
    double cellArea = 0.0;
    double placedArea = 0.0;
    int fov = 0;

    vector<pair<point, cell *>> fieldcells;
    cellpointRtree *cfgRtree = (cellpointRtree *) cell_fg_rtree;
    field *theField = &Fields[fieldId];

    box queryBox(point(theField->xLL, theField->yLL), point(theField->xUR, theField->yUR));
    cfgRtree->query(bgi::intersects(queryBox), back_inserter(fieldcells));
    //theField->fieldCells.reserve(fieldcells.size());
    for (auto cell: fieldcells) {
        // skip for FR-cell or block-cell (what about isPlaced cell?)
        if ((cell.second)->inGroup || (cell.second)->isFixed || macros[(cell.second)->type].type == "BLOCK") continue;
        (cell.second)->fieldId = fieldId;
        theField->fieldCells.push_back(cell.second);
        cellArea += (cell.second)->width * (cell.second)->height;
        if ((cell.second)->isPlaced)
            placedArea += (cell.second)->width * (cell.second)->height;
        fov += (cell.second)->ovcells.size();   // ovcells were updated in get_gcell_cells
    }
    theField->cellArea = cellArea;
    theField->placedArea = placedArea;
    theField->fieldOverlap = fov;
}

void circuit::field_overlapnum_update(int fieldId) {
    int fov = 0;
    double placedArea = 0.0;
    field *theField = &Fields[fieldId];

    for (auto cell: theField->fieldCells) {
        if (cell->isPlaced)
            placedArea += cell->width * cell->height;
        get_ovcells(cell);
        fov += cell->ovcells.size();
    }
    theField->placedArea = placedArea;
    theField->fieldOverlap = fov;
}

// Get field cells (field->fieldCells) & cell-overlap-num (fieldOverlap)
void circuit::get_gcell_cells(int gcellId) {
    vector<pair<point, cell *>> gcellcells;
    auto *cfgRtree = (cellpointRtree *) cell_fg_rtree;
    field *theGcell = &Gcells[gcellId];

    box queryBox(point(theGcell->xLL, theGcell->yLL), point(theGcell->xUR, theGcell->yUR));
    cfgRtree->query(bgi::intersects(queryBox), back_inserter(gcellcells));

    theGcell->fieldCells.reserve(gcellcells.size());
    for (auto cell: gcellcells) {
        cell.second->gcellId = gcellId;
        theGcell->fieldCells.push_back(cell.second);
        get_ovcells(cell.second);
    }
    sort(theGcell->fieldCells.begin(), theGcell->fieldCells.end(), SortId);
}

//pair<double, double> circuit::field_intersect_area_with_fr_or_block(int fieldId, string mode)
void circuit::field_area_n_block_area(int fieldId) {
    field *theField = &Fields[fieldId];
    double frArea = 0.0;
    double blockArea = 0.0;
    box fieldBox(point(theField->xLL, theField->yLL), point(theField->xUR, theField->yUR));
    polygon fieldPoly;
    bg::assign(fieldPoly, fieldBox);

    // FR area
    for (int i = 0; i < (int) groups.size(); i++) {
        group *theGroup = &groups[i];
        for (int j = 0; j < (int) theGroup->regions.size(); j++) {
            deque<polygon> frOutput;
            polygon poly;
            rect *R = &theGroup->regions[j];
            box b(point(R->xLL, R->yLL), point(R->xUR, R->yUR));
            bg::assign(poly, b);
            bg::intersection(fieldPoly, poly, frOutput);
            int num = 0;

            BOOST_FOREACH(polygon const &p, frOutput) {
                            frArea += bg::area(p);
                        }
            if (frOutput.size() > 1) {
                cout << "[WARNING] Multiple intersection exists with a rect-region and field-area!" << endl;
                cout << "frArea: " << frArea << endl;
            }
        }
    }

    // Block area
    for (int i = 0; i < (int) cells.size(); i++) {
        cell *theCell = &cells[i];
        if (macros[theCell->type].type == "BLOCK") {
            deque<polygon> blkOutput;
            polygon poly;
            box b(point((double)theCell->x_coord, (double) theCell->y_coord), point(
                    (double) (theCell->x_coord + theCell->width), (double) (theCell->y_coord + theCell->height)));
            bg::assign(poly, b);
            bg::intersection(fieldPoly, poly, blkOutput);
            int num = 0;
            BOOST_FOREACH(polygon const &p, blkOutput) {
                            blockArea += bg::area(p);
                        }
        }
    }

    //theField->fieldArea -= frArea;
    //theField->fieldArea -= (frArea + blockArea);
    theField->fieldArea = 1.0 * (theField->xUR - theField->xLL) * (theField->yUR - theField->yLL)
                          - (frArea + blockArea);
    //theField->blockArea = blockArea;
}


// not used
void circuit::get_field_density(int fieldId, string mode) {
    field *theField = &Fields[fieldId];
    double cellArea = 0.0;

    if (mode == "init") {
        field_area_n_block_area(fieldId);
    }

    for (int i = 0; i < (int) theField->fieldCells.size(); i++) {
        cell *theCell = theField->fieldCells[i];
        cellArea += (theCell->width) * (theCell->height);
    }

    //theField->fieldDensity =  1.0 * (cellArea + theField->blockArea) / theField->fieldArea;
}




