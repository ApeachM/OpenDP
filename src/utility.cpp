///////////////////////////////////////////////////////////////////////////////
//// Authors: SangGi Do(sanggido@unist.ac.kr), Mingyu Woo(mwoo@eng.ucsd.edu)
////          (respective Ph.D. advisors: Seokhyeong Kang, Andrew B. Kahng)
////
////          Original parsing structure was made by Myung-Chul Kim (IBM).
////
//// BSD 3-Clause License
////
//// Copyright (c) 2018, SangGi Do and Mingyu Woo
//// All rights reserved.
////
//// Redistribution and use in source and binary forms, with or without
//// modification, are permitted provided that the following conditions are met:
////
//// * Redistributions of source code must retain the above copyright notice, this
////   list of conditions and the following disclaimer.
////
//// * Redistributions in binary form must reproduce the above copyright notice,
////   this list of conditions and the following disclaimer in the documentation
////   and/or other materials provided with the distribution.
////
//// * Neither the name of the copyright holder nor the names of its
////   contributors may be used to endorse or promote products derived from
////   this software without specific prior written permission.
////
//// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/////////////////////////////////////////////////////////////////////////////////

#include "circuit.h"

#define _DEBUG
#define SOFT_IGNORE true

void circuit::power_mapping() {
    for (int i = 0; i < rows.size(); i++) {
        row *theRow = &rows[i];
        if (initial_power == VDD) {
            if (i % 2 == 0)
                theRow->top_power = VDD;
            else
                theRow->top_power = VSS;
        } else {
            if (i % 2 == 0)
                theRow->top_power = VSS;
            else
                theRow->top_power = VDD;
        }
    }
    return;
}

void circuit::evaluation() {

    double avg_displacement = 0;
    double sum_displacement = 0;
    double max_displacement = 0;
    int count_displacement = 0;
    double violated_const = 0;

    int H1 = 0;
    int H2 = 0;
    int H3 = 0;
    int H4 = 0;
    int h1 = 0;
    int h2 = 0;
    int h3 = 0;
    int h4 = 0;
    double disp_H1 = 0;
    double disp_H2 = 0;
    double disp_H3 = 0;
    double disp_H4 = 0;


    cell *maxCell = NULL;
    double displacement;
    for (int i = 0; i < cells.size(); i++) {
        cell *theCell = &cells[i];
        displacement = abs(theCell->init_x_coord - theCell->x_coord) + abs(theCell->init_y_coord - theCell->y_coord);
        sum_displacement += displacement;
        if (displacement > max_displacement) {
            max_displacement = displacement;
            maxCell = theCell;
        }
        count_displacement++;
        if ((displacement / rowHeight) > max_disp_const) {
            violated_const += static_cast<double>(displacement / rowHeight);
        }

        if (theCell->height == 2000) {
            H1++;
            disp_H1 += displacement;
            h1 = 1;
        } else if (theCell->height == 4000) {
            H2++;
            disp_H2 += displacement;
            h2 = 1;
        } else if (theCell->height == 6000) {
            H3++;
            disp_H3 += displacement;
            h3 = 1;
        } else if (theCell->height == 8000) {
            H4++;
            disp_H4 += displacement;
            h4 = 1;
        }
        //int num_group_cells = 0;
        //int num_nongroup_cells = 0;
        //if(displacement > max_disp_const){
//		if(theCell->isGroup)
//			num_group_cells++;
//		else
//			num_nongroup_cells++;
//	}
    }
//	cout<<"cell anlaysis: " <<num_group_cells << endl;
//	cout<<"cell analysis2: " << num_nongroup_cells << endl;
    avg_displacement = sum_displacement / count_displacement;
    double sam;

    double avg_disp_H1 = 0;
    double avg_disp_H2 = 0;
    double avg_disp_H3 = 0;
    double avg_disp_H4 = 0;
    if (h1) avg_disp_H1 = static_cast<double>(disp_H1 / H1);
    if (h2) avg_disp_H2 = static_cast<double>(disp_H2 / H2);
    if (h3) avg_disp_H3 = static_cast<double>(disp_H3 / H3);
    if (h4) avg_disp_H4 = static_cast<double>(disp_H4 / H4);
    sam = static_cast<double>((avg_disp_H1 + avg_disp_H2 + avg_disp_H3 + avg_disp_H4) / (h1 + h2 + h3 + h4)) /
          rowHeight; // should modify H1 H2 H3 H4 count


    double max_util = 0.0;
    for (int i = 0; i < groups.size(); i++) {
        group *theGroup = &groups[i];
        if (max_util < theGroup->util)
            max_util = theGroup->util;
    }

    double smm;
    if ((violated_const / max_disp_const) > 1)
        smm = 1 + static_cast<double>((violated_const / max_disp_const) *
                                      (max_displacement / (100 * rowHeight)));
    else smm = 1 + static_cast<double>(max_displacement / (100 * rowHeight));

    //double shpwl = std::max((HPWL("")-HPWL("INIT"))/HPWL("INIT"), 0.0) * (1 + std::max(calc_density_factor(8.0), 0.2));
    double shpwl = std::max((HPWL("") - HPWL("INIT")) / HPWL("INIT"), 0.0) * (1 + 0.2);

    cout << " - - - - - EVALUATION - - - - - " << endl;
    cout << " [EPI] AVG_displacement : " << avg_displacement << endl;
    cout << " [EPI] SUM_displacement : " << sum_displacement << endl;
    cout << " [EPI] MAX_displacement : " << max_displacement << endl;
    cout << " - - - - - - - - - - - - - - - - " << endl;
    cout << " [EPI] GP HPWL          : " << HPWL("INIT") / 1000000 << endl;
    cout << " [EPI] HPWL             : " << HPWL("") / 1000000 << endl;
    cout << " [EPI] avg_Disp_site    : " << Disp() / cells.size() / wsite << endl;
    cout << " [EPI] avg_Disp_row     : " << Disp() / cells.size() / rowHeight << endl;
    cout << " [EPI] delta_HPWL       : " << (HPWL("") - HPWL("INIT")) / HPWL("INIT") * 100 << endl;
    cout << " [EPI] Smm              : " << smm << endl;
    cout << " [EPI] Sam              : " << sam << endl;
    cout << " [EPI] Shpwl            : " << shpwl << endl;
    cout << " [EPI] Stotal           : " << smm * sam * (1 + shpwl) << endl;
    //   cout << " Violated cells   : " << H1 << "  " << H2 << "  " << H3 << "  " << H4 << endl;
    return;
}

double circuit::Disp() {
    double result = 0.0;
    for (int i = 0; i < cells.size(); i++) {
        cell *theCell = &cells[i];
        if (theCell->x_coord == 0 && theCell->y_coord == 0)
            continue;
        result += abs(theCell->init_x_coord - theCell->x_coord) + abs(theCell->init_y_coord - theCell->y_coord);
    }
    return result;
}

double circuit::HPWL(string mode) {
    double hpwl = 0;

    double x_coord = 0;
    double y_coord = 0;

    for (int i = 0; i < nets.size(); i++) {
        rect box;
        net *theNet = &nets[i];
        //cout << " net name : " << theNet->name << endl;
        pin *source = &pins[theNet->source];

        if (source->type == NONPIO_PIN) {
            cell *theCell = &cells[source->owner];
            if (mode == "INIT") {
                x_coord = theCell->init_x_coord;
                y_coord = theCell->init_y_coord;
            } else {
                x_coord = theCell->x_coord;
                y_coord = theCell->y_coord;
            }
            box.xLL = box.xUR = x_coord + source->x_offset * DEFdist2Microns;
            box.yLL = box.yUR = y_coord + source->y_offset * DEFdist2Microns;
        } else {
            box.xLL = box.xUR = source->x_coord;
            box.yLL = box.yUR = source->y_coord;
        }

        for (int j = 0; j < theNet->sinks.size(); j++) {
            pin *sink = &pins[theNet->sinks[j]];
            //cout << " sink name : " << sink->name << endl;
            if (sink->type == NONPIO_PIN) {
                cell *theCell = &cells[sink->owner];
                if (mode == "INIT") {
                    x_coord = theCell->init_x_coord;
                    y_coord = theCell->init_y_coord;
                } else {
                    x_coord = theCell->x_coord;
                    y_coord = theCell->y_coord;
                }
                box.xLL = min(box.xLL, x_coord + sink->x_offset * DEFdist2Microns);
                box.xUR = max(box.xUR, x_coord + sink->x_offset * DEFdist2Microns);
                box.yLL = min(box.yLL, y_coord + sink->y_offset * DEFdist2Microns);
                box.yUR = max(box.yUR, y_coord + sink->y_offset * DEFdist2Microns);
            } else {
                box.xLL = min(box.xLL, sink->x_coord);
                box.xUR = max(box.xUR, sink->x_coord);
                box.yLL = min(box.yLL, sink->y_coord);
                box.yUR = max(box.yUR, sink->y_coord);
            }
        }

        double box_boundary = (box.xUR - box.xLL + box.yUR - box.yLL);

        hpwl += box_boundary;
    }
    return hpwl / static_cast<double>(DEFdist2Microns);
}

double circuit::calc_density_factor(double unit) {

    double gridUnit = unit * rowHeight;
    int x_gridNum = (int) ceil((rx - lx) / gridUnit); //rx--> right / lx --> left
    int y_gridNum = (int) ceil((ty - by) / gridUnit); //ty --> top / by --> bottom
    int numBins = x_gridNum * y_gridNum;

    // Initialize density map
    vector<density_bin> bins(numBins);
    for (int j = 0; j < y_gridNum; j++) {
        for (int k = 0; k < x_gridNum; k++) {
            unsigned binId = j * x_gridNum + k;
            bins[binId].lx = lx + k * gridUnit;
            bins[binId].ly = by + j * gridUnit;
            bins[binId].hx = bins[binId].lx + gridUnit;
            bins[binId].hy = bins[binId].ly + gridUnit;

            bins[binId].hx = min(bins[binId].hx, rx);
            bins[binId].hy = min(bins[binId].hy, ty);

            bins[binId].area = max((bins[binId].hx - bins[binId].lx) * (bins[binId].hy - bins[binId].ly), 0.0);
            bins[binId].m_util = 0.0;
            bins[binId].f_util = 0.0;
            bins[binId].free_space = 0.0;
            bins[binId].overflow = 0.0;
        }
    }

    /* 1. build density map */
    /* (a) calculate overlaps with row sites, and add them to free_space */
    for (vector<row>::iterator theRow = rows.begin(); theRow != rows.end(); ++theRow) {
        int lcol = max((int) floor((theRow->origX - lx) / gridUnit), 0);
        int rcol = min((int) floor((theRow->origX + theRow->numSites * theRow->stepX - lx) / gridUnit), x_gridNum - 1);
        int brow = max((int) floor((theRow->origY - by) / gridUnit), 0);
        int trow = min((int) floor((theRow->origY + rowHeight - by) / gridUnit), y_gridNum - 1);

        for (int j = brow; j <= trow; j++) {
            for (int k = lcol; k <= rcol; k++) {
                unsigned binId = j * x_gridNum + k;

                /* get intersection */
                double lx = max(bins[binId].lx, (double) theRow->origX);
                double hx = min(bins[binId].hx, (double) theRow->origX + theRow->numSites * theRow->stepX);
                double ly = max(bins[binId].ly, (double) theRow->origY);
                double hy = min(bins[binId].hy, (double) theRow->origY + rowHeight);

                if ((hx - lx) > 1.0e-5 && (hy - ly) > 1.0e-5) {
                    double common_area = (hx - lx) * (hy - ly);
                    bins[binId].free_space += common_area;
                    bins[binId].free_space = min(bins[binId].free_space, bins[binId].area);
                }
            }
        }
    }

    /* (b) add utilization by fixed/movable objects */
    for (vector<cell>::iterator theCell = cells.begin(); theCell != cells.end(); ++theCell) {
        int lcol = max((int) floor((theCell->init_x_coord - lx) / gridUnit), 0);
        int rcol = min((int) floor((theCell->init_x_coord + theCell->width - lx) / gridUnit), x_gridNum - 1);
        int brow = max((int) floor((theCell->init_y_coord - by) / gridUnit), 0);
        int trow = min((int) floor((theCell->init_y_coord + theCell->height - by) / gridUnit), y_gridNum - 1);

        for (int j = brow; j <= trow; j++) {
            for (int k = lcol; k <= rcol; k++) {
                unsigned binId = j * x_gridNum + k;

                if (theCell->inGroup)
                    bins[binId].density_limit = max(bins[binId].density_limit, groups[group2id[theCell->group]].util);

                /* get intersection */
                double lx = max(bins[binId].lx, (double) theCell->init_x_coord);
                double hx = min(bins[binId].hx, (double) theCell->init_x_coord + theCell->width);
                double ly = max(bins[binId].ly, (double) theCell->init_y_coord);
                double hy = min(bins[binId].hy, (double) theCell->init_y_coord + theCell->height);
                double x_center = (double) theCell->init_x_coord + (double) theCell->width / 2;
                double y_center = (double) theCell->init_y_coord + (double) theCell->height / 2;

                if (bins[binId].lx <= x_center && x_center < bins[binId].hx)
                    if (bins[binId].ly < y_center && y_center < bins[binId].hy)
                        theCell->binId = binId;

                if ((hx - lx) > 1.0e-5 && (hy - ly) > 1.0e-5) {
                    double common_area = (hx - lx) * (hy - ly);
                    if (theCell->isFixed)
                        bins[binId].f_util += common_area;
                    else
                        bins[binId].m_util += common_area;
                }
            }
        }
    }
    //////////////////////////////// SLOW CALC in CTYPES //////
    int num_cells = 0;
    double totalArea = 0.0;
    double singleBinArea = 64 * rowHeight * rowHeight;
    double totalOverflow = 0.0;
    double totalMovableArea = 0.0;
    for (int i = 0; i < y_gridNum * x_gridNum; i++) {
        //double freeSpace = 0.0;
        double movableArea = 0.0;
        for (vector<cell>::iterator theCell = cells.begin(); theCell != cells.end(); ++theCell) {
            //double movableArea = 0.0;
            //if(theCell->binId == i){
            double lx = max(bins[i].lx, (double) theCell->init_x_coord);
            double hx = min(bins[i].hx, (double) theCell->init_x_coord + theCell->width);
            double ly = max(bins[i].ly, (double) theCell->init_y_coord);
            double hy = min(bins[i].hy, (double) theCell->init_y_coord + theCell->height);
            double common_area = (hx - lx) * (hy - ly);
            if ((hx - lx) > 1.0e-5 && (hy - ly) > 1.0e-5) {
                // totalArea += bins[i].area;
                if (theCell->isFixed) {
                    bins[i].free_space -= common_area;
                } else
                    movableArea += common_area;

                // totalMovableArea += movableArea;
                // double overflow = movableArea - freeSpace * max_utilization/100;
                // totalOverflow += std::max(overflow, 0.0);
            }
            //double overflow = movableArea - freeSpace * max_utilization/100;
            //}
            num_cells++;
        }
        totalMovableArea += movableArea;
        double overflow = movableArea - bins[i].free_space * max_utilization / 100;
        totalOverflow += overflow;
        //cout<<"free space     : " <<bins[i].free_space<<endl;
        //cout<<"movable cells  : " <<movableArea << endl;
    }
    //cout<<num_cells/(y_gridNum*x_gridNum)<<endl;
    //cout << totalOverflow << endl;
    //cout << totalMovableArea << endl;


    for (vector<cell>::iterator theCell = cells.begin(); theCell != cells.end(); ++theCell) {
        if (theCell->binId == UINT_MAX)
            continue;
        density_bin *theBin = &bins[theCell->binId];
        theCell->dense_factor += theBin->m_util / (theBin->free_space - theBin->f_util);
        theCell->binId = UINT_MAX;
    }

    return totalOverflow * singleBinArea * (max_utilization / 100) / totalMovableArea;
}

void circuit::group_analyze() {
    for (int i = 0; i < groups.size(); i++) {
        group *theGroup = &groups[i];
        double region_area = 0;
        double avail_region_area = 0;
        double cell_area = 0;
        for (int j = 0; j < theGroup->regions.size(); j++) {
            rect *theRect = &theGroup->regions[j];
            region_area += (theRect->xUR - theRect->xLL) * (theRect->yUR - theRect->yLL);
            avail_region_area +=
                    (theRect->xUR - theRect->xLL - (int) theRect->xUR % 200 + (int) theRect->xLL % 200 - 200)
                    * (theRect->yUR - theRect->yLL - (int) theRect->yUR % 2000 + (int) theRect->yLL % 2000 - 2000);
        }
        for (int k = 0; k < theGroup->siblings.size(); k++) {
            cell *theCell = theGroup->siblings[k];
            cell_area += theCell->width * theCell->height;
        }
        cout << " GROUP : " << theGroup->name << endl;
        cout << " region count : " << theGroup->regions.size() << endl;
        cout << " cell count : " << theGroup->siblings.size() << endl;
        cout << " region area : " << region_area << endl;
        cout << " avail region area : " << avail_region_area << endl;
        cout << " cell area : " << cell_area << endl;
        cout << " utilization : " << cell_area / region_area << endl;
        cout << " avail util : " << cell_area / avail_region_area << endl;
        cout << " - - - - - - - - - - - - - - - - - - - - " << endl;
    }
    return;
}

pair<int, int> circuit::nearest_coord_to_rect_boundary(cell *theCell, rect *theRect, string mode) {
    int x = INT_MAX;
    int y = INT_MAX;
    int size_x = (int) floor(theCell->width / wsite + 0.5);
    int size_y = (int) floor(theCell->height / rowHeight + 0.5);
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::nearest_coord_to_rect_boundary == invalid mode!" << endl;
        exit(2);
    }
    int temp_x = x;
    int temp_y = y;

    if (check_overlap(theCell, theRect, "init_coord") == true) {
        int dist_x = 0;
        int dist_y = 0;
        if (abs(x - theRect->xLL + theCell->width) > abs(theRect->xUR - x)) {
            dist_x = abs(theRect->xUR - x);
            temp_x = theRect->xUR;
        } else {
            dist_x = abs(x - theRect->xLL);
            temp_x = theRect->xLL - theCell->width;
        }
        if (abs(y - theRect->yLL + theCell->height) > abs(theRect->yUR - y)) {
            dist_y = abs(theRect->yUR - y);
            temp_y = theRect->yUR;
        } else {
            dist_y = abs(y - theRect->yLL);
            temp_y = theRect->yLL - theCell->height;
        }
        assert (dist_x > -1);
        assert (dist_y > -1);
        if (dist_x < dist_y)
            return make_pair(temp_x, y);
        else
            return make_pair(x, temp_y);
    }

    if (x < theRect->xLL)
        temp_x = theRect->xLL;
    else if (x + theCell->width > theRect->xUR)
        temp_x = theRect->xUR - theCell->width;

    if (y < theRect->yLL)
        temp_y = theRect->yLL;
    else if (y + theCell->height > theRect->yUR)
        temp_y = theRect->yUR - theCell->height;

#ifdef DEBUG
    cout << " - - - - - - - - - - - - - - - " << endl;
    cout << " input x_coord : " << x << endl;
    cout << " input y_coord : " << y << endl;
    cout << " found x_coord : " << temp_x << endl;
    cout << " found y_coord : " << temp_y << endl;
#endif

    return make_pair(temp_x, temp_y);
}

int circuit::dist_for_rect(cell *theCell, rect *theRect, string mode) {
    int x = INT_MAX;
    int y = INT_MAX;
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::dist_for_rect == invalid mode!" << endl;
        exit(2);
    }
    int temp_x = 0;
    int temp_y = 0;

    if (x < theRect->xLL)
        temp_x = theRect->xLL - x;
    else if (x + theCell->width > theRect->xUR)
        temp_x = x + theCell->width - theRect->xUR;

    if (y < theRect->yLL)
        temp_y = theRect->yLL - y;
    else if (y + theCell->height > theRect->yUR)
        temp_y = y + theCell->height - theRect->yUR;

    assert(temp_y > -1);
    assert(temp_x > -1);

    return temp_y + temp_x;
}

bool circuit::check_overlap(rect cell, rect box) {
    if (box.xLL >= cell.xUR || box.xUR <= cell.xLL)
        return false;
    if (box.yLL >= cell.yUR || box.yUR <= cell.yLL)
        return false;
    return true;
}

bool circuit::check_overlap(cell *theCell, rect *theRect, string mode) {

    int x = INT_MAX;
    int y = INT_MAX;
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::check_overlap == invalid mode!" << endl;
        exit(2);
    }

    if (theRect->xUR <= x || theRect->xLL >= x + theCell->width)
        return false;
    if (theRect->yUR <= y || theRect->yLL >= y + theCell->height)
        return false;

    return true;
}

bool circuit::check_inside(rect cell, rect box) {
    if (box.xLL > cell.xLL || box.xUR < cell.xUR)
        return false;
    if (box.yLL > cell.yLL || box.yUR < cell.yUR)
        return false;
    return true;
}

bool circuit::check_inside(cell *theCell, rect *theRect, string mode) {

    int x = INT_MAX;
    int y = INT_MAX;
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::check_inside == invalid mode!" << endl;
        exit(2);
    }

    if (theRect->xUR < x + theCell->width || theRect->xLL > x)
        return false;
    if (theRect->yUR < y + theCell->height || theRect->yLL > y)
        return false;

    return true;
}

pair<bool, pair<int, int> > circuit::bin_search(int x_pos, cell *theCell, int x, int y) {
    pair<int, int> pos;
    macro *theMacro = &macros[theCell->type];

    // EDGETYPE 1 - 1 : 400, 1 - 2 : 400, 2 - 2 : 0
    int edge_left = 0;
    int edge_right = 0;
    if (theMacro->edgetypeLeft == 1)
        edge_left = 2;
    if (theMacro->edgetypeRight == 1)
        edge_right = 2;

    int x_step = (int) ceil(theCell->width / wsite) + edge_left + edge_right;
    int y_step = (int) ceil(theCell->height / rowHeight);

    // IF y is out of border
    if (y + y_step > (die.yUR / rowHeight))
        return make_pair(false, pos);

    // If even number multi-deck cell -> check top power
    if (y_step % 2 == 0) {
        if (rows[y].top_power == theMacro->top_power)
            return make_pair(false, pos);
    }

#ifdef DEBUG
    cout << " - - - - - - - - - - - - - - - - - " << endl;
    cout << " Start Bin Search " << endl;
    cout << " cell name : " << theCell->name << endl;
    cout << " target x : " << x << endl;
    cout << " target y : " << y << endl;
#endif

    if (x_pos > x) {
        for (int i = 9; i > -1; i--) {

            // check all grids are empty
            bool available = true;

            if (x + i + x_step > (int) (die.xUR / wsite)) {
                available = false;
            } else {
                for (int k = y; k < y + y_step; k++) {
                    for (int l = x + i; l < x + i + x_step; l++) {
                        if (grid[k][l].linked_cell != NULL || grid[k][l].isValid == false) {
                            available = false;
                            break;
                        }
                        // check group regions
                        if (theCell->inGroup == true) {
                            if (grid[k][l].group != group2id[theCell->group])
                                available = false;
                        } else {
                            if (grid[k][l].group != UINT_MAX)
                                available = false;
                        }
                    }
                    if (available == false)
                        break;
                }
            }
            if (available == true) {
                if (edge_left == 0)
                    pos = make_pair(y, x + i);
                else
                    pos = make_pair(y, x + i + edge_left);

                return make_pair(available, pos);
            }
        }
    } else if (x_pos < x) {
        for (int i = 0; i < 10; i++) {
            bool available = true;
            if (x + i + x_step > (int) (die.xUR / wsite)) {
                available = false;
            } else {
                for (int k = y; k < y + y_step; k++) {
                    for (int l = x + i; l < x + i + x_step; l++) {
                        if (grid[k][l].linked_cell != NULL || !grid[k][l].isValid) {
                            available = false;
                            break;
                        }
                        // check group regions
                        if (theCell->inGroup) {
                            if (grid[k][l].group != group2id[theCell->group])
                                available = false;
                        } else {
                            if (grid[k][l].group != UINT_MAX)
                                available = false;
                        }
                    }
                    if (!available) break;
                }
            }
            if (available) {
                pos = make_pair(y, x + i + edge_left);
                return make_pair(available, pos);
            }
        }
    } else {
        //////////////////////////// ADDED ALGORITHM WITH SIGNFACTOR //////
        for (int i = 0; i < 10; i++) {
            int signFactor;
            for (int j = 0; j < 2; ++j) {
                if (j == 0) signFactor = 1;
                else signFactor = -1;

                // check all grids are empty
                bool available = true;
                if (x + i * signFactor + x_step > (int) (die.xUR / wsite)) {
                    available = false;
                } else {
                    for (int k = y; k < y + y_step; k++) {
                        for (int l = x + i * signFactor; l < x + i * signFactor + x_step; l++) {
                            if (grid[k][l].linked_cell != NULL || grid[k][l].isValid == false) {
                                available = false;
                                break;
                            }
                            // check group regions
                            if (theCell->inGroup == true) {
                                if (grid[k][l].group != group2id[theCell->group])
                                    available = false;
                            } else {
                                if (grid[k][l].group != UINT_MAX)
                                    available = false;
                            }
                        }
                        if (available == false)
                            break;
                    }
                }
                if (available == true) {
#ifdef DEBUG
                    cout << " found pos x - y : " << x << " - " << y << " Finish Search " << endl;
                    cout << " - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
#endif
                    if (edge_left == 0)
                        pos = make_pair(y, x + i * signFactor);
                    else
                        pos = make_pair(y, x + i * signFactor + edge_left);

                    return make_pair(available, pos);
                }
            }
        }
    }
    return make_pair(false, pos);
}

pair<bool, pixel *> circuit::diamond_search(cell *theCell, int x_coord, int y_coord) {
    pixel *myPixel = NULL;
    pair<bool, pair<int, int> > found;
    int x_pos = (int) floor(x_coord / wsite + 0.5);
    int y_pos = (int) floor(y_coord / rowHeight + 0.5);

    int x_start = 0;
    int x_end = 0;
    int y_start = 0;
    int y_end = 0;


    // set search boundary max / min
    if (theCell->inGroup == true) {
        group *theGroup = &groups[group2id[theCell->group]];
        x_start = max(x_pos - (int) (displacement * 5), (int) floor(theGroup->boundary.xLL / wsite));
        x_end = min(x_pos + (int) (displacement * 5),
                    (int) floor(theGroup->boundary.xUR / wsite) - (int) floor(theCell->width / rowHeight + 0.5));
        y_start = max(y_pos - (int) displacement, (int) ceil(theGroup->boundary.yLL / rowHeight));
        y_end = min(y_pos + (int) displacement,
                    (int) ceil(theGroup->boundary.yUR / rowHeight) - (int) floor(theCell->height / rowHeight + 0.5));
    } else {
        x_start = max(x_pos - (int) (displacement * 5), 0);
        x_end = min(x_pos + (int) (displacement * 5),
                    (int) floor(rx / wsite) - (int) floor(theCell->width / wsite + 0.5));
        y_start = max(y_pos - (int) displacement, 0);
        y_end = min(y_pos + (int) displacement,
                    (int) floor(ty / rowHeight) - (int) floor(theCell->height / rowHeight + 0.5));
    }
#ifdef DEBUG
    cout << " == Start Diamond Search ==  " << endl;
    cout << " cell_name : " << theCell->name << endl;
    cout << " cell width : " << theCell->width << endl;
    cout << " cell height : " << theCell->height << endl;
    cout << " cell x step : " << (int)floor(theCell->width/wsite+0.5) << endl;
    cout << " cell y step : " << (int)floor(theCell->height/rowHeight+0.5) << endl;
    cout << " input x : " << x_coord << endl;
    cout << " inpuy y : " << y_coord << endl;
    cout << " x_pos : " << x_pos << endl;
    cout << " y_pos : " << y_pos << endl;
    cout << " x bound ( " << x_start << ") - (" << x_end << ")" << endl;
    cout << " y bound ( " << y_start << ") - (" << y_end << ")" << endl;
#endif
    found = bin_search(x_pos, theCell, min(x_end, max(x_start, x_pos)), max(y_start, min(y_end, y_pos)));
    if (found.first == true) {
        myPixel = &grid[found.second.first][found.second.second];
        return make_pair(found.first, myPixel);
    }

    int div = 4;
    if (design_util > 0.6 || num_fixed_nodes > 0)
        div = 1;

    for (int i = 1; i < (int) (displacement * 2) / div; i++) {
        vector<pixel *> avail_list;
        avail_list.reserve(i * 4);
        pixel *myPixel = NULL;

        int x_offset = 0;
        int y_offset = 0;

        // right side
        for (int j = 1; j < i * 2; j++) {
            x_offset = -((j + 1) / 2);
            if (j % 2 == 1)
                y_offset = (i * 2 - j) / 2;
            else
                y_offset = -(i * 2 - j) / 2;
            found = bin_search(x_pos, theCell, min(x_end, max(x_start, (x_pos + x_offset * 10))),
                               min(y_end, max(y_start, (y_pos + y_offset))));
            if (found.first == true) {
                myPixel = &grid[found.second.first][found.second.second];
                avail_list.push_back(myPixel);
            }
        }

        // left side
        for (int j = 1; j < (i + 1) * 2; j++) {
            x_offset = (j - 1) / 2;
            if (j % 2 == 1)
                y_offset = ((i + 1) * 2 - j) / 2;
            else
                y_offset = -((i + 1) * 2 - j) / 2;
            found = bin_search(x_pos, theCell, min(x_end, max(x_start, (x_pos + x_offset * 10))),
                               min(y_end, max(y_start, (y_pos + y_offset))));
            if (found.first == true) {
                myPixel = &grid[found.second.first][found.second.second];
                avail_list.push_back(myPixel);
            }
        }

        // check from vector
        unsigned dist = UINT_MAX;
        int best = INT_MAX;
        for (int j = 0; j < avail_list.size(); j++) {
            int temp_dist =
                    abs(x_coord - avail_list[j]->x_pos * wsite) + abs(y_coord - avail_list[j]->y_pos * rowHeight);
            if (temp_dist < dist) {
                dist = temp_dist;
                best = j;
            }
        }
        if (best != INT_MAX) {
            return make_pair(true, avail_list[best]);
        }
    }
    return make_pair(false, myPixel);
}

bool circuit::direct_move(cell *theCell, string mode) {
    int x = INT_MAX;
    int y = INT_MAX;
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::check_inside == invalid mode!" << endl;
        exit(2);
    }
    return direct_move(theCell, x, y);
}

bool circuit::direct_move(cell *theCell, int x_coord, int y_coord) {
    pixel *myPixel = NULL;
    pair<bool, pair<int, int> > found;
    int x_pos = (int) floor(x_coord / wsite + 0.5);
    int y_pos = (int) floor(y_coord / rowHeight + 0.5);

    int x_step = (int) ceil(theCell->width / wsite);
    int y_step = (int) ceil(theCell->height / rowHeight);
    int x_end = x_pos + x_step;
    int y_end = y_pos + y_step;

    if (x_pos < (int) ceil(die.xLL / wsite) || x_end > (int) floor(die.xUR / wsite))
        return false;
    else if (y_pos < (int) ceil(die.yLL / rowHeight) || y_end > (int) floor(die.yUR / rowHeight))
        return false;
    else {
        for (int i = y_pos; i < y_end; i++) {
            for (int j = x_pos; j < x_end; j++) {
                if (grid[i][j].linked_cell != NULL || grid[i][j].isValid == true)
                    return false;
            }
        }
        paint_pixel(theCell, x_pos, y_pos);
        return true;
    }
}


bool circuit::shift_move(cell *theCell, int x, int y) {
//	cout << " shift_move start " << endl;

    // set resgion boundary
    rect theRect;
    theRect.xLL = max(die.xLL, x - theCell->width * 3);
    theRect.xUR = min(die.xUR, x + theCell->width * 3);
    theRect.yLL = max(die.yLL, y - theCell->height * 3);
    theRect.yUR = min(die.yUR, y + theCell->height * 3);

    //vector<cell*> overlap_region_cells = get_cells_from_boundary(&theRect);
    //overlap_region_cells.clear();
    //vector<cell*>().swap(overlap_region_cells);
    overlap_region_cells = get_cells_from_boundary(&theRect);

    // erase region cells
    for (int i = 0; i < overlap_region_cells.size(); i++) {
        cell *around_cell = overlap_region_cells[i];
        if (theCell->inGroup == around_cell->inGroup) {
            //assert ( check_inside(around_cell,&theRect,"coord") == true );
            erase_pixel(around_cell);
        }
    }

    // place target cell
    if (map_move(theCell, x, y) == false) {
        cout << " can't insert center cell !! " << endl;
        cout << " cell_name : " << theCell->name << endl;
        return false;
    }

    // rebuild erased cells
    for (int i = 0; i < overlap_region_cells.size(); i++) {
        cell *around_cell = overlap_region_cells[i];
        if (theCell->inGroup == around_cell->inGroup) {
            if (map_move(around_cell, around_cell->init_x_coord, around_cell->init_y_coord) == false) {
#ifdef DEBUG
                cout << " Shift move fail !!" << endl;
                cout << " cell name : " << around_cell->name << endl;
                cout << " x_coord : " << around_cell->init_x_coord << endl;
                cout << " y_coord : " << around_cell->init_y_coord << endl;
#endif
                return false;
            }
        }
    }
    return true;
}

bool circuit::shift_move(cell *theCell, string mode) {
    int x = INT_MAX;
    int y = INT_MAX;
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::shift_move() == invalid mode!" << endl;
        exit(2);
    }
    return shift_move(theCell, x, y);
}

bool circuit::map_move(cell *theCell, string mode) {
    int x = INT_MAX;
    int y = INT_MAX;
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::check_inside == invalid mode!" << endl;
        exit(2);
    }
    return map_move(theCell, x, y);
}

bool circuit::map_move(cell *theCell, int x, int y) {
    pair<bool, pixel *> myPixel = diamond_search(theCell, x, y);
    if (myPixel.first == true) {
        //pair<bool,pixel*> nearPixel = diamond_search(theCell,myPixel.second->x_pos*wsite,myPixel.second->y_pos*rowHeight);
        //if( nearPixel.first == true ) {
        //    paint_pixel(theCell,nearPixel.second->x_pos,nearPixel.second->y_pos);
        //    //cout << " near found!! " << endl;
        //}
        //else
        paint_pixel(theCell, myPixel.second->x_pos, myPixel.second->y_pos);

        /********/
        double X = abs(theCell->init_x_coord - theCell->x_coord);
        double Y = abs(theCell->init_y_coord - theCell->y_coord);
        double DISP = X + Y;

        // reward option (1): -disp
        //REWARD = -1.0 * DISP;

        // reward option (2): k/disp
        if (DISP < 200) DISP = 200.0;
        REWARD += 200.0 / (DISP);
        /********/
        return true;

    } else {
#ifdef DEBUG
        cout << " Map move fail !!" << endl;
        cout << " cell name : " << theCell->name << endl;
        cout << " init_x_coord : " << theCell->init_x_coord << endl;
        cout << " init_y_coord : " << theCell->init_y_coord << endl;
#endif
        return false;
    }
}

vector<cell *> circuit::overlap_cells(cell *theCell) {
    vector<cell *> list;
    int step_x = (int) ceil(theCell->width / wsite);
    int step_y = (int) ceil(theCell->height / rowHeight);

    dense_hash_map<unsigned, cell *> cell_list;
    cell_list.set_empty_key(UINT_MAX);

    for (int i = theCell->y_pos; i < theCell->y_pos + step_y; i++) {
        for (int j = theCell->x_pos; j < theCell->y_pos + step_x; j++) {
            if (grid[i][j].linked_cell != NULL) {
                cell_list[grid[i][j].linked_cell->id] = grid[i][j].linked_cell;
            }
        }
    }
    list.reserve(cell_list.size());
    for (auto &currCell: cell_list) {
        list.push_back(currCell.second);
    }

    return list;
}

// rect should be position
vector<cell *> circuit::get_cells_from_boundary(rect *theRect) {

    assert(theRect->xLL >= die.xLL);
    assert(theRect->yLL >= die.yLL);
    assert(theRect->xUR <= die.xUR);
    assert(theRect->yUR <= die.yUR);

    int x_start = (int) floor(theRect->xLL / wsite + 0.5);
    int y_start = (int) floor(theRect->yLL / rowHeight + 0.5);
    int x_end = (int) floor(theRect->xUR / wsite + 0.5);
    int y_end = (int) floor(theRect->yUR / rowHeight + 0.5);

    vector<cell *> list;

    dense_hash_map<unsigned, cell *> cell_list;
    cell_list.set_empty_key(UINT_MAX);

    for (int i = y_start; i < y_end; i++) {
        for (int j = x_start; j < x_end; j++) {
            if (grid[i][j].linked_cell != NULL) {
                if (grid[i][j].linked_cell->isFixed == false)
                    cell_list[grid[i][j].linked_cell->id] = grid[i][j].linked_cell;
            }
        }
    }

    list.reserve(cell_list.size());
    for (auto &currCell: cell_list) {
        list.push_back(currCell.second);
    }
    return list;
}

double circuit::dist_benefit(cell *theCell, int x_coord, int y_coord) {
    double curr_dist = abs(theCell->x_coord - theCell->init_x_coord) + abs(theCell->y_coord - theCell->init_y_coord);
    double new_dist = abs(theCell->init_x_coord - x_coord) + abs(theCell->init_y_coord - y_coord);
    return new_dist - curr_dist;
}

bool circuit::swap_cell(cell *cellA, cell *cellB) {
    if (cellA == cellB)
        return false;
    else if (cellA->type != cellB->type)
        return false;
    else if (cellA->isFixed == true || cellB->isFixed == true)
        return false;

    double benefit =
            dist_benefit(cellA, cellB->x_coord, cellB->y_coord) + dist_benefit(cellB, cellA->x_coord, cellA->y_coord);

    if (benefit < 0) {
        int A_x_pos = cellB->x_pos;
        int A_y_pos = cellB->y_pos;
        int B_x_pos = cellA->x_pos;
        int B_y_pos = cellA->y_pos;

        erase_pixel(cellA);
        erase_pixel(cellB);
        paint_pixel(cellA, A_x_pos, A_y_pos);
        paint_pixel(cellB, B_x_pos, B_y_pos);
        //cout << "swap benefit : " << benefit << endl;
        //save_score();
        return true;
    }
    return false;
}

bool circuit::refine_move(cell *theCell, string mode) {
    int x = INT_MAX;
    int y = INT_MAX;
    if (mode == "init_coord") {
        x = theCell->init_x_coord;
        y = theCell->init_y_coord;
    } else if (mode == "coord") {
        x = theCell->x_coord;
        y = theCell->y_coord;
    } else if (mode == "pos") {
        x = theCell->x_pos * wsite;
        y = theCell->y_pos * rowHeight;
    } else {
        cerr << "circuit::refine_move == invalid mode!" << endl;
        exit(2);
    }
    return refine_move(theCell, x, y);
}

//
bool circuit::refine_move(cell *theCell, int x_coord, int y_coord) {
    pair<bool, pixel *> myPixel = diamond_search(theCell, x_coord, y_coord);
    if (myPixel.first == true) {
        double new_dist = abs(theCell->init_x_coord - myPixel.second->x_pos * wsite) +
                          abs(theCell->init_y_coord - myPixel.second->y_pos * rowHeight);
        if (new_dist / rowHeight > max_disp_const)
            return false;

        double benefit = dist_benefit(theCell, myPixel.second->x_pos * wsite, myPixel.second->y_pos * rowHeight);
        //if( benefit < 2001-sum_displacement/20 ) {
        if (benefit < 0) {
            //cout << " refine benefit : " << benefit << " : " << 2001 - sum_displacement/10 << endl;
            sum_displacement++;
            erase_pixel(theCell);
            paint_pixel(theCell, myPixel.second->x_pos, myPixel.second->y_pos);
            //save_score();
            return true;
        } else
            return false;
    } else
        return false;
}

pair<bool, cell *> circuit::nearest_cell(int x_coord, int y_coord) {
    bool found = false;
    cell *nearest_cell = NULL;
    double nearest_dist = 99999999999;
    for (int i = 0; i < cells.size(); i++) {
        cell *theCell = &cells[i];
        if (theCell->isPlaced == false)
            continue;

        double dist = abs(theCell->x_coord - x_coord) + abs(theCell->y_coord - y_coord);

        if (dist < rowHeight * 2)
            if (nearest_dist > dist) {
                nearest_dist = dist;
                nearest_cell = theCell;
                found = true;
            }
    }
    return make_pair(found, nearest_cell);
}

void circuit::avg_min2_search_distance() {
    rect designRect;
    designRect.xLL = lx;
    designRect.yLL = by;
    designRect.xUR = rx;
    designRect.yUR = ty;

    for (int i = 0; i < (int) cells.size(); i++) {
        cell *theCell = &cells[i];
        if (theCell->inGroup || theCell->isPlaced || theCell->isFixed) continue;

        pair<double, double> coord;
        pair<double, double> left;
        pair<double, double> right;
        pair<double, double> lower;
        pair<double, double> upper;

        vector<double> dists;

        for (int j = 0; j < (int) groups.size(); j++) {
            group *theGroup = &groups[j];
            for (int k = 0; k < (int) theGroup->regions.size(); k++) {
                rect *theRect = &theGroup->regions[k];
                coord = nearest_coord_to_rect_boundary(theCell, theRect, "coord");
                double dist =
                        abs((double) theCell->x_coord - coord.first) + abs((double) theCell->y_coord - coord.second);
                dists.push_back(dist);
            }
        }
        for (int j = 0; j < num_fixed_nodes; j++) {
            rect *theRect = &fixedCells[j];
            coord = nearest_coord_to_rect_boundary(theCell, theRect, "coord");
            if (theRect->xLL < theCell->x_coord && theRect->xUR > theCell->x_coord + theCell->width &&
                theRect->yLL < theCell->y_coord && theRect->yUR > theCell->y_coord + theCell->width) {
                double dist = 0.0;
                dists.push_back(dist);
            } else {
                double dist =
                        abs((double) theCell->x_coord - coord.first) + abs((double) theCell->y_coord - coord.second);
                dists.push_back(dist);
            }
        }

        //coord = nearest_coord_to_rect_boundary(theCell, &designRect, "coord");
        //double dist  = abs((double)theCell->x_coord - coord.first) + abs((double)theCell->y_coord - coord.second);
        double dist = 0.0;
        dist = ty - (theCell->y_coord + theCell->height);
        dists.push_back(dist);
        dist = (theCell->y_coord) - by;
        dists.push_back(dist);
        dist = rx - (theCell->x_coord + theCell->width);
        dists.push_back(dist);
        dist = theCell->x_coord - lx;
        dists.push_back(dist);

        sort(dists.begin(), dists.end());

        double avg2Dist = (dists[0] + dists[1]) / 2.0;
        theCell->min2SrchDist = avg2Dist;
    }
}


void clear_agent(circuit *ck, Agent *agent) {
    agent->state.cell_list.clear();
    vector<cell *>().swap(agent->state.cell_list);
    for (auto it = agent->state.features.begin(); it != agent->state.features.end(); it++) {
        (*it).clear();
        vector<double>().swap(*it);
    }
    agent->state.features.clear();
    vector<vector<double>>().swap(agent->state.features);
    agent->effCells.clear();
    vector<int>().swap(agent->effCells);
    agent->tarCell = nullptr;
    agent->moveFailCnt = 0;
    agent->Pov = 0;

    cout << "Gcell Reset!" << endl;
}

void clear_memory(circuit *ck, Agent *agent) {
    agent->state.cell_list.clear();
    vector<cell *>().swap(agent->state.cell_list);
    for (auto it = agent->state.features.begin(); it != agent->state.features.end(); it++) {
        (*it).clear();
        vector<double>().swap(*it);
    }
    agent->state.features.clear();
    vector<vector<double>>().swap(agent->state.features);
    agent->state.features.clear();
    agent->effCells.clear();
    vector<int>().swap(agent->effCells);
    agent->tarCell = nullptr;
    agent->moveFailCnt = 0;
    agent->Pov = 0;

    for (int i = 0; i < (int) ck->cells.size(); i++) {
        cell *theCell = &(ck->cells[i]);
        //if(theCell->isFixed) continue;
        if (theCell->isFixed || theCell->inGroup) continue;      /// add inGroup continue for saving runtime
        ck->erase_pixel(theCell);
        theCell->isPlaced = false;
        theCell->moveTry = false;
        theCell->fieldId = UINT_MAX;
        theCell->gcellId = UINT_MAX;
        theCell->sindex = -1;

        theCell->overlapNum = 0;

        // initalize x,y-coord of cells to init-coord
        if (theCell->init_x_coord < ck->lx) {
            theCell->x_coord = ck->lx;
        } else if (theCell->x_coord > ck->rx) {
            cout << "[WARNING] The cell-" << theCell->id << " is out of x-range of design!" << endl;
            theCell->x_coord = ck->rx;
        } else {
            theCell->x_coord = theCell->init_x_coord;
        }
        if (theCell->init_y_coord < ck->by) {
            theCell->y_coord = ck->by;
        } else if (theCell->init_y_coord > ck->ty) {
            cout << "[WARNING] The cell-" << theCell->id << " is out of y-range of design!" << endl;
            theCell->y_coord = ck->ty;
        } else {
            theCell->y_coord = theCell->init_y_coord;
        }
    }
    ck->Rtree_clear();
    ck->Gcells.clear();
    ck->Fields.clear();


    //cout << "Memory cleared!" << endl;
    cout << "Episode Reset!" << endl;
}

void circuit::gcell_unit_calc(int avgCellNum) {

    int cell_num = non_group_cell_cnt();
    double total_FR = 0.0;
    for (int i = 0; i < groups.size(); i++) {
        group *theGroup = &groups[i];
        for (int j = 0; j < theGroup->regions.size(); j++) {
            rect *R = &theGroup->regions[j];
            total_FR += (R->xUR - R->xLL) * (R->yUR - R->yLL);
        }
    }

    double free_space = designArea - total_fArea - total_FR;
    double cell_fs = (double) (free_space / cell_num);
    double length = sqrt(avgCellNum * cell_fs);

    length /= rowHeight;
    length = round(length / field_unit) * (int) field_unit;

    gcell_unit = length;
}


int circuit::non_group_cell_cnt() {

    int cell_num = 0;
    double cell_area = 0;
    for (auto theCell: cells) {
        if (theCell.inGroup || theCell.isPlaced || theCell.isFixed) continue;
        cell_area += theCell.width * theCell.height;
        cell_num++;
    }
    return cell_num;
}


