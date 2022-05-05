#include "circuit.h"


int circuit::overlap_num_calc(cell *theCell) {
    vector<cell *> ovcells = overlap_cells(theCell);
    ovcells.erase(unique(ovcells.begin(), ovcells.end()), ovcells.end());
    int num = (int) ovcells.size();
    /////////// num = 0 for all cells...

    if (num != 0)
        cout << "CELL " << theCell->name << " has " << num << "overlap cells" << endl;

    return num;

}

/*
int circuit::local_overlap_num_calc(cell* theCell)
{
    int overNum = 0;

    rect theRect; // local_region
    theRect.xLL = max(die.xLL, theCell->init_x_coord - theCell->width*3);
    theRect.yLL = max(die.yLL, theCell->init_y_coord - theCell->height*3);
    theRect.xUR = min(die.xUR, theCell->init_x_coord + theCell->width*3);
    theRect.yUR = min(die.yUR, theCell->init_y_coord + theCell->height*3);
    
    vector<cell*> list = get_cells_from_boundary(&theRect);
    
    for(int i=0; i<(int)list.size(); i++) {
        for(int j=i+1; j<list.size(); j++) {
            rect rec1 = rect(list[i]->x_coord, list[i]->y_coord, list[i]->x_coord + list[i]->width, 
                        list[i]->y_coord + list[i]->height);
            rect rec2 = rect(list[j]->x_coord, list[j]->y_coord, list[j]->x_coord + list[j]->width, 
                        list[j]->y_coord + list[j]->height);

            if(check_overlap(rec1, rec2)) overNum++;
        }
    }
    return overNum;
}


double circuit::local_util_calc(cell* theCell)
{
    rect theRect; // local_region
    theRect.xLL = max(die.xLL, theCell->init_x_coord - theCell->width*3);
    theRect.yLL = max(die.yLL, theCell->init_y_coord - theCell->height*3);
    theRect.xUR = min(die.xUR, theCell->init_x_coord + theCell->width*3);
    theRect.yUR = min(die.yUR, theCell->init_y_coord + theCell->height*3);
    
    vector<cell*> list = get_cells_from_boundary(&theRect);

    double sumSize = 0;
    double rectArea = (theRect.xUR - theRect.xLL) * (theRect.yUR - theRect.yLL);
    
    for(int i=0; i<(int)list.size(); i++) {
        sumSize += list[i]->width * list[i]->height;
    }
    
    return 1.0 * sumSize / rectArea ;
}
*/

double circuit::cell_hpwl_calc(cell *theCell, string mode) {
    // Similar to circuit::HPWL(string mode)
    double hpwl = 0;
    double x_coord = 0;
    double y_coord = 0;

    vector<net *> cellNets;
    for (auto port: theCell->ports) {
        unsigned pidx = port.second;
        unsigned nidx = (&pins[pidx])->net;
        net *cellNet = &nets[nidx];
        cellNets.push_back(cellNet); // Need check (pin->net이 nets의 index 맞는지)
    }
    cout << "net size of cell " << theCell->name << ": " << cellNets.size() << endl;
    //////////////////////// EVERY CELLNETS has 0 element /////////////

    for (int i = 0; i < cellNets.size(); i++) {
        rect box;
        net *theNet = cellNets[i];

        pin *source = &pins[theNet->source];
        if (source->type == NONPIO_PIN) {
            cell *srcCell = &cells[source->owner];
            if (mode == "INIT") {
                x_coord = srcCell->init_x_coord;
                y_coord = srcCell->init_y_coord;
            } else {
                x_coord = srcCell->x_coord;
                y_coord = srcCell->y_coord;
            }
            box.xLL = box.xUR = x_coord + source->x_offset * DEFdist2Microns;
            box.yLL = box.yUR = y_coord + source->y_offset * DEFdist2Microns;
        } else {
            box.xLL = box.xUR = source->x_coord;
            box.yLL = box.yUR = source->y_coord;
        }

        for (int j = 0; j < theNet->sinks.size(); j++) {
            pin *sink = &pins[theNet->sinks[j]];
            if (sink->type == NONPIO_PIN) {
                cell *sinkCell = &cells[sink->owner];
                if (mode == "INIT") {
                    x_coord = sinkCell->init_x_coord;
                    y_coord = sinkCell->init_y_coord;
                } else {
                    x_coord = sinkCell->x_coord;
                    y_coord = sinkCell->y_coord;
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


double circuit::reward_calc(int movefailcnt, int gcellid, int pov) {
    //Disp calc start
    double avg_displacement = 0;
    double sum_displacement = 0;
    double max_displacement = 0;
    int count_displacement = 0;
    double violated_const = 0;

    int H1_count = 0;
    int H2_count = 0;
    int H3_count = 0;
    int H4_count = 0;

    int exist_H1 = 0;
    int exist_H2 = 0;
    int exist_H3 = 0;
    int exist_H4 = 0;

    double disp_H1 = 0;
    double disp_H2 = 0;
    double disp_H3 = 0;
    double disp_H4 = 0;

    field *theGcell = &Gcells[gcellid];

    //for(int i = 0; i < cells.size(); i++){
    //    cell* theCell = &cells[i];
    for (int i = 0; i < (int) theGcell->fieldCells.size(); i++) {
        cell *theCell = theGcell->fieldCells[i];
        double displacement =
                abs(theCell->init_x_coord - theCell->x_coord) + abs(theCell->init_y_coord - theCell->y_coord);
        sum_displacement += displacement;
        if (displacement > max_displacement) {
            max_displacement = displacement;
        }
        count_displacement++;

        if ((displacement / rowHeight) > max_disp_const)
            violated_const += static_cast<double>(displacement / rowHeight);

        if (theCell->height == 2000) {
            H1_count++;
            disp_H1 += displacement;
            exist_H1 = 1;
        }
        if (theCell->height == 4000) {
            H2_count++;
            disp_H2 += displacement;
            exist_H2 = 1;
        }
        if (theCell->height == 6000) {
            H3_count++;
            disp_H3 += displacement;
            exist_H3 = 1;
        }
        if (theCell->height == 8000) {
            H4_count++;
            disp_H4 += displacement;
            exist_H4 = 1;
        }

        // Not updated after placed (only effectedCells are updated in get_ovcells)
        //Pov += theCell->overlapNum; 

    }
    avg_displacement = sum_displacement / count_displacement;

    //Smm calc start
    double Smm;
    if ((violated_const / max_disp_const) > 1)
        Smm = 1 + static_cast<double>((violated_const / max_disp_const) * (max_displacement / (100 * rowHeight)));
    else
        Smm = 1 + static_cast<double>(max_displacement / (100 * rowHeight));

    //Sam calc start
    double Sam;
    double avg_disp_H1 = 0;
    double avg_disp_H2 = 0;
    double avg_disp_H3 = 0;
    double avg_disp_H4 = 0;
    if (exist_H1) avg_disp_H1 = static_cast<double>(disp_H1 / H1_count);
    if (exist_H2) avg_disp_H2 = static_cast<double>(disp_H2 / H2_count);
    if (exist_H3) avg_disp_H3 = static_cast<double>(disp_H3 / H3_count);
    if (exist_H4) avg_disp_H4 = static_cast<double>(disp_H4 / H4_count);

    Sam = static_cast<double>((avg_disp_H1 + avg_disp_H2 + avg_disp_H3 + avg_disp_H4) /
                              (exist_H1 + exist_H2 + exist_H3 + exist_H4)) / rowHeight;
    //cout << "Here" << endl;    


    double shpwl = std::max((HPWL("CUR") - HPWL("INIT")) / HPWL("INIT"), 0.0) * (1 + 0.2);
    //double shpwl = std::max((HPWL("CUR") - HPWL("INIT")) / HPWL("INIT"), 0.0) * (1 + std::max(calc_density_factor(8.0), 0.2));  

    //cout << " AVG_displacement : " << avg_displacement << endl;
    //cout << " SUM_displacement : " << sum_displacement << endl;
    //cout << " MAX_displacement : " << max_displacement << endl;
    //cout << " Smm              : " << Smm << endl;
    //cout << " Sam              : " << Sam << endl;
    //cout << " Shpwl            : " << shpwl << endl;
    //cout << " HPWL             : " << HPWL("CUR") << "    " << HPWL("INIT") << endl;

    if (movefailcnt > 0) cout << "[WARNING] Move fail count (x5.0 penalty): " << movefailcnt << endl;

    double S_total = Smm * Sam * (1 + shpwl) + 5.0 * movefailcnt;// + 0.01*pov;
    //cout << " [Gcell-" << gcellid << "] Stotal           : " << S_total << endl;
    return S_total;
}
