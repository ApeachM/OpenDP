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
#include "time.h"

using namespace std;

//extern "C"
//{
circuit *ckt_new() {
    cout << "[CTYPES] New circuit!" << endl;
    return new circuit();
}

Agent *agent_new() {
    cout << "[CTYPES] New agent!" << endl;
    return new Agent();
}

void ckt_read_files(circuit *ck, int argc, char *argv[]) {
    ck->read_files(argc, argv);
    cout << "[CTYPES] Read files done!" << endl;
}

int ckt_region_assn(circuit *ck) {
    ck->region_assign();
    ck->avg_min2_search_distance();

    ck->field_unit = 5.0;
    ck->gcell_unit_calc(300);       // determine (ck->gcellUnit)
    ck->gcell_box_init();
    ck->field_box_init();

    cout << " Rtree initialized!" << endl;
    cout << "  - Gcells size: " << ck->Gcells.size() << endl;
    cout << "  - Fields size: " << ck->Fields.size() << endl;

    ck->Rtree_init_gcells_n_fields();
    ck->Rtree_init();
    ck->gcell_post_init();
    ck->get_fields_in_gcells();

    return (int) ck->Gcells.size();
}

void ckt_rtree_init(circuit *ck) {
    //cout << "[CTYPES] Rtree init!" << endl;
    //ck->gcell_box_init(ck->gcellUnit);
    //ck->field_box_init(5.0);
    ck->Rtree_init_gcells_n_fields();
    ck->Rtree_init();
    ck->get_cell_overlaps();
    //ck->gcell_post_init();              // get Gcell cells & Gcell sorting by cell-size
    //ck->get_fields_in_gcells();
    //cout << " Rtree initialized!" << endl;
    //cout << "  - Gcells size: " << ck->Gcells.size() << endl;
    //cout << "  - Fields size: " << ck->Fields.size() << endl;

    //return (int)ck->Gcells.size();
}

void ckt_simple_placement(circuit *ck)//, int gcellid)
{
    //ck->simple_placement(gcellid);
    ck->simple_placement();
    //cout << "Gcell-" << gcellid << " simple place done!" << endl;
    cout << "[CTYPES] simple place done!" << endl;
}

void ckt_write_def(circuit *ck) {
    //cout << "[CTYPES] Training episode done! Write DEF & Evaluation!" << endl;
    ck->write_def(ck->out_def_name);
    ck->evaluation();
    ck->check_legality();
    cout << " - - - - - < Program END > - - - - - - " << endl;
}

double ckt_agent_clear(circuit *ck, Agent *agent) {
    double GScore = ck->reward_calc(agent->moveFailCnt, agent->targetGcell, agent->Pov);
    clear_agent(ck, agent);
    return GScore;
}

void ckt_memory_clear(circuit *ck, Agent *agent) {
    ck->evaluation();
    //cout << "[CTYPES] Clear memory!" << endl;
    clear_memory(ck, agent);
}

// RL-placement
int ckt_state_init(circuit *ck, Agent *agent, int gcellid) {
    //cout << "[CTYPES] State init start!" << endl;
    vector<cell *> cells_ = ck->non_group_movable_cells(gcellid);

    cout << "  - (" << ck->Gcells[gcellid].xLL << " " << ck->Gcells[gcellid].yLL
         << ") (" << ck->Gcells[gcellid].xUR << " " << ck->Gcells[gcellid].yUR << ")" << endl;
    cout << "[CTYPES] Target cell size (gcell-" << gcellid << "): " << cells_.size() << endl;
    //cout << "[CTYPES] Num of cells in gcell-" << gcellid << ": " << ck->Gcells[gcellid].fieldCells.size() << endl;
    agent->state_init(ck, cells_, gcellid);
    //cout << "[CTYPES] Cell size: " << agent->state.cell_list.size() << endl;
    return (int) agent->state.cell_list.size();
}

double feature_get(Agent *agent, int cellidx, int fidx) {
    return agent->state.features[cellidx][fidx];
}

bool ckt_is_done(Agent *agent) {
    //cout << "[CTYPES] Isdone?" << endl;
    return agent->is_done();
}

int ckt_action(circuit *ck, Agent *agent, int tarID) {
//    cout << "[CTYPES] Action! Target cell index (id): " << tarID << endl;
    return agent->action(ck, tarID);
}

void ckt_feature_update(circuit *ck, Agent *agent, int tarID, int moveType) {
    //cout << "[CTYPES] Feature update!" << endl;
    agent->feature_update(ck, tarID, moveType);
}

int effected_cell_num(Agent *agent) {
    return (int) agent->effCells.size();
}

int effected_cell_sidx(Agent *agent, int ecellidx) {
    // sindex of effected cells
    return agent->effCells[ecellidx];
}

int effected_cell_id(Agent *agent, int cellidx) {
    int idx = agent->effCells[cellidx];
    cell *theCell = agent->state.cell_list[idx];
    return theCell->id;
}

double ckt_reward_calc(circuit *ck, Agent *agent) {
    //cout << "[CTYPES] Reward calculation!" << endl;
    //cout << "Gcell Stotal : " << ck->reward_calc(agent->moveFailCnt, agent->targetGcell, agent->Pov) << endl;
    //return -1.0*ck->REWARD;
    return ck->REWARD;
}

double ckt_Gcell_Stotal(circuit *ck, Agent *agent) {
    return ck->reward_calc(agent->moveFailCnt, agent->targetGcell, agent->Pov);
}
//}

int main(int argc, char *argv[]) {
    /*
     * example arguments
     * -tech_lef ../bench/benchmarks/edit_dist_a_md3/tech.lef -cell_lef ../bench/benchmarks/edit_dist_a_md3/cells_modified.lef -input_def ../bench/benchmarks/edit_dist_a_md3/placed.def -cpu 1 -placement_constraints ../bench/benchmarks/edit_dist_a_md3/placement.constraints -output_def ../output/edit_dist_a_md3_2022-05-04_22:05:59.def
     * Below main function is for virtual one of agent.py
     * */
    int featureNum = 9;
    auto start_time = clock();
    vector<vector<double>> timeset;

    // env = RLegalizer()
    circuit *ck = ckt_new();  // self.obj = self.lib.ckt_new()
    ckt_read_files(ck, argc, argv);  // self.ck.parse(sys.argv)
    ckt_region_assn(ck);  // self.ck.placeinit()
    Agent *agent = agent_new();  // self.ag = self.ck.agent()
    ckt_simple_placement(ck);  // self.ck.place()

    double pre_initial_time = double(clock() - start_time) / CLOCKS_PER_SEC;
    auto time_log = clock();

    vector<double> pre_station;
    pre_station.push_back(pre_initial_time);
    timeset.push_back(pre_station);


    int max_train_ep = 3;
    for (int episode = 0; episode < max_train_ep; ++episode) {
        vector<double> episode_times;
        // env.init()
        cout << "episode: " << episode << endl;
        ckt_rtree_init(ck);

        episode_times.push_back(double(clock() - time_log) / CLOCKS_PER_SEC);
        time_log = clock();

        // s = env.s_init(gcell) ->
        // self.ck.rl_init(self.ag, gcell)
        int gcell = 30;
        int cellNum = ckt_state_init(ck, agent, gcell);  // cellNum = self.lib.ckt_state_init(self.obj, agent, gcell)

        // feature extraction
        vector<vector<double>> cellFeatures;
        for (int cellIdx = 0; cellIdx < cellNum; ++cellIdx) {
            vector<double> cellFeature;
            for (int featureIdx = 0; featureIdx < featureNum; ++featureIdx) {
                double feature;
                feature = feature_get(agent, cellIdx, featureIdx);
                cellFeature.push_back(feature);
            }
            cellFeatures.push_back(cellFeature);
        }

        episode_times.push_back(double(clock() - time_log) / CLOCKS_PER_SEC);
        time_log = clock();

        // step
        for (auto & cellFeature : cellFeatures) {
            int targetID = cellFeature[0];
            int moveType = ckt_action(ck, agent, targetID);
            ckt_feature_update(ck, agent, targetID, moveType);
            ckt_reward_calc(ck, agent);
        }
        episode_times.push_back(double(clock() - time_log) / CLOCKS_PER_SEC);
        time_log = clock();

        double gScore = ckt_agent_clear(ck, agent);

        // episode end
        ckt_memory_clear(ck, agent);
        episode_times.push_back(double(clock() - time_log) / CLOCKS_PER_SEC);
        time_log = clock();

        timeset.push_back(episode_times);
    }
    // print time set
    cout << endl << endl;
    for (int i = 0; i < timeset.size(); ++i) {
        cout << "episode " << i << endl;
        for (int j = 0; j < timeset[i].size(); ++j) {
            if (timeset[i].size() == 1)
                cout << "preinit: ";
            else{
                if (j == 0) { cout << "rtree initial: "; }
                else if (j == 1) { cout << "feature extraction: "; }
                else if (j == 2) { cout << "Do actions: "; }
                else if (j == 3) { cout << "clearing: "; }
            }
            cout << timeset[i][j] << endl;

        }
        cout << endl;
    }

    delete ck;
    delete agent;
    cout << "train finished." << endl;
}
//}
