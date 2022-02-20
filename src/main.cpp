/////////////////////////////////////////////////////////////////////////////
// Authors: SangGi Do(sanggido@unist.ac.kr), Mingyu Woo(mwoo@eng.ucsd.edu)
//          (respective Ph.D. advisors: Seokhyeong Kang, Andrew B. Kahng)
//
//          Original parsing structure was made by Myung-Chul Kim (IBM).
//
// BSD 3-Clause License
//
// Copyright (c) 2018, SangGi Do and Mingyu Woo
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#include "circuit.h"
#include <time.h>
#include "mymeasure.h"

using opendp::circuit;
using opendp::cell;
using opendp::row;
using opendp::pixel;
using opendp::rect;

using std::cout;
using std::endl;
using std::cerr;

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////////////////////////////////
    // arguments for debugging
    argc = 13;
    argv[0] = "../build/opendp";
    argv[1] = "-lef";
    argv[2] = "../bench/benchmarks/pci_bridge32_b_md1/tech.lef";
    argv[3] = "-lef";
    argv[4] = "../bench/benchmarks/pci_bridge32_b_md1/cells_modified.lef";
    argv[5] = "-def";
    argv[6] = "../bench/benchmarks/pci_bridge32_b_md1/placed.def";
    argv[7] = "-placement_constraints";
    argv[8] = "../bench/benchmarks/pci_bridge32_b_md1/placement.constraints";
    argv[9] = "-cpu";
    argv[10] = "1";
    argv[11] = "-output_def";
    argv[12] = "../output/pci_bridge32_b_md1_2022-02-18_13:07:04.def";
    // check arguments
    for (int i = 0; i < 3; ++i) {
        cout << endl;
    }
    for (int i = 0; i < argc; ++i) {
        cout << i << endl;
        cout << argv[i] << endl;
        cout << endl;
    }
    for (int i = 0; i < 3; ++i) {
        cout << endl;
    }
    //////////////////////////////////////////////////////////////////////////

    cout << "===================================================================="
            "======="
         << endl;
    cout << "   Open Source Mixed-Height Standard Cell Detail Placer < "
            "OpenDP_v1.0 >    "
         << endl;
    cout << "   Developers : SangGi Do, Mingyu Woo                               "
            "       "
         << endl;
    cout << "===================================================================="
            "======="
         << endl;

    CMeasure measure;
    measure.start_clock();

    circuit ckt;

    // READ input files - parser.cpp
    ckt.read_files(argc, argv);
    measure.stop_clock("Parser");

    ckt.simple_placement(measure);
    ckt.calc_density_factor(4);

    measure.stop_clock("All");
    ckt.write_def(ckt.out_def_name);

    measure.print_clock();

    // EVALUATION - utility.cpp
    ckt.evaluation();

    // CHECK LEGAL - check_legal.cpp
    ckt.check_legality();
    cout << " - - - - - < Program END > - - - - - " << endl;
    return 0;
}
