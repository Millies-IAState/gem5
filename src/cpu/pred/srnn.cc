/*
 * Copyright (c) 2004-2006 The Regents of The University of Michigan
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "cpu/pred/srnn.hh"

#include "base/intmath.hh"
#include "base/logging.hh"
#include "base/trace.hh"
#include "debug/Fetch.hh"

#include <stdlib.h> // rand
#include <time.h> // time

namespace gem5
{

namespace branch_prediction
{

SrnnBP::SrnnBP(const SrnnBPParams &params)
    : BPredUnit(params),
      localGHRSize(params.localGHRSize),
      localPHTSize(params.localPHTSize),
      GHR(localGHRSize, 1)
{
    if (!isPowerOf2(localGHRSize)) {
        fatal("Invalid GHR size!\n");
    }

    if (!isPowerOf2(localPHTSize)) {
        fatal("Invalid PHT size!\n");
    }

    //Initialize Weights
    for (size_t i = 0; i < localPHTSize; i++)
    {
        //Setup rand seed
        srand(time(NULL));

        //Push a value anywhere between -15 and 15.
        PHT_w.push_back((rand() % 31) - 15);

        //U can be 1 or 0, use random to decide which to push for each loop.
        if((rand() % 2) == 1)
        {
            PHT_u.push_back(1);
        }
        else
        {
            PHT_u.push_back(-1);
        }
    }
    
}

void
SrnnBP::btbUpdate(ThreadID tid, Addr branch_addr, void * &bp_history)
{
// Place holder for a function that is called to update predictor history when
// a BTB entry is invalid or not found.
}


bool
SrnnBP::lookup(ThreadID tid, Addr branch_addr, void * &bp_history)
{
    bool taken;
    unsigned local_predictor_idx = getLocalIndex(branch_addr);

    DPRINTF(Fetch, "Looking up index %#x\n",
            local_predictor_idx);

    uint8_t counter_val = localCtrs[local_predictor_idx];

    DPRINTF(Fetch, "prediction is %i.\n",
            (int)counter_val);

    taken = getPrediction(counter_val);

    return taken;
}

void
SrnnBP::update(ThreadID tid, Addr branch_addr, bool taken, void *bp_history,
                bool squashed, const StaticInstPtr & inst, Addr corrTarget)
{
    assert(bp_history == NULL);

    updateGHR(taken);
    // No state to restore, and we do not update on the wrong
    // path.
    if (squashed) {
        return;
    }


    
}

void
SrnnBP::updateGHR(bool taken)
{
    //Remove the oldest element, and add on the newest.
    GHR.pop_back();
    if(taken)
    {
        GHR.push(1);
    }
    else
    {
        GHR.push(-1);
    }
}

void
SrnnBP::updatePHT()
{

}

inline
bool
SrnnBP::getPrediction(uint8_t &count)
{
    // Get the MSB of the count
    return (count >> (localCtrBits - 1));
}

inline
unsigned
SrnnBP::getLocalIndex(Addr &branch_addr)
{
    return (branch_addr >> instShiftAmt) & indexMask;
}

void
SrnnBP::uncondBranch(ThreadID tid, Addr pc, void *&bp_history)
{
}

} // namespace branch_prediction
} // namespace gem5
