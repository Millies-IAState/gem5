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

#define UNSIGNED_BIT_COUNT 32
#define GHR_LENGTH UNSIGNED_BIT_COUNT
#define U8_BIT_COUNT 8
#define BYTES_PER_INT 4
#define U8_MAX 0xFF
#define WEIGHT_MAX 255
#define WEIGHT_MIN -255
#define WEIGHT_MOD (WEIGHT_MAX + WEIGHT_MAX + 1)
#define PHT_U_COUNT_OFFSET 1

SrnnBP::SrnnBP(const SrnnBPParams &params)
    : BPredUnit(params),
      localGHRSize(GHR_LENGTH), //Constant 32 Bit GHR for now
      localPHTSize(params.localPHTSize), //Number of PHT Rows (Each row has GHRSize w and u weights)
      localPHTUpdateWeight(params.localPHTUpdateWeight),
      PHT_w(localPHTSize),
      PHT_u(localPHTSize)
{
    DPRINTF(Fetch, "Flag1");
    if (!isPowerOf2(localGHRSize) || (localGHRSize > 32)) {
        fatal("Invalid GHR size! Must be power of 2 and 32 or less\n");
    }

    if (!isPowerOf2(localPHTSize)) {
        fatal("Invalid PHT size! Must be power of 2\n");
    }
    DPRINTF(Fetch, "Flag2");
    /** GHR must be a power of 2, so reducing the value by one should populate the lower bits.
     * As the PHT is indexed 0 to size - 1, this should be the valid mask.
    */
    PHT_index_mask = localGHRSize - 1;
    DPRINTF(Fetch, "Flag3");
    //Generate a random set of GHR Bits
    GHR = 0;
    for(size_t i = 0; i < BYTES_PER_INT; i++)
    {
        GHR = GHR | ((unsigned)rand() & U8_MAX);
        GHR = GHR << U8_BIT_COUNT;
    }
    GHR = GHR | ((unsigned)rand() & U8_MAX);

    DPRINTF(Fetch, "Flag4");

    //Initialize Random Weights
    for(size_t i = 0; i < localPHTSize; i++)
    {
        for (size_t j = 0; j < localGHRSize; j++)
        {
            //Setup rand seed
            srand(time(NULL));

            //Push a value anywhere between Weight Max and Min
            PHT_w[i].push_back((rand() % WEIGHT_MOD) - WEIGHT_MAX);
            PHT_u[i].push_back((rand() % WEIGHT_MOD) - WEIGHT_MAX);
        }
    }
    DPRINTF(Fetch, "Flag5");
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
    /** Create a new History Object to store */
    BPHistory *history = new BPHistory();
    bp_history = (void *)history;

    uint32_t local_predictor_idx = branch_addr | PHT_index_mask;
    std::vector<int32_t> weights = PHT_w[local_predictor_idx];
    std::vector<int32_t> uValues = PHT_u[local_predictor_idx];
    std::vector<int64_t> sValues(GHR_LENGTH,0);

    DPRINTF(Fetch, "Looking up index %#x\n",
            local_predictor_idx);

    //Initialize S Values
    for(size_t i = 0; i < GHR_LENGTH; i++)
    {
        //+ Weight on Taken, - Weight on not taken.
        if(((GHR >> i) & 0x01) > 0)
        {
            sValues[i] = weights[i];
        }
        else
        {
            sValues[i] = -1 * weights[i];
        }
    }

    int32_t sCount = GHR_LENGTH >> 1;
    int32_t uIndex = 0;
    while (sCount > 0)
    {
        for(size_t i = 0; i < sCount; i++)
        {
            sValues[i] = (sValues[(i << 1) + 1] + (uValues[uIndex] * sValues[i << 1]));
            uIndex = uIndex + 1;
        }
        sCount = sCount >> 1;
    }


    int64_t predictionValue = sValues[0];

    DPRINTF(Fetch, "prediction value (Y) is %lld.\n",
            predictionValue);

    taken = predictionValue > 0;

    history->prediction = taken;
    history->yValue = predictionValue;
    history->globalHistoryReg = GHR;

    unsigned takenValue = (taken) ? 1 : 0;
    //Update the GHR based on the taken value or not.
    GHR = (GHR << 1) | takenValue;

    return taken;
}

void
SrnnBP::update(ThreadID tid, Addr branch_addr, bool taken, void *bp_history,
                bool squashed, const StaticInstPtr & inst, Addr corrTarget)
{
    assert(bp_history == NULL);
    BPHistory *history = static_cast<BPHistory*>(bp_history);

    unsigned takenValue = (taken) ? 1 : 0;
    // If the taken value was invalid, restore the GHR and commit the correct value.
    if (squashed) {
        GHR = (history->globalHistoryReg << 1) | takenValue;
        return;
    }

    updatePHT(branch_addr, bp_history, taken);
    updateGHR(taken);
}

void
SrnnBP::updateGHR(bool taken)
{
    /** Update the GHR with the newest event, remove the oldest.. */
    unsigned takenValue = (taken) ? 1 : 0;
    GHR = (GHR << 1) | takenValue;
}

#define update_thresh 80
#define u_index 0
#define u_increment 2

void
SrnnBP::updatePHT(Addr pc, void *bp_history, bool actual)
{
    BPHistory *history = static_cast<BPHistory*>(bp_history);

    std::vector<int32_t> weights = PHT_w[pc | PHT_index_mask];
    std::vector<int32_t> uValues = PHT_u[pc | PHT_index_mask];

    if((abs(history->yValue) < update_thresh) || (history->prediction != actual))
    {
        for(size_t i = 0; i < localGHRSize; i++)
        {
            /** Check if the GHR and Prediction were equal, and we can improve the weights */
            if((((GHR >> i) & 0x01) > 0) == (history->prediction))
            {
                if(weights[i] < WEIGHT_MAX)
                {
                    /** If Equal, improve the weight*/
                    weights[i] = weights[i] + localPHTUpdateWeight;
                }
                /** There is one less u variable than w variables. */
                if((i < localGHRSize - PHT_U_COUNT_OFFSET) && (uValues[i] < WEIGHT_MAX))
                {
                    uValues[i] = uValues[i] + localPHTUpdateWeight;
                }
            }
            else //The prodiction was wrong
            {
                if(weights[i] > WEIGHT_MIN)
                {
                    /** If not equal, don't improve the weight.*/
                    weights[i] = weights[i] - localPHTUpdateWeight;
                }
                //Again, one fewer u value
                if((i < (localGHRSize - PHT_U_COUNT_OFFSET)) && (uValues[i] < WEIGHT_MAX))
                {
                    uValues[i] = 1;
                }
            }
        }
    }
}

void
SrnnBP::uncondBranch(ThreadID tid, Addr pc, void *&bp_history)
{
    //Unconditional Branches are ignored.
}

} // namespace branch_prediction
} // namespace gem5
