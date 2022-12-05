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
#include "debug/SrnnBPDB.hh"

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
#define WEIGHT_MAX 15
#define WEIGHT_MIN -15
#define WEIGHT_MOD (WEIGHT_MAX + WEIGHT_MAX + 1)
#define PHT_U_COUNT_OFFSET 1
#define update_thresh 10
#define u_index 0

SrnnBP::SrnnBP(const SrnnBPParams &params)
    : BPredUnit(params),
      localGHRSize(GHR_LENGTH), //Constant 32 Bit GHR for now
      localPHTSize(params.localPHTSize), //Number of PHT Rows (Each row has GHRSize w and u weights)
      localPHTUpdateWeight(params.localPHTUpdateWeight),
      PHT_w(localPHTSize),
      PHT_u(localPHTSize)
{
    DPRINTF(SrnnBPDB, "localGHRSize %u\r\n",localGHRSize);
    DPRINTF(SrnnBPDB, "localPHTSize %u\r\n",localPHTSize);
    DPRINTF(SrnnBPDB, "localPHTUpdateWeight %u\r\n", localPHTUpdateWeight);

    if (!isPowerOf2(localGHRSize) || (localGHRSize > 32)) {
        fatal("Invalid GHR size! Must be power of 2 and 32 or less\n");
    }

    if (!isPowerOf2(localPHTSize)) {
        fatal("Invalid PHT size! Must be power of 2\n");
    }
    firstPrediction = false;
    /** GHR must be a power of 2, so reducing the value by one should populate the lower bits.
     * As the PHT is indexed 0 to size - 1, this should be the valid mask.
    */
    PHT_index_mask = localPHTSize - 1;
    //Generate a random set of GHR Bits
    GHR = 0;
    for(size_t i = 0; i < BYTES_PER_INT; i++)
    {
        GHR = GHR | ((unsigned)rand() & U8_MAX);
        GHR = GHR << U8_BIT_COUNT;
    }
    GHR = GHR | ((unsigned)rand() & U8_MAX);

    //Initialize Random Weights
    for(size_t i = 0; i < localPHTSize; i++)
    {
        DPRINTF(SrnnBPDB, "Indexing PHT: %i\r\n",i);
        for (size_t j = 0; j < localGHRSize; j++)
        {
            //Setup rand seed
            srand(time(NULL));


            int32_t randW = (rand() % WEIGHT_MOD) - WEIGHT_MAX;
            int32_t randU = (rand() % WEIGHT_MOD) - WEIGHT_MAX;

            DPRINTF(SrnnBPDB, "Indexing W and U: %i\r\n",j);
            DPRINTF(SrnnBPDB, "Adding Random U Value %li\r\n",randU);
            DPRINTF(SrnnBPDB, "Adding Random W Value %li\r\n",randW);

            //Push a value anywhere between Weight Max and Min
            PHT_w[i].push_back(randW);
            PHT_u[i].push_back(randU);
        }
        DPRINTF(SrnnBPDB, "\r\n");
    }
}

void
SrnnBP::btbUpdate(ThreadID tid, Addr branch_addr, void * &bp_history)
{
    DPRINTF(SrnnBPDB, "Enter btbUpdate\r\n");
// Place holder for a function that is called to update predictor history when
// a BTB entry is invalid or not found.
    DPRINTF(SrnnBPDB, "Exiting btbUpdate\r\n");
}


bool
SrnnBP::lookup(ThreadID tid, Addr branch_addr, void * &bp_history)
{
    DPRINTF(SrnnBPDB, "Enter Lookup\r\n");
    bool taken;
    /** Create a new History Object to store */
    BPHistory *history = new BPHistory();
    bp_history = (void *)history;

    DPRINTF(SrnnBPDB, "Initializing Indexes and weights\r\n");
    uint64_t local_predictor_idx = (branch_addr >> 2) & PHT_index_mask;
    DPRINTF(SrnnBPDB, "Looking up index %#x\n",
            local_predictor_idx);

    std::vector<int32_t> weights = PHT_w[local_predictor_idx];
    DPRINTF(SrnnBPDB, "Completed Weights");
    std::vector<int32_t> uValues = PHT_u[local_predictor_idx];
    DPRINTF(SrnnBPDB, "Completed uValues");
    std::vector<int64_t> sValues(GHR_LENGTH,0);

    
    DPRINTF(SrnnBPDB, "Initializing SValues\r\n");
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
        DPRINTF(SrnnBPDB, "SValue Set: Index %li - Value: %lli\n",
        i,sValues[i]);
    }

    int32_t sCount = GHR_LENGTH >> 1;
    int32_t uIndex = 0;
    DPRINTF(SrnnBPDB, "Iterating internal nodes\r\n");
    while (sCount > 0)
    {
        DPRINTF(SrnnBPDB, "sCount Loop Value:%#x\n",
            sCount);
        for(size_t i = 0; i < sCount; i++)
        {
            int32_t index1 = i << 1;
            int32_t index2 = index1 + 1;
            DPRINTF(SrnnBPDB, "S Calc Inputs:\nIndex: %li Value:%lli\nIndex:%li Value:%lli\n",
            index1, sValues[index1],index2,sValues[index2]);
            DPRINTF(SrnnBPDB, "U Calc Input:\nIndex: %li Value:%lli\n",
            uIndex, uValues[uIndex]);
            sValues[i] = (sValues[index2] + (uValues[uIndex] * sValues[index1]));
            uIndex = uIndex + 1;
        }
        sCount = sCount >> 1;
    }


    int64_t predictionValue = sValues[0];

    DPRINTF(SrnnBPDB, "prediction value (Y) is %lli.\n",
            predictionValue);

    taken = predictionValue > 0;

    history->prediction = taken;
    history->yValue = predictionValue;
    history->globalHistoryReg = GHR;
    history->unconditionalBranch = false;

    unsigned takenValue = (taken) ? 1 : 0;
    //Update the GHR based on the taken value or not.
    GHR = (GHR << 1) | takenValue;
    DPRINTF(SrnnBPDB, "Exiting Lookup\r\n");
    firstPrediction = true;
    return taken;
    
}

void
SrnnBP::update(ThreadID tid, Addr branch_addr, bool taken, void *bp_history,
                bool squashed, const StaticInstPtr & inst, Addr corrTarget)
{
    DPRINTF(SrnnBPDB, "Entering Update\r\n");
    assert(bp_history);
    BPHistory *history = static_cast<BPHistory*>(bp_history);

    DPRINTF(SrnnBPDB, "Entering If\r\n");
    if (history->unconditionalBranch) {
        delete history;
        return;
    }

    DPRINTF(SrnnBPDB, "Entering squashed\r\n");
    // If the taken value was invalid, restore the GHR and commit the correct value.
    if (squashed) {
        unsigned takenValue = (taken) ? 1 : 0;
        GHR = (history->globalHistoryReg << 1) | takenValue;
            return;
    }

    updatePHT(branch_addr, bp_history, taken);
    updateGHR(taken);    
    DPRINTF(SrnnBPDB, "Exiting Update\r\n");
}

void
SrnnBP::updateGHR(bool taken)
{
    DPRINTF(SrnnBPDB, "Entering updateGHR\r\n");
    /** Update the GHR with the newest event, remove the oldest.. */
    unsigned takenValue = (taken) ? 1 : 0;
    GHR = (GHR << 1) | takenValue;
    DPRINTF(SrnnBPDB, "Exiting updateGHR\r\n");
}

void
SrnnBP::updatePHT(Addr pc, void *bp_history, bool actual)
{
    DPRINTF(SrnnBPDB, "Entering updatePHT\r\n");
    BPHistory *history = static_cast<BPHistory*>(bp_history);

    uint64_t local_predictor_idx = (pc >> 2) & PHT_index_mask;

    std::vector<int32_t> weights = PHT_w[local_predictor_idx];
    std::vector<int32_t> uValues = PHT_u[local_predictor_idx];

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
    DPRINTF(SrnnBPDB, "Exiting updatePHT\r\n");
}

void
SrnnBP::uncondBranch(ThreadID tid, Addr pc, void *&bp_history)
{
    //Unconditional Branches are ignored.
    DPRINTF(SrnnBPDB, "Entering unconditional Branch\r\n");
    BPHistory *history = new BPHistory();
    bp_history = (void *)history;

    /** Some default values for unconditional branch*/
    history->prediction = true;
    history->yValue = 1;
    history->globalHistoryReg = GHR;
    history->unconditionalBranch = true;

    DPRINTF(SrnnBPDB, "Exiting unconditional Branch\r\n");
}

} // namespace branch_prediction
} // namespace gem5
