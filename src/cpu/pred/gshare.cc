#include "cpu/pred/gshare.hh"

#include "base/intmath.hh"
#include "base/logging.hh"
#include "base/trace.hh"
#include "debug/Fetch.hh"

namespace gem5
{

namespace branch_prediction
{

GShareBP::GShareBP(const GShareBPParams &params)
    : BPredUnit(params),
      localPredictorSize(params.localPredictorSize), // set local predictor size
      localCtrBits(params.localCtrBits), // set localCtrBits
      localCtrs(localPredictorSize, SatCounter8(localCtrBits)),
      globalHistory(params.numThreads, 0),
      globalHistoryBits(ceilLog2(params.localPredictorSize)), // Set global instory bits to log2(localPredictorSize)
      localPredictorBits(ceilLog2(localPredictorSize))
{
    if (!isPowerOf2(localPredictorSize)) {
        fatal("Invalid local predictor size!\n");
    }

    DPRINTF(Fetch, "local predictor size: %i\n",
            localPredictorSize);

    DPRINTF(Fetch, "local counter bits: %i\n", localCtrBits);

    DPRINTF(Fetch, "Global History Bits $i\n", globalHistoryBits);

    DPRINTF(Fetch, "instruction shift amount: %i\n",
            instShiftAmt);

    globalHistoryMask = mask(this->globalHistoryBits);
    localPredictorMask = mask(localPredictorBits);
    localThreshold  = (1ULL << (localCtrBits  - 1)) - 1;
}

/**
* If a branch is not taken, because the BTB address is invalid or missing,
* this function sets the appropriate counter in the global and local
* predictors to not taken.
* @param inst_PC The PC to look up the local predictor.
* @param bp_history Pointer that will be set to an object that
* has the branch predictor state associated with the lookup.
*/
void
GShareBP::btbUpdate(ThreadID tid, Addr branch_addr, void * &bp_history)
{
    //Update Global History to Not Taken (clear LSB)
    globalHistory[tid] &= (globalHistoryMask & ~1ULL);
}


/**
 * Looks up a given PC in the BP to see if it is taken or not taken.
 * @param inst_PC The PC to look up.
 * @param bp_history Pointer that will be set to an object that
 * has the branch predictor state associated with the lookup.
 * @return Whether the branch is taken or not taken.
 */
bool
GShareBP::lookup(ThreadID tid, Addr branch_addr, void * &bp_history)
{
    // TODO
    unsigned tableIndex = calcLocHistIdx(tid, branch_addr);
    assert(tableIndex < localPredictorSize);

    bool prediction = getPrediction((uint8_t)localCtrs[tableIndex]);

    BPHistory *history = new BPHistory;
    history->globalHistory = globalHistory[tid];
    history->prediction = prediction;
    bp_history = static_cast<void*>(history);
    updateGlobalHistory(tid,true);
    return prediction  ;
}

/**
* Updates the BP with taken/not taken information.
* @param inst_PC The branch's PC that will be updated.
* @param taken Whether the branch was taken or not taken.
* @param bp_history Pointer to the branch predictor state that is
* associated with the branch lookup that is being updated.
* @param squashed Set to true when this function is called during a
* squash operation.
* @param inst Static instruction information
* @param corrTarget The resolved target of the branch (only needed
* for squashed branches)
*/
void
GShareBP::update(ThreadID tid, Addr branch_addr, bool taken, void *bp_history,
                bool squashed, const StaticInstPtr & inst, Addr corrTarget)
{
    if(bp_history) {
        BPHistory *history = static_cast<BPHistory *>(bp_history);
        unsigned tableIndex = calcLocHistIdx(tid, branch_addr);
        assert(tableIndex < localPredictorSize);

        if(taken) { localCtrs[tableIndex]++; }
        else { localCtrs[tableIndex]--; }

        if(squashed)
        {
            globalHistory[tid] = history->globalHistory << 1;
            if(taken) { globalHistory[tid] |= 1; }
            globalHistory[tid] &= globalHistoryMask;
        }

        delete history;
    }
}

/**
 * Returns if the branch should be taken or not, given a counter
 * value.
 * @param count The counter value.
 */
inline
bool
GShareBP::getPrediction(uint8_t &count) {
    return count > localThreshold;
}

/**
 * Unconditional Branch
 *
 * @param tid Thread ID
 * @param pc PC of branch
 * @param bp_history update branch predictor history
 */
void
GShareBP::uncondBranch(ThreadID tid, Addr pc, void *&bp_history)
{
    BPHistory *history = new BPHistory;
    history->globalHistory = globalHistory[tid];
    history->prediction = true;

    bp_history = static_cast<void*>(history);
    updateGlobalHistory(true);
    return;
}

/**
 * Updates the global history register.
 *
 * @param taken True if branch was taken, else false
 */
void
GShareBP::updateGlobalHistory(ThreadID tid, bool taken)
{
    // Update the global history register
    globalHistory[tid] = globalHistory[tid] << 1;
    if(taken) { globalHistory[tid] |= 1; }
    globalHistory[tid] &= globalHistoryMask;
}

/**
 * Returns the local history index, given a branch address.
 * @param branch_addr The branch's PC address.
 */
inline
unsigned
GShareBP::calcLocHistIdx(ThreadID tid, Addr &branch_addr)
{
    // Calculate the local predictor index
    // NOTE: Be sure to use instShiftAmt to shift the address to remove LSB
    // zeros
    unsigned returnValue = (branch_addr >> instShiftAmt);
    returnValue ^= globalHistory[tid];
    returnValue &= localPredictorMask;
    return returnValue;
}


} // namespace branch_prediction
} // namespace gem5
