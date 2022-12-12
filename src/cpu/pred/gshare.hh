#ifndef __CPU_PRED_GSHARE_HH__

#define __CPU_PRED_GSHARE_HH__

#include <vector>

#include "base/sat_counter.hh"
#include "base/types.hh"
#include "cpu/pred/bpred_unit.hh"
#include "params/GShareBP.hh"

namespace gem5
{

namespace branch_prediction
{

/**
 * Implements a local predictor that uses the PC to index into a table of
 * counters.  Note that any time a pointer to the bp_history is given, it
 * should be NULL using this predictor because it does not have any branch
 * predictor state that needs to be recorded or updated; the update can be
 * determined solely by the branch being taken or not taken.
 */
class GShareBP : public BPredUnit
{
  public:
    /**
     * Default branch predictor constructor.
     */
    GShareBP(const GShareBPParams &params);

    virtual void uncondBranch(ThreadID tid, Addr pc, void * &bp_history);

    /**
     * Returns if the branch should be taken or not, given a counter
     * value.
     * @param count The counter value.
     */
    inline bool getPrediction(uint8_t &count);

    /**
     * Looks up the given address in the branch predictor and returns
     * a true/false value as to whether it is taken.
     * @param branch_addr The address of the branch to look up.
     * @param bp_history Pointer to any bp history state.
     * @return Whether or not the branch is taken.
     */
    bool lookup(ThreadID tid, Addr branch_addr, void * &bp_history);

    /**
     * Updates the branch predictor to Not Taken if a BTB entry is
     * invalid or not found.
     * @param branch_addr The address of the branch to look up.
     * @param bp_history Pointer to any bp history state.
     * @return Whether or not the branch is taken.
     */
    void btbUpdate(ThreadID tid, Addr branch_addr, void * &bp_history);

    /**
     * Updates the branch predictor with the actual result of a branch.
     * @param branch_addr The address of the branch to update.
     * @param taken Whether or not the branch was taken.
     */
    void update(ThreadID tid, Addr branch_addr, bool taken, void *bp_history,
                bool squashed, const StaticInstPtr & inst, Addr corrTarget);

    void squash(ThreadID tid, void *bp_history)
    { } //assert(bp_history == NULL); }

  private:
    /**
     * Updates the global history register.
     *
     * @param taken True if branch was taken, else false
     */
    void updateGlobalHistory(ThreadID tid, bool taken);

    /**
     * Returns the local history index, given a branch address.
     * @param branch_addr The branch's PC address.
     */
    inline unsigned calcLocHistIdx(ThreadID tid, Addr &branch_addr);

    struct BPHistory {
        unsigned globalHistory;  /**< Global history register */
        bool prediction;            /**< True if taken, else false */
    };

    /** Number of counters in the local predictor. */
    unsigned localPredictorSize;

    /** Number of bits of the local predictor's counters. */
    const unsigned localCtrBits;

    /** Array of counters that make up the local predictor. */
    std::vector<SatCounter8> localCtrs;

    /** Global history register. Contains as much history as specified by
     *  globalHistoryBits. Actual number of bits used is determined by
     *  globalHistoryMask and choiceHistoryMask. */
    std::vector<unsigned> globalHistory;

    /** Number of bits for the global history. Determines maximum number of
        entries in global and choice predictor tables. */
    unsigned globalHistoryBits;

    /** Mask to apply to globalHistory to access global history table.
     *  Based on globalPredictorSize.*/
    unsigned globalHistoryMask;

    unsigned localPredictorBits;
    unsigned localPredictorMask;


    /** Thresholds for the counter value; above the threshold is taken,
     *  equal to or below the threshold is not taken.
     */
    unsigned localThreshold;

};

} // namespace branch_prediction
} // namespace gem5

#endif // __CPU_PRED_GSHARE_PRED_HH__
