/*
 * Copyright (c) 2011, 2014 ARM Limited
 * All rights reserved
 *
 * The license below extends only to copyright in the software and shall
 * not be construed as granting a license to any other intellectual
 * property including but not limited to intellectual property relating
 * to a hardware implementation of the functionality of the software
 * licensed hereunder.  You may use the software subject to the license
 * terms below provided that you ensure that this notice is replicated
 * unmodified and in its entirety in all distributions of the software,
 * modified or unmodified, in source code or in binary form.
 *
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

#ifndef __CPU_PRED_2BIT_LOCAL_PRED_HH__
#define __CPU_PRED_2BIT_LOCAL_PRED_HH__

#include <vector>

#include "base/sat_counter.hh"
#include "base/types.hh"
#include "cpu/pred/bpred_unit.hh"
#include "params/SrnnBP.hh"

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
class SrnnBP : public BPredUnit
{
  public:
    /**
     * Default branch predictor constructor.
     */
    SrnnBP(const SrnnBPParams &params);

    /** Called on an unconditional branch execution. */
    virtual void uncondBranch(ThreadID tid, Addr pc, void * &bp_history) override;

    /**
     * Looks up the given address in the branch predictor and returns
     * a true/false value as to whether it is taken.
     * @param branch_addr The address of the branch to look up.
     * @param bp_history Pointer to any bp history state.
     * @return Whether or not the branch is taken.
     */
    bool lookup(ThreadID tid, Addr branch_addr, void * &bp_history) override;

    /**
     * Updates the branch predictor to Not Taken if a BTB entry is
     * invalid or not found.
     * @param branch_addr The address of the branch to look up.
     * @param bp_history Pointer to any bp history state.
     * @return Whether or not the branch is taken.
     */
    void btbUpdate(ThreadID tid, Addr branch_addr, void * &bp_history) override;

    /**
     * Updates the branch predictor with the actual result of a branch.
     * @param branch_addr The address of the branch to update.
     * @param taken Whether or not the branch was taken.
     */
    void update(ThreadID tid, Addr branch_addr, bool taken, void *bp_history,
                bool squashed, const StaticInstPtr & inst, Addr corrTarget) override;

    void squash(ThreadID tid, void *bp_history) override;

  protected:
  /**
   * PC Hash functions
   */
  static inline unsigned int hash1(unsigned int a)
  {
      a = (a ^ 0xdeadbeef) + (a<<4);
      a = a ^ (a>>10);
      a = a + (a<<7);
      a = a ^ (a>>13);
      return a;
  }
  static inline unsigned int hash2(unsigned int key)
  {
      int c2 = 0x27d4eb2d; // a prime or an odd constant
      key = (key ^ 61) ^ (key >> 16);
      key = key + (key << 3);
      key = key ^ (key >> 4);
      key = key * c2;
      key = key ^ (key >> 15);
      return key;
  }
  static inline unsigned int hash(unsigned int key, unsigned int i)
  {
      return hash2(key) * i + hash1(key);
  }
  static inline unsigned int hashPC(unsigned int pc, int pcshift)
  {
      if (pcshift < 0) {
          return hash(pc, -pcshift);
      } else if (pcshift < 11) {
          unsigned int x = pc;
          x ^= (pc >> pcshift);
          return x;
      } else {
          return pc >> (pcshift-11);
      }
  }

  private:
  class BPHistory
    {
      public:
      unsigned globalHistoryReg;

      Addr address;

      /** The Result of the SRNN network
      */
      int64_t yValue;
        
      /** Was the prediction taken or not taken
      * true: taken
      * false: not taken
      */
      bool prediction;

      /** Value indicating if the reason for the branch was unconditional */
      bool unconditionalBranch;

    };
    
    int32_t generateRandomInt();
    uint32_t generateRandomUnsignedInt();

    /** Updates the GHR Register*/
    void updateGHR(bool taken);

    /** Update the PHT */
    void updatePHT(Addr pc, void *bp_history, bool actual);

    /** Size of the GHR. */
    const unsigned localGHRSize;

    /** Size of the PHT */
    const unsigned localPHTSize;

    const unsigned localPHTUpdateWeight;

    /** The Global History Registers*/
    unsigned GHR;

    /** Bits available for the PHT */
    int32_t PHT_index_mask;

    /** PHT Weight Values*/
    std::vector<std::vector<int32_t>> PHT_w;

    /** PHT U Values*/
    std::vector<std::vector<int32_t>> PHT_u;

    /** Number of PHT Precision Bits*/
    unsigned localPHTBits;

    /** Max Weight */
    int32_t weightMax;

    /** Min Weight*/
    int32_t weightMin;

    int32_t updateThreshold;
};

} // namespace branch_prediction
} // namespace gem5

#endif // __CPU_PRED_2BIT_LOCAL_PRED_HH__
