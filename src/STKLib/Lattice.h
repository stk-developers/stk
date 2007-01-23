/***************************************************************************
 *   copyright            : (C) 2004 by Lukas Burget,UPGM,FIT,VUT,Brno     *
 *   email                : burget@fit.vutbr.cz                            *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



//#############################################################################
//# STK includes                                                             ##
//#############################################################################
#include "Net.h"
#include "DecoderNetwork.h"


//#############################################################################
//# general includes                                                         ##
//#############################################################################


#ifndef STK_Lattice_h
#define STK_Lattice_h

namespace STK
{
  /** 
   * @brief 
   */
  class Lattice : public DecoderNetwork
  {
  public:
    /** 
     * @brief Builds the lattice from the network of word link records 
     * 
     * @param pWlr pointer to the "ending" node of the word link record
     * network
     */
    //void
    ////BuildFromWlr(WordLinkRecord* pWlr);
        


    /** 
     * @brief Performs forward-backward to get posterior probabilities
     * 
     * @return 
     */
    FLOAT
    ForwardBackward();
    
    /** 
     * @brief Free records with posterior probabilities
     * 
     * @return 
     */
    FLOAT
    Lattice::
    FreePosteriors();

    /** 
     * @brief Prunes the lattice based on the posterior probabilities
     * 
     * @param thresh 
     */
    void
    PosteriorPrune(const FLOAT& thresh);


  private:
    //static void 
    //AllocateNodesForWordLinkRecords(WordLinkRecord* pWlr, 
    //    NodeType*& rpNode);

    //static void 
    //EstablishLinks(WordLinkRecord* pWlr);
  }; 
  // class Lattice
  //***************************************************************************

};
// namespace STK
//*****************************************************************************

#endif // #ifndef STK_Lattice_h

//#############################################################################
//# EOF                                                                      ##
//#############################################################################

