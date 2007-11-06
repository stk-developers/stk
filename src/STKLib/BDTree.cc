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

#include "BDTree.h"
#include "Matrix.h"
#include "common.h"

#include <vector>
#include <iostream>


namespace STK
{
  // //***************************************************************************/
  // //***************************************************************************/
  // void
  // MapDistribution::
  // Dump(std::ostream& rStream, const std::string& rPrefix) const
  // {
  //   rStream << rPrefix;

  //   MapDistribution::Container::const_iterator i;
  //   for (i = mVec.begin(); i!=mVec.end(); ++i) {
  //     rStream << i->first << ": " << i->second << "\t";
  //   }
  //   rStream << std::endl;
  // }


  // //***************************************************************************/
  // //***************************************************************************/
  // const MapDistribution::ProbType
  // MapDistribution::
  // Entropy() const
  // {
  //   double e = 0;
  //   MapDistribution::Container::const_iterator i;

  //   for (i=mVec.begin(); i!=mVec.end(); ++i) {
  //     e -= i->second * log2(i->second);
  //   }
  //   
  //   return e;
  // }


  // //***************************************************************************/
  // //***************************************************************************/
  // const MapDistribution::ProbType
  // MapDistribution::
  // SplitEntropy(const MapDistribution& rD1, const MapDistribution& rD2)
  // {
  //   double w1 = rD1.mN / (rD1.mN + rD2.mN);
  //   double w2 = rD2.mN / (rD1.mN + rD2.mN);

  //   return w1*rD1.Entropy() + w2*rD2.Entropy();
  // }


  // //***************************************************************************/
  // //***************************************************************************/
  // void
  // MapDistribution::
  // ComputeFromCounts(const Data& rData)
  // {
  //   Data::const_iterator  i;
  //   MapDistribution::Container::iterator  j;
  //   double                n = 0;
  //     
  //   // set probs to 0
  //   for (j=mVec.begin(); j!=mVec.end(); ++j) {
  //     j->second = 0.0;
  //   }

  //   for (i=rData.begin(); i!=rData.end(); ++i) {
  //     n += *i;
  //   }

  //   for (i=rData.begin(), j=mVec.begin(); j!=mVec.end(); ++i, ++j) {
  //     j->second /= n;
  //   }
  // }


  // //***************************************************************************/
  // //***************************************************************************/
  // void
  // MapDistribution::
  // ComputeFromNGrams(const NGramSubset& rData)
  // {
  //   mVec.clear();
  //   mN = 0;

  //   // collect the counts
  //   NGramSubset::NGramContainer::const_iterator i;

  //   for (i=rData.mData.begin(); i!=rData.mData.end(); ++i) {
  //     NGram::TokenType* p_token = &((**i)[0]);
  //     NGram::ProbType   counts = (*i)->Counts();

  //     if (mVec.find(*p_token) == mVec.end()) {
  //       mVec[*p_token] = counts;
  //     }
  //     else {
  //       mVec[*p_token] += counts;
  //     }

  //     mN += counts;
  //   }
  //     
  //   // normalize counts to probs
  //   MapDistribution::Container::iterator  j;
  //   
  //   for (j=mVec.begin(); j!=mVec.end(); ++j) {
  //     j->second /= mN;
  //   }
  // }


  // //***************************************************************************/
  // //***************************************************************************/
  // MapDistribution::ProbType 
  // MapDistribution::
  // operator [] (const NGram::TokenType& rToken)
  // {
  //   if(mVec.find(rToken) != mVec.end()) {
  //     return mVec[rToken];
  //   }
  //   else {
  //     return 0.0;
  //   }
  // }


  // //***************************************************************************/
  // //***************************************************************************/
  // void
  // MapDistribution::
  // Smooth(const MapDistribution& rDistr, double r)
  // {
  //   MapDistribution::Container::iterator i_this;
  //   MapDistribution::Container::const_iterator i_parent;

  // } //Smooth(const MapDistribution& rDistr, double r)
  

  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Fix()
  {
    if (mN > 0) {
      double n=0.0;
      VecDistribution::Container::iterator i;

      for (i = mVec.begin(); i!=mVec.end(); ++i) {
        n += *i;
      }

      for (i = mVec.begin(); i!=mVec.end(); ++i) {
        *i /= n;
      }
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Dump(std::ostream& rStream, const std::string& rPrefix) const
  {
    rStream << rPrefix << "Data count: " << mN << std::endl;
    rStream << rPrefix;

    VecDistribution::Container::const_iterator i;
    for (i = mVec.begin(); i!=mVec.end(); ++i) {
      rStream << *i << "  ";
    }
    rStream << std::endl;
  }


  //***************************************************************************/
  //***************************************************************************/
  const VecDistribution::ProbType
  VecDistribution::
  Entropy() const
  {
    double e = 0;
    VecDistribution::Container::const_iterator i;

    for (i=mVec.begin(); i!=mVec.end(); ++i) {
      if (*i > 0) {
        e -= static_cast<double>(*i) * log(static_cast<double>(*i));
      }
    }
    
    return e;
  }


  //***************************************************************************/
  //***************************************************************************/
  const VecDistribution::ProbType
  VecDistribution::
  SplitEntropy(const VecDistribution& rD1, const VecDistribution& rD2)
  {
    double w1 = rD1.mN / (rD1.mN + rD2.mN);
    double w2 = rD2.mN / (rD1.mN + rD2.mN);

    return w1*rD1.Entropy() + w2*rD2.Entropy();
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  ComputeFromNGramSubsets(const NGramSubsets& rSubsets)
  {
    NGramSubsets::const_iterator i_subset;
    Container::iterator j;
    double mass = 0.0;
      
    // set probs to 0
    for (j=mVec.begin(); j!=mVec.end(); ++j) {
      *j = 0.0;
    }

    for (i_subset=rSubsets.begin(); i_subset!=rSubsets.end(); ++i_subset) {
      NGramSubset::NGramContainer::const_iterator i_ngram;

      // accumulate each dist and sum overal data mass
      for (i_ngram=i_subset->mData.begin(); 
           i_ngram!=i_subset->mData.end(); 
           ++i_ngram) 
      {
        NGram::TokenType* p_token = &((**i_ngram)[0]);
        NGram::ProbType   counts = (*i_ngram)->Counts();

        mVec[*p_token] += counts;
        mass           += counts;
      }
    }

    // normalize
    for (j=mVec.begin(); j!=mVec.end(); ++j) {
      *j = *j / mass;
    }

    mN = mass;
  }


  //***************************************************************************/
  //***************************************************************************/
  // void
  // VecDistribution::
  // ComputeFromNGrams(const NGramSubset& rData)
  // {
  //   mVec.clear();
  //   mN = 0;

  //   // collect the counts
  //   NGramSubset::NGramContainer::const_iterator i;

  //   for (i=rData.mData.begin(); i!=rData.mData.end(); ++i) {
  //     NGram::TokenType* p_token = &((**i)[0]);
  //     NGram::ProbType   counts = (*i)->Counts();

  //     mVec[*p_token]  += counts;
  //     mN              += counts;
  //   }
  //     
  //   // normalize counts to probs
  //   VecDistribution::Container::iterator  j;
  //   
  //   for (j=mVec.begin(); j!=mVec.end(); ++j) {
  //     *j /= mN;
  //   }
  // }


  //***************************************************************************/
  //***************************************************************************/
  VecDistribution::ProbType 
  VecDistribution::
  operator [] (const NGram::TokenType& rToken)
  {
    assert (rToken < mVec.size());
    return mVec[rToken];
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Merge(const VecDistribution& rDistr)
  {
    double norm = 0.0;
    VecDistribution::Container::iterator i_this = mVec.begin();
    VecDistribution::Container::const_iterator i_parent = rDistr.mVec.begin();

    assert(mVec.size() == rDistr.mVec.size());

    while (i_this != mVec.end()) {
      *i_this = mN * (*i_this) + rDistr.mN * (*i_parent);
      norm += *i_this;

      ++i_this;
      ++i_parent;
    } // while (i_this != mVec.end()) {

    // normalize
    if (norm > 0) {
      for (i_this = mVec.begin(); i_this != mVec.end(); ++i_this) {
        *i_this /= norm;
      }
    }

    mN += rDistr.mN;
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Reset() 
  {
    VecDistribution::Container::iterator i_this;

    for (i_this = mVec.begin(); i_this != mVec.end(); ++i_this) {
      *i_this = 0.0;
    }

    mN = 0.0;
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Smooth(const VecDistribution& rDistr, double r)
  {
    double norm   = 0; //DBL_EPSILON;
    double new_mn = 0; //DBL_EPSILON;
    VecDistribution::Container::iterator        i_this   = mVec.begin();
    VecDistribution::Container::const_iterator  i_parent = rDistr.mVec.begin();

    assert(mVec.size() == rDistr.mVec.size());

    while (i_this != mVec.end()) {
      // compute weight
      double b =  *i_this*mN / (*i_this*mN + r) ;

      *i_this = b * (*i_this) + (1-b) * (*i_parent);
      norm    += *i_this;
      new_mn  += b * (*i_this * mN) + (1-b) * (*i_parent * rDistr.mN);

      ++i_this;
      ++i_parent;
    } // while (i_this != mVec.end() && i_parent != rDistr.mVec.end()) {

    if (norm > 0) {
      // normalize
      for (i_this = mVec.begin(); i_this != mVec.end(); ++i_this) {
        *i_this /= norm;
      }
    }
    //mN = new_mn;
  } //Smooth(const VecDistribution& rDistr, double r)
  

  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  MapAdapt(const VecDistribution& rDistr, double r)
  {
    double norm   = 0; //DBL_EPSILON;
    double new_mn = 0; //DBL_EPSILON;
    VecDistribution::Container::iterator        i_this   = mVec.begin();
    VecDistribution::Container::const_iterator  i_parent = rDistr.mVec.begin();

    assert(mVec.size() == rDistr.mVec.size());

    while (i_this != mVec.end()) {
      // compute weight
      double b =  mN / (mN + r) ;

      *i_this = b * (*i_this) + (1-b) * (*i_parent);
      norm    += *i_this;

      ++i_this;
      ++i_parent;
    } // while (i_this != mVec.end() && i_parent != rDistr.mVec.end()) {

    if (norm > 0) {
      // normalize
      for (i_this = mVec.begin(); i_this != mVec.end(); ++i_this) {
        *i_this /= norm;
      }
    }
  } //MapAdapt(const VecDistribution& rDistr, double r)
  

  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  BDTree(NGramSubsets& rData, BDTreeBuildTraits& rTraits,
      BDTree* pParent, const std::string& rPrefix)
  : mDist(rData[0].Parent().pTargetTable()->Size()),
    mpBackoffDist(NULL)
  {
    BuildFromNGrams(rData, rTraits, pParent, rPrefix);
  }


  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  BDTree(const BDTree& rOrig)
  : mDist(rOrig.mDist)
  {
    if (!rOrig.IsLeaf()) {
      mpTree0     = rOrig.mpTree0->Clone();
      mpTree1     = rOrig.mpTree1->Clone();
      mpQuestion  = rOrig.mpQuestion->Clone();
    }
    else {
      mpTree0     = NULL;
      mpTree1     = NULL;
      mpQuestion  = NULL;
    }
    
    if (NULL == rOrig.mpBackoffDist) {
      mpBackoffDist = NULL;
    }
    else {
      mpBackoffDist = new VecDistribution(*(rOrig.mpBackoffDist));
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  BDTree(std::istream& rStream, BDTreeHeader& rHeader)
  : mDist(rHeader.mVocabSize)
  {
    Read(rStream, rHeader);
  }


  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  ~BDTree()
  {
    if (NULL != mpTree0) {
      delete mpTree0;
    }
    if (NULL != mpTree1) {
      delete mpTree1;
    }
    if (NULL != mpQuestion) {
      delete mpQuestion;
    }
    if (NULL != mpBackoffDist) {
      delete mpBackoffDist;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  BuildFromNGrams(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits,
      BDTree* pParent, const std::string& rPrefix)
  {
    double    total_entropy;
    double    split_entropy;
    double    entropy_red;

    double    total_mmi;
    double    split_mmi;
    double    mmi_inc;
    BQuestion* p_new_question = NULL;
    int       depth = rTraits.mCurDepth;

    // compute distribution
    //mDist.ComputeFromNGrams(rNGrams[0].mData.begin(), rNGrams[0].mData.end());
    mDist.ComputeFromNGramSubsets(rNGrams);

    if (NULL != pParent ) {
      if (rTraits.mSmoothR > 0) {
        mDist.Smooth(pParent->mDist, rTraits.mSmoothR);
      }
    }


    if (rTraits.mVerbosity > 0) {
      std::cout << "BDTree::Info      " << rPrefix << "Data mass: " << rNGrams.Mass() << std::endl;
      std::cout << "BDTree::Info      " << rPrefix << "Token count: " << rNGrams.TokenCount() << std::endl;
      std::cout << "BDTree::Info      " << rPrefix << "Depth:     " << depth << std::endl;
    }

    // check for max depth criterion
    if (rTraits.mMaxDepth <= rTraits.mCurDepth) {
      if (rTraits.mVerbosity > 0) {
        std::cout << "BDTree::Info      " << rPrefix << "Leaf due to maximum depth criterion" << std::endl;
      }

      mpQuestion = NULL;
      mpTree0    = NULL;
      mpTree1    = NULL;
      return;
    }

    // check for data sufficiency
    if (mDist.Counts() < rTraits.mMinInData) {
      if (rTraits.mVerbosity > 0) {
        std::cout << "BDTree::Info      " << rPrefix << "Leaf due to minimum input data criterion" << std::endl;
      }

      mpQuestion = NULL;
      mpTree0    = NULL;
      mpTree1    = NULL;
      return;
    }

    if (rTraits.mMMItoggle) {
      // compute mmi
      total_mmi       = rNGrams.MMI(rTraits.mMMIEta);
      p_new_question  = FindSimpleQuestion_MMI(rNGrams, rTraits, &split_mmi);
      mmi_inc         = split_mmi - total_mmi;

      if (NULL != p_new_question) {
        p_new_question->Dump(std::cout, std::string("BDTree::Info      ")+rPrefix);
        // compute entropy
        total_entropy = rNGrams.ParallelEntropy();
        split_entropy = rNGrams.ParallelSplitEntropy(*p_new_question, NULL, NULL);
        entropy_red   = total_entropy - split_entropy;

        std::cout << "BDTree::Info      " << rPrefix << "Data Entropy:      " << total_entropy << std::endl;
        std::cout << "BDTree::Info      " << rPrefix << "Split Entropy:     " << split_entropy << std::endl;
        std::cout << "BDTree::Info      " << rPrefix << "Entropy reduction: " << entropy_red << std::endl;
        std::cout << "BDTree::Info      " << rPrefix << "Data MMI:           " << total_mmi << std::endl;
        std::cout << "BDTree::Info      " << rPrefix << "Data split MMI:     " << split_mmi << std::endl;
        std::cout << "BDTree::Info      " << rPrefix << "MMI Increase:       " << mmi_inc << std::endl;

        if (isnan(entropy_red)) {
          throw std::runtime_error("Entropy is unexpectedly NaN");
        }
      }

      if (rTraits.mVerbosity > 0) {
      }
    }
    else {
      // compute entropy
      total_entropy = rNGrams.ParallelEntropy();
      // find best predictor and return split entropy
      p_new_question = FindSimpleQuestion(rNGrams, rTraits, &split_entropy);

      entropy_red = total_entropy - split_entropy;

      if (rTraits.mVerbosity > 0) {
        if (entropy_red == 0) {
          std::cout << "BDTree::Info      " << rPrefix << "No predictor found " << std:: endl;
          std::cout << "BDTree::Info      " << rPrefix << "Data Entropy:      " << total_entropy << std::endl;
        }
        else {
          p_new_question->Dump(std::cout, std::string("BDTree::Info      ")+rPrefix);
          std::cout << "BDTree::Info      " << rPrefix << "Data Entropy:      " << total_entropy << std::endl;
          std::cout << "BDTree::Info      " << rPrefix << "Split Entropy:     " << split_entropy << std::endl;
          std::cout << "BDTree::Info      " << rPrefix << "Entropy reduction: " << entropy_red << std::endl;
        }
      }
    }

    if (rTraits.mMMItoggle && mmi_inc < rTraits.mMinMMIIncrease) {
      if (rTraits.mVerbosity > 0) {
        std::cout << "BDTree::Info      " << rPrefix << "Leaf due to minimum MMI increase criterion" << std::endl;
      }

      if (NULL != p_new_question) {
        delete p_new_question;
      }

      mpQuestion = NULL;
      mpTree0    = NULL;
      mpTree1    = NULL;
    }
    else if (entropy_red < rTraits.mMinReduction) {
      if (rTraits.mVerbosity > 0) {
        std::cout << "BDTree::Info      " << rPrefix << "Leaf due to minimum entropy reduction criterion" << std::endl;
      }

      if (NULL != p_new_question) {
        delete p_new_question;
      }

      mpQuestion = NULL;
      mpTree0    = NULL;
      mpTree1    = NULL;
    }
    // if entropy reduction && MI increase still big enough, build each subtree
    else {
      NGramSubsets* p_split_data0 = new NGramSubsets;
      NGramSubsets* p_split_data1 = new NGramSubsets;

      mpQuestion  = p_new_question;
      rNGrams.Split(*mpQuestion, *p_split_data0, *p_split_data1);

      rTraits.mCurDepth = depth + 1;
      mpTree0     = new BDTree(*p_split_data0, rTraits, this, rPrefix + "-  ");
      rTraits.mCurDepth = depth + 1;
      mpTree1     = new BDTree(*p_split_data1, rTraits, this, rPrefix + "+  ");

      delete p_split_data0;
      delete p_split_data1;
    } // else declare leaf
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  GetInfo(BDTreeInfo& rInfo, const BDTree* pParent) const
  {
    rInfo.mTotalNodes++;

    if (!IsLeaf()) {
      mpTree0->GetInfo(rInfo, this);
      mpTree1->GetInfo(rInfo, this);

      if (NULL == pParent) {
        rInfo.mAverageEntropy /= rInfo.mTotalData;
      }
    }
    else {
      rInfo.mTotalLeaves++;
      rInfo.mTotalData += mDist.Counts();
      rInfo.mAverageEntropy += mDist.Counts() * mDist.Entropy();
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion*
  BDTree::
  FindSimpleQuestion(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, double* pSplitEnt) 
  {
    double          e = 0.0;
    double          best_e = -1.0;
    BSetQuestion*   p_best_question = NULL;


    // go through all predictors and find subset
    for (int i=1; i<rTraits.mOrder; ++i) {
      BSetQuestion* p_tmp_question = FindSubset_Greedy(rNGrams, rTraits, &e, i);

      if (e < best_e || best_e < 0) {
        if (NULL != p_best_question) {
          delete p_best_question;
        }
        best_e = e;
        p_best_question = p_tmp_question;
      } 
    }

    if (NULL != pSplitEnt ) {
      *pSplitEnt = best_e;
    }

    return p_best_question;
  }


  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion*
  BDTree::
  FindSimpleQuestion_MMI(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, double* pSplitMMI) 
  {
    double          mmi = -10000;
    double          best_mmi = mmi - 1.0;
    BSetQuestion*   p_best_question = NULL;


    // go through all predictors and find subset
    for (int i=1; i<rTraits.mOrder; ++i) {
      BSetQuestion* p_tmp_question = FindSubset_Greedy_MMI(rNGrams, rTraits, &mmi, i);

      if (mmi > best_mmi) {
        if (NULL != p_best_question) {
          delete p_best_question;
        }
        best_mmi = mmi;
        p_best_question = p_tmp_question;
      } 
    }

    if (NULL != pSplitMMI ) {
      *pSplitMMI = best_mmi;
    }

    return p_best_question;
  }


  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion*
  BDTree::
  FindSubset_Greedy(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, double* pSplitEnt,
      int pred) 
  {
    bool    inserted  = true;
    bool    deleted   = true;
    double  e         = rNGrams.ParallelEntropy();
    int     vocab_size = rNGrams[0].Parent().pTargetTable()->Size();
    int     i_vocab;
    int     i_best_change;
    double  mass0;
    double  mass1;
    BSetQuestion* p_new_question = new BSetQuestion(pred, vocab_size);

    
    // continue probing if any change was made in the previous iteration
    while (inserted || deleted) {
      inserted      = false;
      deleted       = false;

      // probe insertions
      for (i_vocab=0; i_vocab<vocab_size; ++i_vocab) {
        if (!p_new_question->EvalRawToken(i_vocab)) {
          p_new_question->Set(i_vocab);
          double this_e = rNGrams.ParallelSplitEntropy(*p_new_question, &mass0, &mass1);

          if (this_e < e && mass0 >= rTraits.mMinInData && 
              mass1 >= rTraits.mMinInData) 
          {
            e             = this_e;
            i_best_change = i_vocab;
            inserted      = true;
          }

          p_new_question->Unset(i_vocab);
        }
      }

      // commit insertion
      if (inserted) {
        p_new_question->Set(i_best_change);
      }

      // probe deletions
      for (i_vocab=0; i_vocab < vocab_size; ++i_vocab) {
        // only inspect positive questions
        if (p_new_question->EvalRawToken(i_vocab)) {
          p_new_question->Unset(i_vocab);
          double this_e = rNGrams.ParallelSplitEntropy(*p_new_question, &mass0, &mass1);

          if (this_e < e && mass0 >= rTraits.mMinInData && 
              mass1 >= rTraits.mMinInData) 
          {
            e             = this_e;
            i_best_change = i_vocab;
            deleted       = true;
          }

          p_new_question->Set(i_vocab);
        }
      }

      // commit deletion
      if (deleted) {
        p_new_question->Unset(i_best_change);
      }
    }

    // update entropy if possible
    if (NULL != pSplitEnt) {
      *pSplitEnt = e;
    }

    return p_new_question;
  }


  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion*
  BDTree::
  FindSubset_Greedy_MMI(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, 
      double* pSplitMMI, int pred) 
  {
    double  mass0;
    double  mass1;
    bool    inserted  = true;
    bool    deleted   = true;
    double  mmi       = rNGrams.MMI(rTraits.mMMIEta);
    size_t  vocab_size = rNGrams[0].Parent().pTargetTable()->Size();
    BSetQuestion* p_new_question = 
      new BSetQuestion(pred, vocab_size);

    size_t i_vocab;
    size_t i_best_change;
    

    // continue probing if any change was made in the previous iteration
    while (inserted || deleted) {
      inserted      = false;
      deleted       = false;

      // probe insertions
      for (i_vocab=0; i_vocab!=vocab_size; ++i_vocab) {
        if (!p_new_question->EvalRawToken(i_vocab)) {
          p_new_question->Set(i_vocab);
          double this_mmi = rNGrams.SplitMMI(*p_new_question, rTraits.mMMIAlpha, 
              rTraits.mMMIEta, &mass0, &mass1);

          if (this_mmi > mmi && mass0 >= rTraits.mMinInData && 
              mass1 >= rTraits.mMinInData) 
          {
            mmi           = this_mmi;
            i_best_change = i_vocab;
            inserted      = true;
          }

          p_new_question->Unset(i_vocab);
        }
      }

      // commit insertion
      if (inserted) {
        p_new_question->Set(i_best_change);
      }

      // probe deletions
      for (i_vocab=0; i_vocab!=vocab_size; ++i_vocab) {
        // only inspect positive questions
        if (p_new_question->EvalRawToken(i_vocab)) {
          p_new_question->Unset(i_vocab);
          double this_mmi = rNGrams.SplitMMI(*p_new_question, rTraits.mMMIAlpha, 
              rTraits.mMMIEta, &mass0, &mass1);

          if (this_mmi > mmi && mass0 >= rTraits.mMinInData && 
              mass1 >= rTraits.mMinInData) 
          {
            mmi           = this_mmi;
            i_best_change = i_vocab;
            deleted       = true;
          }

          p_new_question->Set(i_vocab);
        }
      }

      // commit deletion
      if (deleted) {
        p_new_question->Unset(i_best_change);
      }
    }

    // update entropy if possible
    if (NULL != pSplitMMI) {
      *pSplitMMI = mmi;
    }

    return p_new_question;
  }


  //***************************************************************************/
  //***************************************************************************/
  double
  BDTree::
  ScoreNGram(const NGram& rNGram)
  {
    BDTree* p_node = this;

    while (!p_node->IsLeaf()) {
      p_node = p_node->mpQuestion->Eval(rNGram) ? p_node->mpTree1 : p_node->mpTree0;
    } 

    if (NULL != p_node->mpBackoffDist) {
      return (*(p_node->mpBackoffDist))[rNGram[0]];
    }
    else {
      return p_node->mDist[rNGram[0]];
    }
  } // ScoreNGram(const NGram& rNGram)


  //***************************************************************************/
  //***************************************************************************/
  double
  BDTree::
  ScoreNGramSubset(const NGramSubset& rNGrams)
  {
    NGramSubset::NGramContainer::const_iterator i;
    double  score = 0.0;

    for (i = rNGrams.mData.begin(); i != rNGrams.mData.end(); ++i) {
      double n_gram_score = ScoreNGram(**i);

      // we don't score the n-grams, whose prob is 0
      if (n_gram_score > 0) {
        score += (*i)->Counts() * log(n_gram_score);
      }
    }

    return score;
  } // ScoreNGramSubset(const NGramSubset& rNGrams, double r)


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  Dump(std::ostream& rStream, const std::string& rPrefix) const
  {
    if (IsLeaf()) {
      if (NULL != mpBackoffDist) {
        mpBackoffDist->Dump(rStream, rPrefix);
      }
      else {
        mDist.Dump(rStream, rPrefix);
      }

    }
    else {
      std::string new_pref;
      mpQuestion->Dump(rStream, rPrefix);

      //mDist.Dump(rStream, rPrefix);

      new_pref = rPrefix + "-  ";
      mpTree0->Dump(rStream, new_pref);

      new_pref = rPrefix + "+  ";
      mpTree1->Dump(rStream, new_pref);
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  RecomputeDists()
  {
    if (!IsLeaf()) {
      mpTree0->RecomputeDists();
      mpTree1->RecomputeDists();
      mDist.Reset();
      mDist.Merge(mpTree0->mDist);
      mDist.Merge(mpTree1->mDist);
      //mDist.Dump(std::cout, "xx");
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  ComputeBackoffDists(BDTree* pParent, double r)
  {
    mpBackoffDist = new VecDistribution(mDist);

    if (NULL != pParent && r > 0) {
      mpBackoffDist->Smooth(*(pParent->mpBackoffDist), r);
    }

    if (!IsLeaf()) {
      mpTree0->ComputeBackoffDists(this, r);
      mpTree1->ComputeBackoffDists(this, r);
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  CollectCounts(const NGramSubset& rNGrams)
  {
    NGramSubset::NGramContainer::const_iterator it;

    // go through all data
    for (it = rNGrams.mData.begin(); it != rNGrams.mData.end(); ++it) {
      BDTree* p_node  = this;
      NGram&  r_ngram = **it;
      
      while (!p_node->IsLeaf()) {
        p_node = p_node->mpQuestion->Eval(r_ngram) ? p_node->mpTree1 : p_node->mpTree0;
      } 

      //std::cout << r_ngram[0] << std::endl;
      p_node->mDist.mVec[r_ngram[0]]  += r_ngram.Counts();
      p_node->mDist.mN                += r_ngram.Counts();
      //p_node->mDist.Dump(std::cout, "");
    }
  }
  
  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  AdaptLeafDists(const NGramSubset& rNGrams, BDTreeBuildTraits& rTraits)
  {
    double r = rTraits.mAdaptR;

    NGramSubset::NGramContainer::const_iterator it;

    std::vector<BDTree*> leaf_map;
    std::vector<BDTree*>::iterator i_leaf;
    GetLeaves(leaf_map);

    // go through all data
    for (it = rNGrams.mData.begin(); it != rNGrams.mData.end(); ++it) {
      BDTree* p_node  = this;
      NGram&  r_ngram = **it;
      
      while (!p_node->IsLeaf()) {
        p_node = p_node->mpQuestion->Eval(r_ngram) ? p_node->mpTree1 : p_node->mpTree0;
      } 

      // if mpBackoff is not initialized yet
      if (NULL == p_node->mpBackoffDist) {
        p_node->mpBackoffDist = new VecDistribution(p_node->mDist);
        p_node->mDist.mN      = 0.0;
        p_node->mDist.Reset();
      }

      p_node->mDist.mVec[r_ngram[0]]  += r_ngram.Counts();
      p_node->mDist.mN                += r_ngram.Counts();
    }

    if (rTraits.mMapAdapt) {
      // recompute all distributions
      for (i_leaf = leaf_map.begin(); i_leaf != leaf_map.end(); ++i_leaf) {
        (*i_leaf)->mDist.Fix();
        if (r > 0) {
          (*i_leaf)->mDist.MapAdapt(*((**i_leaf).mpBackoffDist), r);
        }
      }
    }
    else {
      // recompute all distributions
      for (i_leaf = leaf_map.begin(); i_leaf != leaf_map.end(); ++i_leaf) {
        (*i_leaf)->mDist.Fix();
        if (r > 0) {
          (*i_leaf)->mDist.Smooth(*((**i_leaf).mpBackoffDist), r);
        }
      }
    }
  }
  
  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  ResetLeaves()
  {
    if (!IsLeaf()) {
      mpTree0->ResetLeaves();
      mpTree1->ResetLeaves();
    }
    else {
      mDist.Reset();
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  int
  BDTree::
  GetLeaves(std::vector<BDTree*>& rCollection)
  {
    int n = 0;
    if (IsLeaf()) {
      rCollection.insert(rCollection.end(), this);
      return 1;
    }
    else {
      n += mpTree0->GetLeaves(rCollection);
      n += mpTree1->GetLeaves(rCollection);
      return n;
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  FillLeafSupervector(BasicVector<double>& rVector, bool backoff,
      bool includeCounts)
  {
    size_t vec_size;
    size_t vocab_size;
    size_t n_leaves;
    size_t i = 0;

    std::vector<BDTree*> leaves;
    std::vector<BDTree*>::iterator i_leaf;
    VecDistribution::Container::iterator i_prob;

    // compute sizes
    n_leaves    = GetLeaves(leaves);
    vocab_size  = leaves.front()->mDist.mVec.size();
    vec_size    = includeCounts ? n_leaves * vocab_size + n_leaves 
                                : n_leaves * vocab_size;

    // prepare new vector
    rVector.Destroy();
    rVector.Init(vec_size);
    rVector.Clear();

    // stack the leaves together and copy to the new supervector
    for (i_leaf = leaves.begin(); i_leaf != leaves.end(); ++i_leaf) {
      VecDistribution::Container::iterator i_prob_begin = backoff ? 
        (**i_leaf).mpBackoffDist->mVec.begin() : (**i_leaf).mDist.mVec.begin();

      VecDistribution::Container::iterator i_prob_end   = backoff ? 
        (**i_leaf).mpBackoffDist->mVec.end() : (**i_leaf).mDist.mVec.end();

      for (i_prob  = i_prob_begin; i_prob != i_prob_end; ++i_prob)
      {
        rVector[i] =  *i_prob;
        ++i;
      }

      if (includeCounts) {
        rVector[i] = backoff ? (**i_leaf).mDist.Counts() 
                             : (**i_leaf).mpBackoffDist->Counts() ;
        i++;
      }
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  // static
  void
  BDTree::
  ChCo_AdaptWeightVector(
      BasicVector<double>&        rVector, 
      const Matrix<double>&       rXformMatrix, 
      const BasicVector<double>&  rLM,
      const BasicVector<double>&  rNGrams, 
      int nIter, 
      double eps)
  {
    int                   i;
    BasicVector<double>   dV(rXformMatrix.Cols());
    BasicVector<double>   aux_vec(rVector.Length()); //rLM + V * rVector
    double                ngram_mass = rNGrams.Sum();


    // gradient descent iteration
    for (int i = 0; i < nIter; ++i) {
      // reset grad vector
      dV.Clear();
      aux_vec = rLM;

      // compute the aux vector
      aux_vec.AddCMtVMul(1.0, rXformMatrix, rVector);

      dV.AddCMtVMul(1.0, rXformMatrix, rNGrams);

      // update the vector
      rVector.AddCVMul(eps, dV);
    }
  }  // ChCo_AdaptWeightVector


  //***************************************************************************/
  //***************************************************************************/
  // BSetQuestion
  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion::
  BSetQuestion(std::istream& rStream, BDTreeHeader& rHeader)
  : BQuestion(), mSet(rHeader.mVocabSize), mPred(rHeader.mOrder)
  {
    Read(rStream, rHeader);
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BSetQuestion::
  Set(const NGram::TokenType& rToken)
  {
    mSet[rToken] = 1;
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BSetQuestion::
  Unset(const NGram::TokenType& rToken)
  {
    //BSetQuestion::SetType::const_iterator i_ngram;

    //i_ngram = mSet.find(rToken);
    //if (i_ngram != mSet.end()) {
    //  mSet.erase(i_ngram); 
    //}

    mSet[rToken] = 0;
  }

  //***************************************************************************/
  //***************************************************************************/
  bool
  BSetQuestion::
  Eval(const BQuestionTerm& rTerm) const
  {
    const NGram* p_ngram = reinterpret_cast<const NGram*>(&rTerm);
    //BSetQuestion::SetType::const_iterator i_ngram;

    //i_ngram = mSet.find((*p_ngram)[mPred]);
    //return !(i_ngram == mSet.end() || *i_ngram == 0);
    return mSet[(*p_ngram)[mPred]];
  }

  //***************************************************************************/
  //***************************************************************************/
  bool
  BSetQuestion::
  EvalRawToken(const NGram::TokenType& rToken) const
  {
    //const NGram* p_ngram = static_cast<const NGram*>(&rTerm);
    //BSetQuestion::SetType::const_iterator i_ngram;

    //i_ngram = mSet.find(rToken);
    //return !(i_ngram == mSet.end() || *i_ngram == 0);
    return mSet[rToken];
  }

  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion*
  BSetQuestion::
  Clone()
  {
    return new BSetQuestion(*this);
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BSetQuestion::
  Dump(std::ostream& rStream, const std::string& rPrefix) const
  {
    rStream << rPrefix << "================================" << std::endl;
    rStream << rPrefix << "= Simple question " << std::endl;
    rStream << rPrefix << "=  Predictor: " << mPred << std::endl;
    rStream << rPrefix << "=  Positive: ";

    BSetQuestion::SetType::const_iterator i;
    for (i = mSet.begin(); i!=mSet.end(); ++i) {
      rStream << *i<< "  ";
    }
    rStream << std::endl;
    rStream << rPrefix << "================================" << std::endl;
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BSetQuestion::
  DumpImplicit() const
  { Dump(std::cout, "Implicit dump "); }


} // namespace STK

// EOF
