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
#include "mymath.h"

#include <vector>
#include <iostream>
#include <ctime>


namespace STK
{

  // // virtual
  // void
  // Distribution::
  // Dump(std::ostream& rStream, const std::string& rPrefix) const
  // {
  // }
  
  // //***************************************************************************/
  // //***************************************************************************/
  // virtual
  void
  Distribution::
  DumpImplicit() const
  {
    Dump(std::cout, " ");
  }
  
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


  //***************************************************************************/
  //***************************************************************************/
  MapDistribution::ProbType
  MapDistribution::
  operator [] (const NGram::TokenType& rToken) const
  {
    assert (rToken < static_cast<int>(mVocabSize)); 

    Container::const_iterator i = mMap.find(rToken);
    if(i != mMap.end()) {
      return i->second;
    }
    else {
      return 0.0;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  MapDistribution::ProbType&
  MapDistribution::
  operator [] (const NGram::TokenType& rToken)
  {
    assert (rToken < static_cast<int>(mVocabSize)); 
    return mMap[rToken];
  }


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
  VecDistribution::
  VecDistribution(std::istream& rStream, BDTreeHeader& rHeader)
  {
    Read(rStream, rHeader);
  }


  //***************************************************************************/
  //***************************************************************************/
  VecDistribution*
  VecDistribution::
  Clone()
  {
    return new VecDistribution(*this);
  }


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
        e -= static_cast<double>(*i) * my_log(static_cast<double>(*i));
      }
    }
    
    return e;
  }


  //***************************************************************************/
  //***************************************************************************/
  const Distribution::ProbType
  Distribution::
  SplitEntropy(const Distribution& rD1, const Distribution& rD2)
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
  VecDistribution::ProbType 
  VecDistribution::
  operator [] (const NGram::TokenType& rToken) const
  {
    assert (rToken < static_cast<int>(mVec.size()));
    return mVec[rToken];
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Merge(const Distribution& rDistr)
  {
    const VecDistribution* p_distr = dynamic_cast<const VecDistribution*>(&rDistr);

    double norm = 0.0;
    VecDistribution::Container::iterator i_this = mVec.begin();
    VecDistribution::Container::const_iterator i_parent = p_distr->mVec.begin();

    assert(mVec.size() == p_distr->mVec.size());

    while (i_this != mVec.end()) {
      *i_this = mN * (*i_this) + p_distr->mN * (*i_parent);
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

    mN += p_distr->mN;
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
  Smooth(const Distribution& rDistr, FLOAT r)
  {
    const VecDistribution* p_distr = dynamic_cast<const VecDistribution*>(&rDistr);

    double norm   = 0; //DBL_EPSILON;
    double new_mn = 0; //DBL_EPSILON;
    VecDistribution::Container::iterator        i_this   = mVec.begin();
    VecDistribution::Container::const_iterator  i_parent = p_distr->mVec.begin();

    assert(mVec.size() == p_distr->mVec.size());

    while (i_this != mVec.end()) {
      // compute weight
      double b =  *i_this*mN / (*i_this*mN + r) ;

      *i_this = b * (*i_this) + (1-b) * (*i_parent);
      norm    += *i_this;
      new_mn  += b * (*i_this * mN) + (1-b) * (*i_parent * p_distr->mN);

      ++i_this;
      ++i_parent;
    } // while (i_this != mVec.end() && i_parent != p_distr->mVec.end()) {

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
  Interpolate(const Distribution& rDistr, FLOAT r)
  {
    double norm   = 0; //DBL_EPSILON;
    //double new_mn = 0; //DBL_EPSILON;
    size_t i = 0;

    assert(this->Size() == rDistr.Size());

    for (i = 0; i<mVec.size(); i++) {
      (*this)[i] = r * (*this)[i] + (1-r) * rDistr[i];
      norm    += (*this)[i];

      i++;
    } // while (i_this != mVec.end() && i_parent != rDistr.mVec.end()) {

    if (norm > 0) {
      // normalize
      for (i = 0; i<mVec.size(); i++) {
        (*this)[i] /= norm;
      }
    }
    //mN = new_mn;
  } //Smooth(const VecDistribution& rDistr, double r)
  

  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  MapAdapt(const Distribution& rDistr, FLOAT r)
  {
    const VecDistribution* p_distr = dynamic_cast<const VecDistribution*>(&rDistr);

    double norm   = 0; //DBL_EPSILON;
    //double new_mn = 0; //DBL_EPSILON;
    VecDistribution::Container::iterator        i_this   = mVec.begin();
    VecDistribution::Container::const_iterator  i_parent = p_distr->mVec.begin();

    assert(mVec.size() == p_distr->mVec.size());

    while (i_this != mVec.end()) {
      // compute weight
      double b =  mN / (mN + r) ;

      *i_this = b * (*i_this) + (1-b) * (*i_parent);
      norm    += *i_this;

      ++i_this;
      ++i_parent;
    } // while (i_this != mVec.end() && i_parent != p_distr->mVec.end()) {

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
  BDTree(NGramSubsets& rData, NGramSubsets* rHeldoutData, BDTreeBuildTraits& rTraits,
      BDTree* pParent, const std::string& rPrefix)
  : mpDist(NULL), mpBackoffDist(NULL)
  {
    BuildFromNGrams(rData, rHeldoutData, rTraits, pParent, rPrefix);
  }

  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  BDTree(const BDTree& rOrig)
  {
    if (!rOrig.IsLeaf()) {
      mpDist      = NULL;
      mpTree0     = rOrig.mpTree0->Clone();
      mpTree1     = rOrig.mpTree1->Clone();
      mpQuestion  = rOrig.mpQuestion->Clone();
    }
    else {
      //mpDist      = new VecDistribution(*(rOrig.mpDist));
      mpDist      = rOrig.mpDist->Clone();
      mpTree0     = NULL;
      mpTree1     = NULL;
      mpQuestion  = NULL;
    }
    
    if (NULL == rOrig.mpBackoffDist) {
      mpBackoffDist = NULL;
    }
    else {
      //mpBackoffDist = new VecDistribution(*(rOrig.mpBackoffDist));
      mpBackoffDist = rOrig.mpBackoffDist->Clone();
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  BDTree(std::istream& rStream, BDTreeHeader& rHeader)
  : mpDist(NULL), mpBackoffDist(NULL)
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
    if (NULL != mpDist) {
      delete mpDist;
    }
    if (NULL != mpBackoffDist) {
      delete mpBackoffDist;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  BuildFromNGrams(NGramSubsets& rNGrams, NGramSubsets* rHeldoutNGrams, BDTreeBuildTraits& rTraits,
      BDTree* pParent, const std::string& rPrefix)
  {
    FLOAT     total_entropy, heldout_total_entropy;
    FLOAT     split_entropy, heldout_split_entropy;
    FLOAT     entropy_red, heldout_entropy_red;

    FLOAT     total_mmi;
    FLOAT     split_mmi;
    FLOAT     mmi_inc;
    BSetQuestion* p_new_question = NULL;
    int       depth = rTraits.mCurDepth;
    size_t    vocab_size = rNGrams[0].Parent().pTargetTable()->Size();



    if (NULL != pParent && rTraits.mSmoothR > 0) {
      // compute distribution as we  want to smooth right away, for
      // which we need the parent distribution
      if (rTraits.mUseMapDistribution) {
        //mpDist = new MapDistribution(vocab_size);
      }
      else {
        mpDist = new VecDistribution(vocab_size);
      }
      mpDist->ComputeFromNGramSubsets(rNGrams);
      mpDist->Smooth(*(pParent->mpDist), rTraits.mSmoothR);
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

      if (NULL == mpDist) {
        mpDist = new VecDistribution(vocab_size);
        mpDist->ComputeFromNGramSubsets(rNGrams);
      }
      mpQuestion = NULL;
      mpTree0    = NULL;
      mpTree1    = NULL;
      return;
    }

    // check for data sufficiency
    if (rNGrams.Mass() < rTraits.mMinInData) {
      if (rTraits.mVerbosity > 0) {
        std::cout << "BDTree::Info      " << rPrefix << "Leaf due to minimum input data criterion" << std::endl;
      }

      if (NULL == mpDist) {
        mpDist = new VecDistribution(vocab_size);
        mpDist->ComputeFromNGramSubsets(rNGrams);
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

      if(total_entropy > 0)
      {
        // find best predictor and return split entropy
        p_new_question = FindSimpleQuestion(rNGrams, rTraits, &split_entropy);

        if(rHeldoutNGrams != NULL)
        {
	  double mass0, mass1;
          // compute heldout entropy before split
          heldout_total_entropy = rHeldoutNGrams->ParallelEntropy();
	  heldout_split_entropy = rHeldoutNGrams->ParallelSplitEntropy(*p_new_question, &mass0, &mass1);
	  heldout_entropy_red   = heldout_total_entropy - heldout_split_entropy;
        }

        entropy_red = total_entropy - split_entropy;
      }
      else
      {
        total_entropy = 0;
        entropy_red = 0; 
      }


      if (rTraits.mVerbosity > 0) {
        if (entropy_red == 0) {
          std::cout << "BDTree::Info      " << rPrefix << "No predictor found " << std:: endl;
          std::cout << "BDTree::Info      " << rPrefix << "Data Entropy:      " << total_entropy << std::endl;
        }
        else {
	  if(rTraits.mLVCSR)
	    std::cout << "BDTree::Info      " << rPrefix << "Predictor:      " << p_new_question->Predictor()<< std::endl;	    
	  else
	    p_new_question->Dump(std::cout, std::string("BDTree::Info      ")+rPrefix);

          std::cout << "BDTree::Info      " << rPrefix << "Data Entropy:      " << total_entropy << std::endl;
          std::cout << "BDTree::Info      " << rPrefix << "Split Entropy:     " << split_entropy << std::endl;
          std::cout << "BDTree::Info      " << rPrefix << "Entropy reduction: " << entropy_red << std::endl;
	  if(rHeldoutNGrams != NULL)
	    std::cout << "BDTree::Info      " << rPrefix << "Heldout entropy reduction: " << heldout_entropy_red << std::endl;
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

      if (NULL == mpDist) {
        mpDist = new VecDistribution(vocab_size);
        mpDist->ComputeFromNGramSubsets(rNGrams);
      }
      mpQuestion = NULL;
      mpTree0    = NULL;
      mpTree1    = NULL;
    }
    else if (entropy_red == 0
      || (rHeldoutNGrams != NULL && heldout_entropy_red <= 0) 
      || (rHeldoutNGrams == NULL && entropy_red < rTraits.mMinReduction)) 
    {
      if (rTraits.mVerbosity > 0) 
      {
	if(rHeldoutNGrams != NULL)
	  std::cout << "BDTree::Info      " << rPrefix << "Leaf due to heldout data entropy reduction criterion" << std::endl;
	else
	  std::cout << "BDTree::Info      " << rPrefix << "Leaf due to minimum entropy reduction criterion" << std::endl;
      }

      if (NULL != p_new_question) {
        delete p_new_question;
      }

      if (NULL == mpDist) {
        mpDist = new VecDistribution(vocab_size);
        mpDist->ComputeFromNGramSubsets(rNGrams);
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

      NGramSubsets* p_split_hld_data0 = NULL;
      NGramSubsets* p_split_hld_data1 = NULL;
      if(rHeldoutNGrams != NULL)
      {
        p_split_hld_data0 = new NGramSubsets;
        p_split_hld_data1 = new NGramSubsets;
        rHeldoutNGrams->Split(*mpQuestion, *p_split_hld_data0, *p_split_hld_data1);
      }


      rTraits.mCurDepth = depth + 1;
      mpTree0     = new BDTree(*p_split_data0, p_split_hld_data0, rTraits, this, rPrefix + "-  ");
      rTraits.mCurDepth = depth + 1;
      mpTree1     = new BDTree(*p_split_data1, p_split_hld_data1, rTraits, this, rPrefix + "+  ");

      delete p_split_data0;
      delete p_split_data1;

      delete p_split_hld_data0;
      delete p_split_hld_data1;

      if (NULL != mpDist) {
        delete mpDist;
        mpDist = NULL;
      }
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
      rInfo.mTotalData += mpDist->Counts();
      rInfo.mAverageEntropy += mpDist->Counts() * mpDist->Entropy();
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion*
  BDTree::
  FindSimpleQuestion(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt) 
  {
    FLOAT           e = 0.0;
    FLOAT           best_e = -1.0;
    BSetQuestion*   p_best_question = NULL;

    // choose predictor randomly for a random tree
    if(rTraits.mRandomizeTree && rTraits.mLVCSR)
    {
      int rand_pred;
      srand((unsigned)time(0));


      // In case --morphpred=70 or anything, the value means that we keep those best predictors performing better that percentage to the best one
      // and then randomly choose between them
      if(rTraits.mMorphologicalPredictors > 1)
      {
        double* e_array = new double[rTraits.mOrder];
        int i, counter = 0;
        // go through all predictors and find subset
        for (i=1; i<rTraits.mOrder; ++i) 
        {
          BSetQuestion* p_tmp_question = FindSubset_Greedy(rNGrams, rTraits, &e, i);
          e_array[i-1] = e;
          delete p_tmp_question;
        }
        double min_e, max_e, limit;
        min_e = e_array[0];
        max_e = e_array[0];
        for(i = 0; i < rTraits.mOrder-1; i++)
	{
          if(e_array[i] > max_e)
            max_e = e_array[i];
          else if(e_array[i] < min_e)
            min_e = e_array[i];
        }
        // We're interested in lower entropy, that's why we do 100 - rTraits.mMorphologicalPredictors and then take everything that's lower
        limit = (max_e - min_e) * (100 - rTraits.mMorphologicalPredictors) / 100 + min_e;
        int preds_below_limit = 0;
        for(i = 0; i < rTraits.mOrder-1; i++)
          if(e_array[i] <= limit)
            preds_below_limit++;

        // Now random value is in [1, preds_above_limit+1)
        rand_pred = (int) ( ( (double)rand() / (  (double)(RAND_MAX) + (double)(1) ) ) * preds_below_limit ) + 1;
        for(i = 0; i < rTraits.mOrder-1; i++)
        {
          if(e_array[i] <= limit)
            counter++;
          if(counter == rand_pred)
	  {
            p_best_question = FindSubset_Greedy(rNGrams, rTraits, &e, i+1);
            best_e = e;
            break;
          }
	}
        delete e_array;
      }
      // In case of morph trees with --morphpred=1 , we do always keep word-questions in random shuffle and choose one additional predictor randomly
      else if(rTraits.mMorphologicalPredictors)
      {
        rand_pred = (int) ( ( (double)rand() / (  (double)(RAND_MAX) + (double)(1) ) ) * (rTraits.mOrder-1) ) + 1;

        std::string factor_type, wrd;
        NGram::TokenType ngram_token;
        NGramSubset::NGramContainer::const_iterator it;

        // go through all predictors and find subset
        for (int i=1; i<rTraits.mOrder; ++i) 
        {
          it = rNGrams[0].mData.begin();
          ngram_token = (**it)[i];
          factor_type = rNGrams[0].Parent().pPredictorTable()->IToA(ngram_token);
          if(i == rand_pred || factor_type[0] == 'W')
	  {
            BSetQuestion* p_tmp_question = FindSubset_Greedy(rNGrams, rTraits, &e, i);

            if (e < best_e || best_e < 0) {
              if (NULL != p_best_question) {
                delete p_best_question;
              }
              best_e = e;
              p_best_question = p_tmp_question;
            }
	  } 
        }
      }
      else
      {
        rand_pred = (int) ( ( (double)rand() / (  (double)(RAND_MAX) + (double)(1) ) ) * (rTraits.mOrder-1) ) + 1;
        BSetQuestion* p_tmp_question = FindSubset_Greedy(rNGrams, rTraits, &e, rand_pred);

        if (e < best_e || best_e < 0) {
          if (NULL != p_best_question) {
            delete p_best_question;
          }
          best_e = e;
          p_best_question = p_tmp_question;
        } 
      }
    }
    else if(rTraits.mRandomizePredictors) {
      int rand_pred = (int) (((double)rand() / ((double)(RAND_MAX) + 
              (double)(1.0))) * (rTraits.mOrder - 1)) + 1;

      p_best_question = FindSubset_Greedy(rNGrams, rTraits, &best_e, rand_pred);
    }
    else
    {
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
  FindSimpleQuestion_MMI(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitMMI) 
  {
    FLOAT           mmi = -10000;
    FLOAT           best_mmi = mmi - 1.0;
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
  FindSubset_Greedy(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt,
      int pred) 
  {
    bool    inserted  = true;
    bool    deleted   = true;
    double  e         = rNGrams.ParallelEntropy();
    double  this_e, this_log_likelihood; // log_likelihood
    double  mass0, mass1;

    //int     vocab_size = rNGrams[0].Parent().pTargetTable()->Size();
    int     target_vocab_size = rNGrams[0].Parent().pTargetTable()->Size();
    int     vocab_size = rNGrams[0].Parent().pPredictorTable()->Size();

    int     i_vocab, vocab_start, vocab_end;
    int     i_best_change;
    int     node_vocab_size;
    BSetQuestion* p_new_question = new BSetQuestion(pred, vocab_size);
    
    // It's not only this - we also implement Exchange algorithm there for finding the best question (as described in Peng Xu thesis)
    // It should be pretty fast. Implemented only for trigrams
    if(rTraits.mLVCSR && rTraits.mOrder >= 3)
    {
      NGramSubsets::const_iterator        i_subset;

      // NB!!! if vocab_size and node_vocab_size are different (we reduce the matrices somehow), some things must be changed in the below code
      
      // Wwe detect wich kinf of predictor we deal with and find corresponding vocabulary subset
      if(rTraits.mMorphologicalPredictors)
      {
	MorphVocabSection VocabSection(rNGrams, pred);
	node_vocab_size =  VocabSection.Size();
	vocab_start = VocabSection.FirstElem();
	vocab_end = VocabSection.LastElem();
      }
      else
      {
	node_vocab_size = vocab_size;
	vocab_start = 0;
	vocab_end = vocab_size-1;
      }

      srand((unsigned)time(0));

      // make initial rundom shuffle of the question set if it's a random tree
      if(rTraits.mRandomizeTree)
        for(i_vocab = vocab_start; i_vocab <= vocab_end; ++i_vocab)
          p_new_question->RandomShuffle(i_vocab);

      // construct helper "sparce" matrix
      SparseBigramMatrix BigramMatrix0(target_vocab_size, vocab_size), BigramMatrix1(target_vocab_size, vocab_size);

      // first iterate over ngrams and count, how many units do we have in each row (corresponding to predictor) and allocate memory

      // collect counts for each language and overal count for all langs
      for (i_subset=rNGrams.begin(); i_subset!=rNGrams.end(); ++i_subset) 
      {
        //std::cout << i_subset->Mass() << std::endl;
	BigramMatrix0.CreateSizeVector(*i_subset, *p_new_question, pred, false);
	BigramMatrix1.CreateSizeVector(*i_subset, *p_new_question, pred, true);

	BigramMatrix0.AllocateMem();
	BigramMatrix1.AllocateMem();

	BigramMatrix0.Fill(*i_subset, *p_new_question, pred, false);
	BigramMatrix1.Fill(*i_subset, *p_new_question, pred, true);

        // Ondra asks: I don't understand this line, log_likelihood is not used any longer
        // Why do we have to compute this???
	//log_likelihood = BigramMatrix0.CountLogLikelihood(BigramMatrix1);
      }

      int i_predictor;
      // If we randomly put a 1 for some predictior, but there is no data for it at all - switch it to 0
      for(i_predictor = vocab_start; i_predictor <= vocab_end; ++i_predictor)
      {
	if(p_new_question->EvalRawToken(i_predictor)
            && !BigramMatrix0.GetSizeVectorCell(i_predictor) && !BigramMatrix1.GetSizeVectorCell(i_predictor))
          p_new_question->Unset(i_predictor);
      }


      // Move elements from one set to another to find the best question (Exchange algorithm)
      //double this_e = rNGrams.ParallelSplitEntropy(*p_new_question, &mass0, &mass1);
      
      while (inserted || deleted) 
      {
	inserted      = false;
	deleted       = false;

	// move from A_ to A (insertion) - we try to insert as many as possible
	for(i_predictor = vocab_start; i_predictor <= vocab_end; ++i_predictor) // NB!!! if vocab_size and node_vocab_size are different (we reduce the matrices somehow), then there things must be changed 
	{
	  if (!p_new_question->EvalRawToken(i_predictor) && BigramMatrix0.GetSizeVectorCell(i_predictor)) 
	  {
	    p_new_question->Set(i_predictor);

	    this_log_likelihood = BigramMatrix1.InsertBasicElement(BigramMatrix0, i_predictor);
	    
	    // now we get log-likelihood. We can either maximize this or minimize entropy. The switch to entropy is:
	    mass0 = BigramMatrix0.GetTotalCountSum();
	    mass1 = BigramMatrix1.GetTotalCountSum();
	    this_e = -this_log_likelihood  / (mass0 + mass1);

	    // Let the question be inserted and save changes - i.e. matrices are modified after each insertion
	    if (this_e > 0 && this_e < e && mass0 >= rTraits.mMinInData && mass1 >= rTraits.mMinInData) 
	    {
	      e = this_e;
	      inserted = true;
	    }
	    // If reduction is insignificant, cancel the insertion - i.e. modify back the matrices
	    else
	    {
	      p_new_question->Unset(i_predictor);
	      BigramMatrix0.InsertBasicElement(BigramMatrix1, i_predictor);
	    }
	  }
	}

	// move from A to A_ (deletion) - we try to delete as many as possible
	for(i_predictor = vocab_start; i_predictor < vocab_end; ++i_predictor) // NB!!! if vocab_size and node_vocab_size are different (we reduce the matrices somehow), then there things must be changed 
	{
	  if (p_new_question->EvalRawToken(i_predictor) && BigramMatrix1.GetSizeVectorCell(i_predictor)) 
	  {
	    p_new_question->Unset(i_predictor);

	    this_log_likelihood = BigramMatrix0.InsertBasicElement(BigramMatrix1, i_predictor);
	    
	    // now we get log-likelihood. We can either maximize this or minimize entropy. The switch to entropy is:
	    mass0 = BigramMatrix0.GetTotalCountSum();
	    mass1 = BigramMatrix1.GetTotalCountSum();
	    this_e = -this_log_likelihood  / (mass0 + mass1);

	    // Let the question be inserted and save changes - i.e. matrices are modified after each insertion
	    if (this_e > 0 && this_e < e && mass0 >= rTraits.mMinInData && mass1 >= rTraits.mMinInData) 
	    {
	      e = this_e;
	      deleted = true;
	    }
	    // If reduction is insignificant, cancel the insertion - i.e. modify back the matrices
	    else
	    {
	      p_new_question->Set(i_predictor);
	      BigramMatrix1.InsertBasicElement(BigramMatrix0, i_predictor);
	    }
	  }
	}
      }
    }
    else
    {
      if (rTraits.mRandomizeQuestions) {
        BSetQuestion  aux_question(pred, vocab_size);
        double        aux_mass_set;
        double        aux_mass_unset = rNGrams.Mass();

        for(i_vocab = 0; i_vocab <= vocab_size; ++i_vocab) {
          // probe the question. Does it do any data split???? ...
          aux_question.Set(i_vocab);
          rNGrams.ParallelSplitEntropy(aux_question, &aux_mass_set, &mass0);
          aux_question.Unset(i_vocab);
          
          // ... if yes, then randomize it
          if (aux_mass_set != aux_mass_unset) {
            p_new_question->RandomShuffle(i_vocab);
          }
        }
      }

      // continue probing if any change was made in the previous iteration
      while (inserted || deleted) {
        inserted      = false;
        deleted       = false;

        // probe insertions
        for (i_vocab=0; i_vocab<vocab_size; ++i_vocab) {
          if (!p_new_question->EvalRawToken(i_vocab)) {
            p_new_question->Set(i_vocab);

            double this_e = rNGrams.ParallelSplitEntropy(*p_new_question, &mass0, &mass1);

            if (this_e > 0 && this_e < e && mass0 >= rTraits.mMinInData && 
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


            if (this_e > 0 && this_e < e && mass0 >= rTraits.mMinInData && 
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
      FLOAT* pSplitMMI, int pred) 
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
  FLOAT
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
      return (*(p_node->mpDist))[rNGram[0]];
    }
  } // ScoreNGram(const NGram& rNGram)


  //***************************************************************************/
  //***************************************************************************/
  FLOAT
  BDTree::
  ScoreNGramSubset(const NGramSubset& rNGrams)
  {
    NGramSubset::NGramContainer::const_iterator i;
    double  score = 0.0;

    for (i = rNGrams.mData.begin(); i != rNGrams.mData.end(); ++i) {
      double n_gram_score = ScoreNGram(**i);

      // we don't score the n-grams, whose prob is 0
      if (n_gram_score > 0) {
        score += (*i)->Counts() * my_log(n_gram_score);
      }
    }

    return score;
  } // ScoreNGramSubset(const NGramSubset& rNGrams, double r)


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  OutputNGramSubsetARPA(const NGramSubset& rNGrams, std::ostream& lm_stream, int rOrder)
  {
    const double log_e_10 = 2.30258;
    NGramSubset::NGramContainer::const_iterator i;
    const VocabularyTable *predictor_voc_table, *target_voc_table;
    target_voc_table = rNGrams.Parent().pTargetTable();
    predictor_voc_table = rNGrams.Parent().pPredictorTable();

    double  log_score;
    int j;
    NGram::TokenType voc_index;
    std::string word_string;

    for (i = rNGrams.mData.begin(); i != rNGrams.mData.end(); ++i) {
      double n_gram_score = ScoreNGram(**i);

      if(n_gram_score == 0)
	log_score = -99;
      else
      {
	log_score = my_log(n_gram_score);
	log_score /= log_e_10;
      }

      lm_stream << log_score << " ";
      for(j = 1; j <= rOrder; j++)
      {
	voc_index = (**i)[rOrder - j];
	if((rOrder - j) == 0)
	  word_string = target_voc_table->IToA(voc_index);
	else
	  word_string = predictor_voc_table->IToA(voc_index);

	lm_stream << word_string << " ";
      }

      lm_stream << std::endl;
    }
  }


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
        mpDist->Dump(rStream, rPrefix);
      }

    }
    else {
      std::string new_pref;
      mpQuestion->Dump(rStream, rPrefix);

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
  DumpImplicit() const
  {
    Dump(std::cout, " ");
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
      if (NULL == mpDist) {
        //mpDist = new VecDistribution(mpTree0->mpDist->Size());
        mpDist = mpTree0->mpDist->Clone();
      }
      mpDist->Reset();
      mpDist->Merge(*(mpTree0->mpDist));
      mpDist->Merge(*(mpTree1->mpDist));
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  ComputeBackoffDists(BDTree* pParent, FLOAT r)
  {
    mpBackoffDist = mpDist->Clone();

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

      VecDistribution* p_aux_distr = 
        dynamic_cast<VecDistribution*>(p_node->mpDist);

      if (NULL == p_aux_distr) {
        throw std::runtime_error("Cannot collect couns for non-vector distribution implementations");
      }

      p_aux_distr->mVec[r_ngram[0]]  += r_ngram.Counts();
      p_aux_distr->mN                += r_ngram.Counts();
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

    // fill the leaves collection
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
        p_node->mpBackoffDist = p_node->mpDist->Clone();
        p_node->mpDist->mN      = 0.0;
        p_node->mpDist->Reset();
      }

      (*(p_node->mpDist))[r_ngram[0]]  += r_ngram.Counts();
      p_node->mpDist->mN                += r_ngram.Counts();
    }

    if (rTraits.mMapAdapt) {
      // recompute all distributions
      for (i_leaf = leaf_map.begin(); i_leaf != leaf_map.end(); ++i_leaf) {
        (*i_leaf)->mpDist->Fix();
        if (r > 0) {
          (*i_leaf)->mpDist->MapAdapt(*((**i_leaf).mpBackoffDist), r);
        }
      }
    }
    else {
      // recompute all distributions
      for (i_leaf = leaf_map.begin(); i_leaf != leaf_map.end(); ++i_leaf) {
        (*i_leaf)->mpDist->Fix();
        if (r > 0) {
          (*i_leaf)->mpDist->Smooth(*((**i_leaf).mpBackoffDist), r);
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
      mpDist->Reset();
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
  FillLeafSupervector(BasicVector<FLOAT>& rVector, bool backoff,
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

    VecDistribution* p_aux_distr = 
      dynamic_cast<VecDistribution*>(leaves.front()->mpDist);

    if (NULL == p_aux_distr) {
      throw std::runtime_error("Cannot collect counts for non-vector distribution implementations");
    }

    vocab_size  = p_aux_distr->mVec.size();
    vec_size    = includeCounts ? n_leaves * vocab_size + n_leaves 
                                : n_leaves * vocab_size;

    // prepare new vector
    rVector.Destroy();
    rVector.Init(vec_size);
    rVector.Clear();

    // stack the leaves together and copy to the new supervector
    for (i_leaf = leaves.begin(); i_leaf != leaves.end(); ++i_leaf) {
      if (backoff) {
        p_aux_distr = dynamic_cast<VecDistribution*>((**i_leaf).mpBackoffDist);
      }
      else {
        p_aux_distr = dynamic_cast<VecDistribution*>((**i_leaf).mpDist);
      }


      if (NULL == p_aux_distr) {
        throw std::runtime_error("Cannot collect couns for non-vector distribution implementations");
      }

      VecDistribution::Container::iterator i_prob_begin = p_aux_distr->mVec.begin();
      VecDistribution::Container::iterator i_prob_end   = p_aux_distr->mVec.end();

      for (i_prob  = i_prob_begin; i_prob != i_prob_end; ++i_prob)
      {
        rVector[i] =  *i_prob;
        ++i;
      }

      if (includeCounts) {
        rVector[i] = backoff ? (**i_leaf).mpDist->Counts() 
                             : (**i_leaf).mpBackoffDist->Counts() ;
        i++;
      }
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  PushLeafSupervector(const std::vector<FLOAT>& rVector, bool backoff,
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

    VecDistribution* p_aux_distr = 
      dynamic_cast<VecDistribution*>(leaves.front()->mpDist);

    if (NULL == p_aux_distr) {
      throw std::runtime_error("Cannot collect couns for non-vector distribution implementations");
    }

    vocab_size  = p_aux_distr->mVec.size();
    vec_size    = includeCounts ? n_leaves * vocab_size + n_leaves 
                                : n_leaves * vocab_size;

    // stack the leaves together and copy to the new supervector
    for (i_leaf = leaves.begin(); i_leaf != leaves.end(); ++i_leaf) {
      if (backoff) {
        p_aux_distr = dynamic_cast<VecDistribution*>((**i_leaf).mpBackoffDist);
      }
      else {
        p_aux_distr = dynamic_cast<VecDistribution*>((**i_leaf).mpDist);
      }

      VecDistribution::Container::iterator i_prob_begin = p_aux_distr->mVec.begin();
      VecDistribution::Container::iterator i_prob_end   = p_aux_distr->mVec.end();

      for (i_prob  = i_prob_begin; i_prob != i_prob_end; ++i_prob)
      {
        *i_prob = rVector[i];
        ++i;
      }

      if (includeCounts) {
        if (backoff) {
          (**i_leaf).mpBackoffDist->mN = rVector[i];
        }
        else {
          (**i_leaf).mpDist->mN = rVector[i];
        }

        i++;
      }
    }
  }



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
  void
  BSetQuestion::
  RandomShuffle(const NGram::TokenType& rToken)
  {
    int rand_value;

    rand_value = rand();
    if(rand_value < RAND_MAX/2)
      mSet[rToken] = 0;
    else
      mSet[rToken] = 1;
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



  //***************************************************************************/
  //***************************************************************************/
  SparseBigramMatrix::~SparseBigramMatrix()
  {
    int i;
    for(i = 0; i < mPredictorVocabSize; i++)
      delete [] mPointerVector[i];
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  SparseBigramMatrix::CreateSizeVector(const NGramSubset& rNGrams, const BSetQuestion& rQuestion, int pred, bool YesAnswer)
  {
    int vocab_size = rNGrams.Parent().pPredictorTable()->Size();
    bool already_seen;
    BasicVector< NGram::TokenType> seen_contexts(vocab_size);
    NGram::TokenType prev_present_token = -1;

    // iterate over all N-grams at this particular node
    NGramSubset::NGramContainer::const_iterator it;
    for (it=rNGrams.mData.begin(); it!=rNGrams.mData.end(); ++it) 
    {
      already_seen = false;
      NGram::TokenType* predictor_token = &((**it)[pred]);
      NGram::TokenType* present_token = &((**it)[0]);
      NGram::ProbType   counts = (*it)->Counts();

      if (!(counts > 0)) {
        continue;
      }
      // assert(counts > 0); 

      // The N-grams are all sorted accordind to most current words (first 
      // present words, then -1 context etc.). We take benefit of this since 
      // we must count not different N-grams but N-1 grams with different 
      // predictors

      // we start reading N-grams with different present word
      if(*present_token != prev_present_token) {
	prev_present_token = *present_token;
	seen_contexts.Clear();
      }
      else if(seen_contexts[*predictor_token]) {
	already_seen = true;
      }


      if(!already_seen)
      {
	if((YesAnswer && rQuestion.Eval(**it))
	  || (!YesAnswer && !rQuestion.Eval(**it)))
	{ 
	  mSizeVector[*predictor_token]++;

	  seen_contexts[*predictor_token] = 1;
	}
      }
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  SparseBigramMatrix::AllocateMem()
  {
    int i, size;
    for(i = 0; i < mPredictorVocabSize; i++)  
    {
      size = mSizeVector[i];
      //this->mPointerVector[i] = 0;
      if(size)
	mPointerVector[i] = new SparseBigramCell[size];
      else 
	mPointerVector[i] = NULL;
    }

  }

  //***************************************************************************/
  //***************************************************************************/ 
  void
  SparseBigramMatrix::Fill(const NGramSubset& rNGrams, 
      const BSetQuestion& rQuestion, const int pred, const bool YesAnswer)
  {
    // iterate over all N-grams at this particular node
    NGramSubset::NGramContainer::const_iterator it;
    int i;
    bool bigram_found;

    for (it=rNGrams.mData.begin(); it!=rNGrams.mData.end(); ++it) 
    {
      NGram::TokenType* present_token = &((**it)[0]);
      NGram::TokenType* predictor_token = &((**it)[pred]);
      NGram::ProbType   counts = (*it)->Counts();

      if (!(counts > 0)) {
        continue;
      }
      //assert(counts > 0); 

      if((YesAnswer && rQuestion.Eval(**it))
	  || (!YesAnswer && !rQuestion.Eval(**it)))
      {
	bigram_found = false;
	for(i = 0; i < mSizeVector[*predictor_token]; i++)
	{
	  // do we start from zero or 1????
	  int     current_column = (this->mPointerVector[*predictor_token] + i)->GetColumn();
	  double  current_count  = (this->mPointerVector[*predictor_token] + i)->GetCount();

          // Ondra asks: WHY break ????
	  if (!(current_count > 0)) {
	    break;
          }

	  else if(*present_token == current_column)
	  {
	    (this->mPointerVector[*predictor_token] + i)->IncrementCount(counts);
	    mMarginalCountsPresent[*present_token] += counts;
	    mMarginalCountsPredictor[*predictor_token] += counts;
	    mTotalSum += counts;
	    bigram_found = true;

	    break;
	  }
	}

	if(!bigram_found)
	{
	  assert(i < mSizeVector[*predictor_token]);

	  (this->mPointerVector[*predictor_token] + i)->SetCell(*present_token, counts);
	  mMarginalCountsPresent[*present_token] += counts;
	  mMarginalCountsPredictor[*predictor_token] += counts;
	  mTotalSum += counts;
	}
      }
    }
  }

  //***************************************************************************/
  //***************************************************************************/  
  SparseBigramMatrix::ProbType
  SparseBigramMatrix::InsertBasicElement(SparseBigramMatrix& MoveFromThisSet, 
      const int Predictor)
  {
    int i, helper;
    double log_likelihood_sum = 0;
    double counts_moved = 0;
    SparseBigramCell *AuxCopy;


    // Iterate over all elements in the basic element (i.e. over all present 
    // words) to be inserted
    for(i = 0; i < MoveFromThisSet.mSizeVector[Predictor]; i++)
    {
      // do we start from zero or 1????
      int     current_column = (MoveFromThisSet.mPointerVector[Predictor] + i)->GetColumn();
      double  current_count  = (MoveFromThisSet.mPointerVector[Predictor] + i)->GetCount();

      // assert(current_count > 0);

      MoveFromThisSet.mMarginalCountsPresent[current_column] -= current_count;
      this->mMarginalCountsPresent[current_column] += current_count;
      counts_moved += current_count;
    }

    this->mTotalSum += counts_moved;
    MoveFromThisSet.mTotalSum -= counts_moved;
    // Actually we do not use predctor marginal count vector
    MoveFromThisSet.mMarginalCountsPredictor[Predictor] -= counts_moved;

    // assert(MoveFromThisSet.mMarginalCountsPredictor[Predictor] == 0);

    this->mMarginalCountsPredictor[Predictor] += counts_moved;

    // assert(this->mMarginalCountsPredictor[Predictor] == counts_moved);

    // update size vectors
    helper = this->mSizeVector[Predictor];
    this->mSizeVector[Predictor] =  MoveFromThisSet.mSizeVector[Predictor];
    MoveFromThisSet.mSizeVector[Predictor] = helper;
    // We updated marginal count matrices and TotalSum, now move whole basic elements
    AuxCopy = this->mPointerVector[Predictor];
    this->mPointerVector[Predictor] = MoveFromThisSet.mPointerVector[Predictor];
    MoveFromThisSet.mPointerVector[Predictor] = AuxCopy;
      
    log_likelihood_sum = this->CountLogLikelihood(MoveFromThisSet);

    return log_likelihood_sum;
  }

  //***************************************************************************/
  //***************************************************************************/  
  SparseBigramMatrix::ProbType
  SparseBigramMatrix::CountLogLikelihood(const SparseBigramMatrix& SecondSet) const
  {
    int i;
    double log_likelihood = 0, count_A, count_A_, total_A, total_A_;

    // Formula sum_w[C(w,A)logC(w,A) + C(w,A_complement)logC(w,A_complement)] - C(A)logC(A) - C(A_complement)logC(A_complement)      
    // Iterate over all elements in the basic element (i.e. over all present words) to be inserted
    for(i = 0; i < mPresentVocabSize; i++)
    {
      count_A  = this->mMarginalCountsPresent[i];
      count_A_ = SecondSet.mMarginalCountsPresent[i];

      if(count_A > 0)
        log_likelihood += count_A  * my_log(count_A);
      if(count_A_ > 0)
	log_likelihood += count_A_ * my_log(count_A_);
    }
    total_A  = this->mTotalSum;
    total_A_ = SecondSet.mTotalSum;

    // if (close_enough(total_A, 0.0, 10)) {
    //   total_A = 0.0;
    // }
    // if (close_enough(total_A_, 0.0, 10)) {
    //   total_A_ = 0.0;
    // }

    assert(total_A >= 0 && total_A_ >= 0);

    if(total_A > 0) {
      log_likelihood -= total_A  * my_log(total_A);
    }
    if(total_A_ > 0) {
      log_likelihood -= total_A_ * my_log(total_A_);
    }

    return log_likelihood;
  }

  //***************************************************************************/
  //***************************************************************************/  
  MorphVocabSection:: MorphVocabSection(NGramSubsets& rNGrams, int pred)
  {
    int i, vocab_total;
    std::string factor_type, wrd;
    NGram::TokenType ngram_token;
    NGramSubset::NGramContainer::const_iterator it;

    it = rNGrams[0].mData.begin();
    ngram_token = (**it)[pred];

    vocab_total = rNGrams[0].Parent().pPredictorTable()->Size();
    factor_type = rNGrams[0].Parent().pPredictorTable()->IToA(ngram_token);

    mFirst = 0;

    for(i =0 ; i < vocab_total; i++)
    {
      wrd = rNGrams[0].Parent().pPredictorTable()->IToA(i);
      if(factor_type[0] == wrd[0])
      {
	if(mFirst == 0)
	  mFirst = i;

	mLast  = i;
      }
      else
      {
	if(mFirst != 0)
	  break;
	else
	  continue;
      }
    }
  }


} // namespace STK

// EOF
