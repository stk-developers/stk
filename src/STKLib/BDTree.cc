#include "BDTree.h"

#include <vector>
#include <iostream>

using std::vector;

namespace STK
{
  //***************************************************************************/
  //***************************************************************************/
  void
  Distribution::
  Dump(std::ostream& rStream, const std::string& rPrefix) const
  {
    rStream << rPrefix;

    Distribution::Container::const_iterator i;
    for (i = mVec.begin(); i!=mVec.end(); ++i) {
      rStream << i->first << ": " << i->second << "\t";
    }
    rStream << std::endl;
  }

  //***************************************************************************/
  //***************************************************************************/
  const Distribution::ProbType
  Distribution::
  Entropy() const
  {
    double e = 0;
    Distribution::Container::const_iterator i;

    for (i=mVec.begin(); i!=mVec.end(); ++i) {
      e -= i->second * log2(i->second);
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
  Distribution::
  ComputeFromCounts(const Data& rData)
  {
    Data::const_iterator  i;
    Distribution::Container::iterator  j;
    double                n = 0;
      
    // set probs to 0
    for (j=mVec.begin(); j!=mVec.end(); ++j) {
      j->second = 0.0;
    }

    for (i=rData.begin(); i!=rData.end(); ++i) {
      n += *i;
    }

    for (i=rData.begin(), j=mVec.begin(); j!=mVec.end(); ++i, ++j) {
      j->second /= n;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  Distribution::
  ComputeFromNGrams(const NGramSubset& rData)
  {
    mVec.clear();
    mN = 0;

    // collect the counts
    NGramSubset::NGramContainer::const_iterator i;

    for (i=rData.mData.begin(); i!=rData.mData.end(); ++i) {
      NGram::TokenType* p_token = &((**i)[0]);
      NGram::ProbType   counts = (*i)->Counts();

      if (mVec.find(*p_token) == mVec.end()) {
        mVec[*p_token] = counts;
      }
      else {
        mVec[*p_token] += counts;
      }

      mN += counts;
    }
      
    // normalize counts to probs
    Distribution::Container::iterator  j;
    
    for (j=mVec.begin(); j!=mVec.end(); ++j) {
      j->second /= mN;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  BDTree(NGramSubset& rData, BDTreeBuildTraits& rTraits,
      const std::string& rPrefix)
  {
    BuildFromNGrams(rData, rTraits, rPrefix);
  }


  //***************************************************************************/
  //***************************************************************************/
  BDTree::
  ~BDTree()
  {
    
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  BuildFromNGrams(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits,
      const std::string& rPrefix)
  {
    FLOAT   total_entropy;
    FLOAT   split_entropy;
    double  entropy_red;
    BQuestion* p_new_question = NULL;

    // compute distribution
    mDist.ComputeFromNGrams(rNGrams);
    if (rTraits.mVerbosity > 0) {
      std::cout << "BDTree::Info      " << rPrefix << "Data mass: " << mDist.Counts() << std::endl;
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

    // compute entropy
    total_entropy = rNGrams.Entropy();
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

    // if entropy reduction still big enough, build each subtree
    if (entropy_red >= rTraits.mMinReduction) {
      NGramSubset* p_split_data0 = new NGramSubset;
      NGramSubset* p_split_data1 = new NGramSubset;

      mpQuestion  = p_new_question;
      rNGrams.Split(*mpQuestion, *p_split_data0, *p_split_data1);

      mpTree0     = new BDTree(*p_split_data0, rTraits, rPrefix + "-  ");
      mpTree1     = new BDTree(*p_split_data1, rTraits, rPrefix + "+  ");

      delete p_split_data0;
      delete p_split_data1;
    } // else declare leaf
    else { 
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
  }


  //***************************************************************************/
  //***************************************************************************/
  BSetQuestion*
  BDTree::
  FindSimpleQuestion(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt) 
  {
    FLOAT          e = 0.0;
    double          best_e = -1.0;
    BSetQuestion*   p_best_question = NULL;


    // go through all predictors and find subset
    for (int i=1; i<=rTraits.mNPred; ++i) {
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
  FindSubset_Greedy(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt,
      int pred) 
  {
    bool    inserted  = true;
    bool    deleted   = true;
    double  e         = rNGrams.Entropy();
    BSetQuestion* p_new_question = new BSetQuestion(pred);

    std::map<NGram::TokenType, int> vocab;
    std::map<NGram::TokenType, int>::iterator i_vocab;
    std::map<NGram::TokenType, int>::iterator i_best_change;
    
    // collect vocabulary counts and total number of tokens
    NGramSubset::NGramContainer::const_iterator i;

    for (i=rNGrams.mData.begin(); i!=rNGrams.mData.end(); ++i) {
      NGram::TokenType* p_token = &((**i)[pred]);;
      vocab[*p_token] = 1;
    }

    // continue probing if any change was made in the previous iteration
    while (inserted || deleted) {
      inserted      = false;
      deleted       = false;

      // probe insertions
      for (i_vocab=vocab.begin(); i_vocab!=vocab.end(); ++i_vocab) {
        if (!p_new_question->EvalRawToken(i_vocab->first)) {
          p_new_question->Set(i_vocab->first);
          double this_e = rNGrams.SplitEntropy(*p_new_question);

          if (this_e < e) {
            e             = this_e;
            i_best_change = i_vocab;
            inserted      = true;
          }

          p_new_question->Unset(i_vocab->first);
        }
      }

      // commit insertion
      if (inserted) {
        p_new_question->Set(i_best_change->first);
      }

      // probe deletions
      for (i_vocab=vocab.begin(); i_vocab!=vocab.end(); ++i_vocab) {
        // only inspect positive questions
        if (p_new_question->EvalRawToken(i_vocab->first)) {
          p_new_question->Unset(i_vocab->first);
          double this_e = rNGrams.SplitEntropy(*p_new_question);

          if (this_e < e) {
            e             = this_e;
            i_best_change = i_vocab;
            deleted       = true;
          }

          p_new_question->Set(i_vocab->first);
        }
      }

      // commit deletion
      if (deleted) {
        p_new_question->Unset(i_best_change->first);
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
  void
  BDTree::
  Dump(std::ostream& rStream, const std::string& rPrefix) const
  {
    if (IsLeaf()) {
      mDist.Dump(rStream, rPrefix);
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
  // BSetQuestion
  //***************************************************************************/
  //***************************************************************************/
  void
  BSetQuestion::
  Set(const NGram::TokenType& rToken)
  {
    mSet.insert(rToken);
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BSetQuestion::
  Unset(const NGram::TokenType& rToken)
  {
    BSetQuestion::SetType::const_iterator i_ngram;

    i_ngram = mSet.find(rToken);
    if (i_ngram != mSet.end()) {
      mSet.erase(i_ngram); 
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  bool
  BSetQuestion::
  Eval(const BQuestionTerm& rTerm) const
  {
    const NGram* p_ngram = reinterpret_cast<const NGram*>(&rTerm);
    BSetQuestion::SetType::const_iterator i_ngram;

    i_ngram = mSet.find((*p_ngram)[mPred]);
    return !(i_ngram == mSet.end() || *i_ngram == 0);
  }

  //***************************************************************************/
  //***************************************************************************/
  bool
  BSetQuestion::
  EvalRawToken(const NGram::TokenType& rToken) const
  {
    //const NGram* p_ngram = static_cast<const NGram*>(&rTerm);
    BSetQuestion::SetType::const_iterator i_ngram;

    i_ngram = mSet.find(rToken);
    return !(i_ngram == mSet.end() || *i_ngram == 0);
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
    rStream << rPrefix << "=  Posivive: ";

    BSetQuestion::SetType::const_iterator i;
    for (i = mSet.begin(); i!=mSet.end(); ++i) {
      rStream << *i<< " ";
    }
    rStream << std::endl;
    rStream << rPrefix << "================================" << std::endl;
  }

  void
  BSetQuestion::
  DumpImplicit() const
  { Dump(std::cout, "Implicit dump "); }
} // namespace STK

// EOF
