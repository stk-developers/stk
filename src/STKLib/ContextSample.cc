#include "ContextSample.h"
#include "BasicVector.h"
#include "Matrix.h"
#include "stkstream.h"
#include "mymath.h"

#include <sstream>
#include <exception>
#include <utility>
#include <algorithm>
#include <string>


namespace STK {
  //***************************************************************************/
  //***************************************************************************/
  const std::string&
  VocabularyTable::
  IToA(const int index) const
  {
    if (index >= static_cast<int>(mIntMap.size())) {
      throw std::runtime_error(std::string("Index ") + " out of bound");
    }
    return mIntMap[index][0];
  } // IToA(..)


  //***************************************************************************/
  //***************************************************************************/
  int
  VocabularyTable::
  AToI(const std::string& rString)
  {
    std::map<std::string,int>::iterator it = mStrMap.find(rString);

    if (it == mStrMap.end()) {
      if (mStrictMode) {
        throw std::runtime_error(std::string("Symbol ") + rString + " not found");
      }
      else {
        int new_i = mStrMap.size();
        mStrMap[rString] = new_i;
        return new_i;
      }
    }

    return it->second;
  } // AToI(...)
    
  
  //***************************************************************************/
  //***************************************************************************/
  VocabularyTable&
  VocabularyTable::
  LoadFromFile(const std::string& rFName, bool isExtended) 
  {
    try {
      // open stream
      IStkStream stream(rFName.c_str());

      // testing...
      if (!stream.good()) {
        throw std::runtime_error(std::string("Cannot open vocabulary file ") + rFName);
      }

      if (isExtended) {
        LoadFromStreamExtended(stream);
      }
      else {
        LoadFromStream(stream);
      }

      stream.close();
      return *this;
    }
    catch (std::runtime_error& rError) {
      RethrowMessage("");
      throw;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  VocabularyTable&
  VocabularyTable::
  LoadFromStream(std::istream& rIStream)
  {
    try {
      int i=0;
      std::vector<std::string> tmp_vec;
    
      tmp_vec.push_back("");

      while (!rIStream.eof()) {
        rIStream >> tmp_vec[0] >> std::ws; 

        if (rIStream.bad() || rIStream.fail()) {
          throw std::runtime_error("Cannot read here...");
        }

        mStrMap[tmp_vec[0]] = i;
        mIntMap.push_back(tmp_vec);
        i++;
      }
      return *this;
    }
    catch (std::runtime_error& rError) {
      RethrowMessage("");
      throw;
    }
  }



  //***************************************************************************/
  //***************************************************************************/
  VocabularyTable&
  VocabularyTable::
  LoadFromStreamExtended(std::istream& rIStream)
  {
    try {
      int i=0;
      std::string tmp_line;
      std::string tmp_str;

      std::cout << "Extended... "  << std::endl;

      // skip file initial white spaces
      rIStream >>std::ws;

      // read line by line until the end
      while (!rIStream.eof()) {
        std::vector<std::string> tmp_collection;

        std::getline(rIStream, tmp_line);

        if (rIStream.bad() || rIStream.fail()) {
          throw std::runtime_error("Cannot read here...");
        }
        
        // skip trailing whitespaces so that eof detection works correctly
        rIStream >> std::ws;

        std::stringstream line_stream(tmp_line);
        line_stream >> std::ws;

        while (!line_stream.eof()) {
          line_stream >> tmp_str >> std::ws; 

          if (line_stream.bad() || line_stream.fail()) {
            throw std::runtime_error("Cannot read here...");
          }

          std::cout << tmp_str << "===" << i << std::endl;
          mStrMap[tmp_str] = i;
          tmp_collection.push_back(tmp_str);
        }

        mIntMap.push_back(tmp_collection);
        i++;
      }
      return *this;
    }
    catch (std::runtime_error& rError) {
      RethrowMessage("");
      throw;
    }
  }



  //***************************************************************************/
  //***************************************************************************/
  std::ostream&
  operator << (std::ostream& rOstr, const NGramPool& rWhat) 
  {
    std::vector<NGram*>::const_iterator i;

    for (i=rWhat.mData.begin(); i!=rWhat.mData.end(); ++i) {
      NGram& r_ngram = *(*i);

      for (int j=rWhat.Order()-1; j>=0; j--) {
        rOstr << r_ngram[j] << " ";
      }
      rOstr << r_ngram.Counts() << std::endl;
    }

    return rOstr;
  }


  //***************************************************************************/
  //***************************************************************************/
  // virtual
  void
  NGramPool::
  Clear() 
  {
    NGramContainer::iterator i;

    // go through the data and delete the token arrays inside the NGrams and 
    // then the NGrams themself
    for (i=mData.begin(); i!=mData.end(); ++i) {
      delete [] (*i)->mpTokens;
      delete *i;
    }

    mData.clear();
  }

  //***************************************************************************/
  //***************************************************************************/
  // virtual
  NGramPool::
  ~NGramPool()
  {
    Clear();
  } // ~NGramPool()


  //***************************************************************************/
  //***************************************************************************/
  void
  NGramPool::
  AddFromStream(std::istream& rStream, FLOAT weight, FLOAT threshold)
  {
    try 
    {
      std::string       line_buf;
      std::string       tmp_token;
      std::vector<std::string> token_chunk(Order());
      std::set<NGram*, NGramSubset::NGramCompare> ngram_map;
      std::set<NGram*, NGramSubset::NGramCompare>::iterator i_ngram;
      double            counts=0;
      int               token_counter = 0;
      NGram*            tmp_sample = NULL;

      // this could speed things up
      line_buf.reserve(LINE_BUF_DEFAULT_MIN_LENGTH);
      tmp_token.reserve(CONTEXT_SAMPLE_DEFAULT_MIN_LENGTH);

      // copy the original data to a map container, which is faster for 
      // insertion and update
      ngram_map.insert(mData.begin(), mData.end());

      rStream >> std::ws;

      // read till the end of file
      while (!rStream.eof()) 
      {
        // read one line and store in stringstream for parsing
        std::getline(rStream, line_buf);

        rStream >> std::ws;
        std::stringstream line_buf_stream(line_buf);

        if (rStream.bad() || rStream.fail()) {
          throw std::runtime_error("Cannot read here...");
        }

        // try to read the desired number of tokens and counts
        for (size_t i=0; i<Order(); ++i) {
          line_buf_stream >> token_chunk[i];
        }

        line_buf_stream >> counts;
        line_buf_stream >> std::ws;

        // if we failed reading the stream or we don't have eof yet
        // it means that wrong number of tokens are on line, so resume
        // the error state and continue reading
        if (line_buf_stream.fail() || !line_buf_stream.eof()) {
          line_buf_stream.clear();
          continue;
        }

        assert (counts >= 0);

        // create new NGram
        if (NULL == tmp_sample) {
          tmp_sample = this->pNewNGram();
        }
        
        // copy the chunk to the new NGram Token array
        std::vector<std::string>::reverse_iterator it; 
        size_t token_index = 0;

        for (it = token_chunk.rbegin(); it!=token_chunk.rend(); ++it, ++token_index) 
	{
	  int x;
          if(token_index == 0)
	    x = mpTargetTable->AToI(*it);
	  else
	    x = mpPredictorTable->AToI(*it);
	    
          (*tmp_sample)[token_index] = x;
        }


        if (counts > threshold) {
          tmp_sample->mCounts = counts * weight ;
          // if we have the ngram in the map, just add the counts
          if ((i_ngram = ngram_map.find(tmp_sample)) != ngram_map.end()) {
            (**i_ngram).mCounts += tmp_sample->mCounts;
          }
          else {
            ngram_map.insert(ngram_map.end(),tmp_sample);
            tmp_sample = NULL;
          }
          ++token_counter;
        }
        else {
          delete tmp_sample;
          tmp_sample = NULL;
        }
      } // while (!eof)

      // create new data vector out of the updated map container
      mData.clear();
      mData.insert(mData.begin(), ngram_map.begin(), ngram_map.end());

      // clear the temporary sample if last record was of different order
      if (NULL != tmp_sample) {
        delete tmp_sample;
      }
    }
    catch (std::runtime_error& rError) 
    {
      RethrowMessage("");
      throw;
    }
    catch (...) {
      RethrowMessage("");
      throw;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  NGramPool::
  AddFromFile(const std::string& rFileName, FLOAT weight, FLOAT threshold)
  {
    try {
      if (weight > 0) {
        std::cout << "Processing " << rFileName << std::endl;
        // open stream
        IStkStream stream(rFileName.c_str());

        // testing...
        if (!stream.good()) {
          throw std::runtime_error(std::string("Cannot open N-gram counts file ") +
              rFileName);
        }

        stream >> std::ws;
        // stream.seekg(0);
        AddFromStream(stream, weight, threshold);

        stream.close();
      }
      else {
        std::cout << "Skipping " << rFileName << " because weight is too low" 
          << std::endl;
      }
    }
    catch (std::runtime_error& rError) {
      RethrowMessage("");
      throw;
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  NGram*
  NGramPool::
  pNewNGram() 
  {
    NGram* tmp_sample     = new NGram();
    tmp_sample->mpTokens  = new NGram::TokenType[Order()];
    tmp_sample->mOrder    = Order();

    return tmp_sample;
  }


  //***************************************************************************/
  //***************************************************************************/
  NGram*
  NGramPool::
  FindNGram(const NGram* pNGram)
  {
    NGramContainer::iterator i;

    for (i=mData.begin(); i!=mData.end(); ++i) 
    {
      bool different=true;

      for (size_t j=0; j<Order() && different; ++j) {
        // TODO

      }
      if (!different) {
        return *i;
      }
    }

    return NULL;
  } // FindNGram(const NGram* pNGram)


  //****************************************************************************
  //****************************************************************************
  //****************************************************************************
  //****************************************************************************
  double
  NGramSubset::
  Mass() const 
  {
    double mass = 0.0;
    NGramContainer::const_iterator i;

    for (i=mData.begin(); i!=mData.end(); ++i) {
      mass += (*i)->Counts();
    }

    return mass;
  }


  //****************************************************************************
  //****************************************************************************
  bool
  NGramSubset::NGramCompare::
  operator()(NGram* s1, NGram* s2) const
  {
    int i;
    assert(s1->Order() == s2->Order());

    for (i = 0; i < s1->Order(); ++i) {
      if ((*s1)[i] < (*s2)[i]) {
        return true;
      }
      else if ((*s1)[i] > (*s2)[i]) {
        return false;
      }
    }
    return false;
  }


  //****************************************************************************
  //****************************************************************************
  void
  NGramSubset::
  Split(const BQuestion& rQuestion, NGramSubset& rData0, NGramSubset& rData1)
  {
    rData0.mData.reserve(mData.size());
    rData1.mData.reserve(mData.size());

    rData0.mOrder = mOrder;
    rData0.mpPool = mpPool;
    rData1.mOrder = mOrder;
    rData1.mpPool = mpPool;

    NGramContainer::iterator i;
    for (i = mData.begin(); i != mData.end(); ++i) {
      if (rQuestion.Eval(**i)) {
        rData1.mData.insert(rData1.mData.end(), *i);
      }
      else {
        rData0.mData.insert(rData0.mData.end(), *i);
      }
    }
  }    

  
  //****************************************************************************
  //****************************************************************************
  double
  NGramSubset::
  Entropy() const
  {
    std::map<NGram::TokenType, NGram::ProbType> vocab;
    std::map<NGram::TokenType, NGram::ProbType>::iterator i_vocab;
    double  N = 0;
    
    // collect vocabulary counts and total number of tokens
    NGramContainer::const_iterator i;
    for (i=mData.begin(); i!=mData.end(); ++i) {
      NGram::TokenType* p_token = &((**i)[0]);;
      NGram::ProbType   counts = (*i)->Counts();

      if (vocab.find(*p_token) == vocab.end()) {
        vocab[*p_token] = counts;
      }
      else {
        vocab[*p_token] += counts;
      }

      N += counts;
    }
    
    // accumulate some stats for entropy computation
    // e = -sum(#token/N * log(#token/N)) =
    //   = -1/N*sum(#token*log(#token)) + log(N)
    //   where sums are over all tokens
    double acc = 0;

    for (i_vocab=vocab.begin(); i_vocab!=vocab.end(); ++i_vocab) {
      acc += i_vocab->second * my_log(i_vocab->second);
    }

    return -acc/N + my_log(N);
  }


  //****************************************************************************
  //****************************************************************************
  double
  NGramSubset::
  SplitEntropy(const BQuestion& rQuestion) const
  {
    size_t          vocab_size = mpPool->pTargetTable()->Size();
    BasicVector<double> vocab0(vocab_size);
    BasicVector<double> vocab1(vocab_size);
    size_t          i;
    double          N0 = 0.0;
    double          N1 = 0.0;

    // collect vocabulary counts and total number of tokens
    NGramContainer::const_iterator it;
    for (it=mData.begin(); it!=mData.end(); ++it) {
      NGram::TokenType* p_token = &((**it)[0]);;
      NGram::ProbType   counts = (*it)->Counts();

      if (! rQuestion.Eval(**it)) {
        vocab0[*p_token] += counts;
        N0 += counts;
      }
      else {
        vocab1[*p_token] += counts;
        N1 += counts;
      }
    }
    
    // accumulate some stats for entropy computation
    // e = -sum(#token/N * log(#token/N)) =
    //   = -1/N*sum(#token*log(#token)) + log(N)
    //   where sums are over all tokens
    double acc0 = 0;
    double acc1 = 0;

    for (i=0; i<vocab_size; ++i) {
      if (vocab0[i] > 0) {
        acc0 += vocab0[i] * my_log(vocab0[i]);
      }
    }
    for (i=0; i<vocab_size; ++i) {
      if (vocab1[i] > 0) {
        acc1 += vocab1[i] * my_log(vocab1[i]);
      }
    }
    
    double p0 = N0/(N0+N1);
    double p1 = N1/(N0+N1);

    return p0*(-acc0/N0 + my_log(N0)) + p1*(-acc1/N1 + my_log(N1));
  }


  //****************************************************************************
  //****************************************************************************
  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  SplitMMI(const BQuestion& rQuestion, double alpha, double eta,
      double* pMass0, double* pMass1) const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    double          e0 = 0.0;
    double          e1 = 0.0;
    double          N0_all = 0.0;
    double          N1_all = 0.0;

    NGramSubsets::const_iterator        i_subset;

    Matrix<double> m_vocab0(this->size(), vocab_size);
    Matrix<double> m_vocab1(this->size(), vocab_size);

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i=0; i_subset!=this->end(); ++i_subset, i++) 
    {
      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) 
      {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        if (! rQuestion.Eval(**it)) {
          m_vocab0[i][*p_token]   += counts;
          N0_all                  += counts;
        }
        else {
          m_vocab1[i][*p_token]   += counts;
          N1_all                  += counts;
        }
      }
    }

    // TODO: smoothing here
    Matrix<double> m_probs0(m_vocab0);
    m_probs0.NormalizeRows().Log();

    Matrix<double> m_probs1(m_vocab1);
    m_probs1.NormalizeRows().Log();

    for (i = 0; i < this->size(); ++ i){
      BasicVector<double> aux_vec(this->size());

      double aux_sum0 = 0.0;
      double aux_sum1 = 0.0;
      double mul_const0 = 1.0;
      double mul_const1 = 1.0;
      
      for (size_t j = 0 ; j < m_vocab0.Cols(); ++j) {
        aux_sum0 += m_vocab0[i][j];
        aux_sum1 += m_vocab1[i][j];
      }

      mul_const0 = eta / (1+alpha* (aux_sum0 -1));
      mul_const1 = eta / (1+alpha* (aux_sum1 -1));
      
      aux_vec.AddCMVMul(mul_const0, m_probs0, m_vocab0[i]);
      e0 += aux_vec[i] - aux_vec.LogSumExp();

      aux_vec.Clear();

      aux_vec.AddCMVMul(mul_const1, m_probs1, m_vocab1[i]);
      e1 += aux_vec[i] - aux_vec.LogSumExp();
    }

    if (NULL != pMass0) {
      *pMass0 = N0_all;
    }
    if (NULL != pMass1) {
      *pMass1 = N1_all;
    }

    // return averaged MI's
    return (N0_all * e0 + N1_all * e1) / (N0_all + N1_all);
  }


  // //****************************************************************************
  // //****************************************************************************
  // double
  // NGramSubsets::
  // SplitMMI(const BQuestion& rQuestion, double eta) const
  // {
  //   size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
  //   size_t          i;
  //   double          e0 = 0.0;
  //   double          e1 = 0.0;
  //   double          N0_all = 0.0;
  //   double          N1_all = 0.0;

  //   NGramSubsets::const_iterator        i_subset;
  //   BasicVector<double>                 vocab_all0(vocab_size);
  //   BasicVector<double>                 vocab_all1(vocab_size);
  //   std::vector< BasicVector<double>* > vocabs0(this->size());
  //   std::vector< BasicVector<double>* > vocabs1(this->size());
  //   std::vector< BasicVector<double>* >::iterator i_vocabs0;
  //   std::vector< BasicVector<double>* >::iterator i_vocabs1;

  //   // initialize the vocabs
  //   for (i_vocabs0 = vocabs0.begin(); i_vocabs0 != vocabs0.end(); ++i_vocabs0) {
  //     *i_vocabs0 = new BasicVector<double>(vocab_size);
  //   }

  //   for (i_vocabs1 = vocabs1.begin(); i_vocabs1 != vocabs1.end(); ++i_vocabs1) {
  //     *i_vocabs1 = new BasicVector<double>(vocab_size);
  //   }

  //   // collect counts for each language and overal count for all langs
  //   for (i_subset=this->begin(), i_vocabs0=vocabs0.begin(), i_vocabs1=vocabs1.begin(); 
  //       i_subset!=this->end(); 
  //       ++i_subset, ++i_vocabs0, ++i_vocabs1) 
  //   {
  //     double          N0 = 0.0;
  //     double          N1 = 0.0;

  //     // collect vocabulary counts and total number of tokens
  //     NGramSubset::NGramContainer::const_iterator it;
  //     for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) 
  //     {
  //       NGram::TokenType* p_token = &((**it)[0]);;
  //       NGram::ProbType   counts = (*it)->Counts();

  //       if (! rQuestion.Eval(**it)) {
  //         vocab_all0[*p_token]    += counts;
  //         (**i_vocabs0)[*p_token] += counts;
  //         N0_all                  += counts;
  //         N0                      += counts;
  //       }
  //       else {
  //         (**i_vocabs1)[*p_token] += counts;
  //         vocab_all1[*p_token]    += counts;
  //         N1_all                  += counts;
  //         N1                      += counts;
  //       }
  //     }

  //     // normalize the stats
  //     // for (i=0; i<vocab_size; ++i) {
  //     //   (**i_vocabs0)[i] /= N0;
  //     //   (**i_vocabs1)[i] /= N1;
  //     // }
  //   }

  //   // normalize the overals
  //   for (i=0; i<vocab_size; ++i) {
  //     vocab_all0[i] /= N0_all;
  //     vocab_all1[i] /= N1_all;
  //   }

  //   // compute the MI's
  //   for (i_vocabs0=vocabs0.begin(), i_vocabs1=vocabs1.begin(); 
  //       i_vocabs0!=vocabs0.end(); 
  //       ++i_vocabs0, ++i_vocabs1) 
  //   {
  //     for (i = 0; i < vocab_size; ++i) {
  //       (**i_vocabs0)[i] /= N0_all;
  //       (**i_vocabs1)[i] /= N1_all;

  //       if ((**i_vocabs0)[i] > 0.0 && vocab_all0[i] > 0.0) {
  //         e0 += (**i_vocabs0)[i] * log((**i_vocabs0)[i]) 
  //             - (**i_vocabs0)[i] * log(vocab_all0[i]);
  //       }
  //       if ((**i_vocabs1)[i] > 0.0 && vocab_all1[i] > 0.0) {
  //         e1 += (**i_vocabs1)[i] * log((**i_vocabs1)[i]) 
  //             - (**i_vocabs1)[i] * log(vocab_all1[i]);
  //       }
  //     }
  //   }

  //   for (i_vocabs0 = vocabs0.begin(); i_vocabs0 != vocabs0.end(); ++i_vocabs0) {
  //     delete *i_vocabs0;
  //   }

  //   for (i_vocabs1 = vocabs1.begin(); i_vocabs1 != vocabs1.end(); ++i_vocabs1) {
  //     delete *i_vocabs1;
  //   }

  //   // return averaged MI's
  //   return (N0_all * e0 + N1_all * e1) / (N0_all + N1_all);
  // }


  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  MMI(double eta) const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    double          e = 0.0;

    NGramSubsets::const_iterator        i_subset;

    Matrix<double> m_vocab(this->size(), vocab_size);

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i=0;
         i_subset!=this->end(); 
         ++i_subset, ++i) 
    {
      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        m_vocab[i][*p_token] += counts;
      }
    }

    m_vocab.NormalizeRows();
    Matrix<double> m_probs(m_vocab);
    m_probs.Log();

    for (i = 0; i < this->size(); ++ i){
      BasicVector<double> aux_vec(this->size());

      aux_vec.AddCMVMul(eta, m_probs, m_vocab[i]);
      e += aux_vec[i] - aux_vec.LogSumExp();
    }

    return e;
  }


  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  ParallelEntropy() const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    size_t          j;
    double          n_all = 0.0;
    double          e = 0.0;

    NGramSubsets::const_iterator        i_subset;
    Matrix<double> m_vocab(this->size(), vocab_size);

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i = 0;
         i_subset!=this->end(); 
         ++i_subset,  ++i) 
    {
      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        m_vocab[i][*p_token]   += counts;
        n_all                  += counts;
      }
    }

    for (i = 0; i< m_vocab.Rows(); ++i) {
      for (j = 0; j < m_vocab.Cols(); ++j) {
        double p = m_vocab[i][j] / n_all;
        if (p > 0.0) {
          e += p * my_log(p);
        }
      }
    }

    m_vocab.Destroy();

    return -e;
  }


  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  ParallelSplitEntropy(const BQuestion& rQuestion, double* pMass0, double* pMass1) const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    size_t          j;
    double          e0 = 0.0;
    double          e1 = 0.0;
    double          N0_all = 0.0;
    double          N1_all = 0.0;

    NGramSubsets::const_iterator        i_subset;

    Matrix<double> m_vocab0(this->size(), vocab_size);
    Matrix<double> m_vocab1(this->size(), vocab_size);

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i = 0; i_subset!=this->end(); ++i_subset, ++i) 
    {
      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) 
      {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        assert(counts >= 0);

        if (! rQuestion.Eval(**it)) {
          m_vocab0[i][*p_token]   += counts;
          N0_all                  += counts;
        }
        else {
          m_vocab1[i][*p_token]   += counts;
          N1_all                  += counts;
        }
      }
    }

    for (i = 0; i < m_vocab0.Rows(); ++i) {
      for (j = 0; j < m_vocab0.Cols(); ++j) {
        double p0 = m_vocab0[i][j];
        double p1 = m_vocab1[i][j];
        
        if (p0 > DBL_EPSILON && N0_all > DBL_EPSILON) {
          p0 /= N0_all;
          e0 += p0 * my_log(p0);
        }
        if (p1 > DBL_EPSILON && N1_all > DBL_EPSILON) {
          p1 /= N1_all;
          e1 += p1 * my_log(p1);
        }
      }
    }

    if (NULL != pMass0) {
      *pMass0 = N0_all;
    }
    if (NULL != pMass1) {
      *pMass1 = N1_all;
    }

    return -(N0_all * e0 + N1_all * e1) / (N0_all + N1_all);
  }


  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  Entropy() const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    double          e = 0.0;
    double          N_all = 0.0;

    NGramSubsets::const_iterator        i_subset;
    BasicVector<double>                 vocab_all(vocab_size);

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(); i_subset!=this->end(); ++i_subset) 
    {
      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        vocab_all[*p_token]    += counts;
        N_all                  += counts;
      }
    }

    // normalize the overals
    for (i=0; i<vocab_size; ++i) {
      vocab_all[i] /= N_all;
    }

    // compute the entropy
    for (i = 0; i < vocab_size; ++i) {
      if (vocab_all[i] > 0.0) {
        e -= vocab_all[i] * my_log(vocab_all[i]);
      }
    }

    return e;
  }


  //****************************************************************************
  //****************************************************************************
  void
  NGramSubsets::
  Split(const BQuestion& rQuestion, NGramSubsets& rData0, NGramSubsets& rData1)
  {
    iterator it;
    for (it = this->begin(); it != this->end(); ++it) {
      iterator it0 = rData0.insert(rData0.end(),NGramSubset(it->Order()));
      iterator it1 = rData1.insert(rData1.end(),NGramSubset(it->Order()));
      it0->mpPool = it->mpPool;
      it1->mpPool = it->mpPool;
      it->Split(rQuestion, *it0, *it1);
    }
  }

  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  Mass() const
  {
    double mass =0.0;
    const_iterator it;
    for (it = this->begin(); it != this->end(); ++it) {
      mass += it->Mass();
    }
    return mass;
  }

  //****************************************************************************
  //****************************************************************************
  size_t
  NGramSubsets::
  TokenCount() const
  {
    size_t token_count = 0;
    const_iterator it;

    for (it = this->begin(); it != this->end(); ++it) {
      token_count += it->TokenCount();
    }
    return token_count;
  }

  //****************************************************************************
  //****************************************************************************
  void
  NGramSubsets::
  Destroy()
  {
    // iterator it;
    // for (it = begin(); it != end; ++it) {
    //   delete *it;
    // }
  }

  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  SplitEntropy(const BQuestion& rQuestion) const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    double          e0 = 0.0;
    double          e1 = 0.0;
    double          N0_all = 0.0;
    double          N1_all = 0.0;

    NGramSubsets::const_iterator        i_subset;
    BasicVector<double>                 vocab_all0(vocab_size);
    BasicVector<double>                 vocab_all1(vocab_size);

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(); 
        i_subset!=this->end(); 
        ++i_subset) 
    {
      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) 
      {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        if (! rQuestion.Eval(**it)) {
          vocab_all0[*p_token]    += counts;
          N0_all                  += counts;
        }
        else {
          vocab_all1[*p_token]    += counts;
          N1_all                  += counts;
        }
      }
    }

    // normalize the overals
    for (i=0; i<vocab_size; ++i) {
      vocab_all0[i] /= N0_all;
      vocab_all1[i] /= N1_all;
    }

    // compute the MMI
    for (i = 0; i < vocab_size; ++i) {
      if (vocab_all0[i] > 0.0) {
        e0 += vocab_all0[i] * my_log(vocab_all0[i]);
      }
      if (vocab_all1[i] > 0.0) {
        e1 += vocab_all1[i] * my_log(vocab_all1[i]);
      }
    }

    return - (N0_all * e0 + N1_all * e1) / (N0_all + N1_all);
  }


  //****************************************************************************
  //****************************************************************************
  // double
  // NGramSubset::
  // SplitEntropy(const BQuestion& rQuestion) const
  // {
  //   std::map<NGram::TokenType, NGram::ProbType> vocab0;
  //   std::map<NGram::TokenType, NGram::ProbType> vocab1;
  //   std::map<NGram::TokenType, NGram::ProbType>::iterator i_vocab;
  //   double  N0 = 0;
  //   double  N1 = 0;
  //   
  //   // collect vocabulary counts and total number of tokens
  //   NGramContainer::const_iterator i;
  //   for (i=mData.begin(); i!=mData.end(); ++i) {
  //     NGram::TokenType* p_token = &((**i)[0]);;
  //     NGram::ProbType   counts = (*i)->Counts();

  //     if (! rQuestion.Eval(**i)) {
  //       if (vocab0.find(*p_token) == vocab0.end()) {
  //         vocab0[*p_token] = counts;
  //       }
  //       else {
  //         vocab0[*p_token] += counts;
  //       }

  //       N0 += counts;
  //     }
  //     else {
  //       if (vocab1.find(*p_token) == vocab1.end()) {
  //         vocab1[*p_token] = counts;
  //       }
  //       else {
  //         vocab1[*p_token] += counts;
  //       }

  //       N1 += counts;
  //     }
  //   }
  //   
  //   // accumulate some stats for entropy computation
  //   // e = -sum(#token/N * log(#token/N)) =
  //   //   = -1/N*sum(#token*log(#token)) + log(N)
  //   //   where sums are over all tokens
  //   double acc0 = 0;
  //   double acc1 = 0;

  //   for (i_vocab=vocab0.begin(); i_vocab!=vocab0.end(); ++i_vocab) {
  //     acc0 += i_vocab->second * log(i_vocab->second);
  //   }
  //   for (i_vocab=vocab1.begin(); i_vocab!=vocab1.end(); ++i_vocab) {
  //     acc1 += i_vocab->second * log(i_vocab->second);
  //   }
  //   
  //   double p0 = N0/(N0+N1);
  //   double p1 = N1/(N0+N1);

  //   return p0*(-acc0/N0 + log(N0)) + p1*(-acc1/N1 + log(N1));
  // }

} // namespace STK

// EOF
