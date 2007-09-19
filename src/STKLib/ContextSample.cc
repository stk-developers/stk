#include "ContextSample.h"
#include "BasicVector.h"
#include "stkstream.h"

#include <sstream>
#include <exception>
#include <utility>

namespace STK {
  //***************************************************************************/
  //***************************************************************************/
  const std::string&
  VocabularyTable::
  IToA(const int index) const
  {
    if (index >= mIntMap.size()) {
      throw std::runtime_error("Index out of bound");
    }
    return mIntMap[index];
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
        throw std::runtime_error("Symbol not found");
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
  LoadFromFile(const std::string& rFName) 
  {
    try {
      // open stream
      IStkStream stream(rFName.c_str());

      // testing...
      if (!stream.good()) {
        std::runtime_error("Cannot open vocabulary file");
      }

      LoadFromStream(stream);

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
      std::string tmp_str;

      while (!rIStream.eof()) {
        rIStream >> tmp_str >> std::ws; 

        if (rIStream.bad() || rIStream.fail()) {
          throw std::runtime_error("Cannot read here...");
        }

        mStrMap[tmp_str] = i;
        mIntMap.push_back(tmp_str);
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
  AddFromStream(std::istream& rStream, const double& weight)
  {
    try 
    {
      std::string       line_buf;
      std::string       tmp_token;
      std::vector<std::string> token_chunk(Order());
      double            counts;
      int               token_counter = 0;

      // this could speed things up
      line_buf.reserve(LINE_BUF_DEFAULT_MIN_LENGTH);
      tmp_token.reserve(CONTEXT_SAMPLE_DEFAULT_MIN_LENGTH);

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
        for (int i=0; i<Order(); ++i) {
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

        // this mechanizm should provide some speedup, TODO: check this
        if (mData.size()+1 < mData.capacity()) {
          //mData.reserve(mData.capacity() * 2); 
        }
        
        // create new NGram
        NGram* tmp_sample = this->pNewNGram();
        tmp_sample->mCounts = counts * weight ;
        
        // copy the chunk to the new NGram Token array
        std::vector<std::string>::reverse_iterator it; 
        size_t token_index = 0;

        for (it = token_chunk.rbegin(); it!=token_chunk.rend(); ++it, ++token_index) {
          int x = mpPredictorTable->AToI(*it); 
          (*tmp_sample)[token_index] = x;
        }

        mData.push_back(tmp_sample);

        ++token_counter;
      } // while (!eof)
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
  AddFromFile(const std::string& rFileName, const double& weight)
  {
    try {
      std::cout << "Processing " << rFileName << std::endl;
      // open stream
      IStkStream stream(rFileName.c_str());

      // testing...
      if (!stream.good()) {
        std::runtime_error("Cannot open N-gram counts file");
      }

      stream.seekg(0);
      AddFromStream(stream, weight);

      stream.close();
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
    NGram* tmp_sample = new NGram();
    tmp_sample->mpTokens = new NGram::TokenType[Order()];

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

      for (int j=0; j<Order() && different; ++j) {
        // TODO

      }
      if (!different) {
        return *i;
      }
    }

    return NULL;
  } // FindNGram(const NGram* pNGram)


  const bool
  NGramPool::
  CompareNGram(const NGram* pFirst, const NGram* pSecond) const
  {
    for (int i=0; i<Order(); ++i) {
      // TODO
    }
  }

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
        rData1.mData.push_back(*i);
      }
      else {
        rData0.mData.push_back(*i);
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
      acc += i_vocab->second * log(i_vocab->second);
    }

    return -acc/N + log(N);
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
        acc0 += vocab0[i] * log(vocab0[i]);
      }
    }
    for (i=0; i<vocab_size; ++i) {
      if (vocab1[i] > 0) {
        acc1 += vocab1[i] * log(vocab1[i]);
      }
    }
    
    double p0 = N0/(N0+N1);
    double p1 = N1/(N0+N1);

    return p0*(-acc0/N0 + log(N0)) + p1*(-acc1/N1 + log(N1));
  }


  //****************************************************************************
  //****************************************************************************
  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  SplitMMI(const BQuestion& rQuestion) const
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
    std::vector< BasicVector<double>* > vocabs0(this->size());
    std::vector< BasicVector<double>* > vocabs1(this->size());
    std::vector< BasicVector<double>* >::iterator i_vocabs0;
    std::vector< BasicVector<double>* >::iterator i_vocabs1;

    // initialize the vocabs
    for (i_vocabs0 = vocabs0.begin(); i_vocabs0 != vocabs0.end(); ++i_vocabs0) {
      *i_vocabs0 = new BasicVector<double>(vocab_size);
    }

    for (i_vocabs1 = vocabs1.begin(); i_vocabs1 != vocabs1.end(); ++i_vocabs1) {
      *i_vocabs1 = new BasicVector<double>(vocab_size);
    }

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i_vocabs0=vocabs0.begin(), i_vocabs1=vocabs1.begin(); 
        i_subset!=this->end(); 
        ++i_subset, ++i_vocabs0, ++i_vocabs1) 
    {
      double          N0 = 0.0;
      double          N1 = 0.0;

      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) 
      {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        if (! rQuestion.Eval(**it)) {
          vocab_all0[*p_token]    += counts;
          (**i_vocabs0)[*p_token] += counts;
          N0_all                  += counts;
          N0                      += counts;
        }
        else {
          (**i_vocabs1)[*p_token] += counts;
          vocab_all1[*p_token]    += counts;
          N1_all                  += counts;
          N1                      += counts;
        }
      }

      // normalize the stats
      // for (i=0; i<vocab_size; ++i) {
      //   (**i_vocabs0)[i] /= N0;
      //   (**i_vocabs1)[i] /= N1;
      // }
    }

    // normalize the overals
    for (i=0; i<vocab_size; ++i) {
      vocab_all0[i] /= N0_all;
      vocab_all1[i] /= N1_all;
    }

    // compute the MI's
    for (i_vocabs0=vocabs0.begin(), i_vocabs1=vocabs1.begin(); 
        i_vocabs0!=vocabs0.end(); 
        ++i_vocabs0, ++i_vocabs1) 
    {
      for (i = 0; i < vocab_size; ++i) {
        (**i_vocabs0)[i] /= N0_all;
        (**i_vocabs1)[i] /= N1_all;

        if ((**i_vocabs0)[i] > 0.0 && vocab_all0[i] > 0.0) {
          e0 += (**i_vocabs0)[i] * log((**i_vocabs0)[i]) 
              - (**i_vocabs0)[i] * log(vocab_all0[i]);
        }
        if ((**i_vocabs1)[i] > 0.0 && vocab_all1[i] > 0.0) {
          e1 += (**i_vocabs1)[i] * log((**i_vocabs1)[i]) 
              - (**i_vocabs1)[i] * log(vocab_all1[i]);
        }
      }
    }

    for (i_vocabs0 = vocabs0.begin(); i_vocabs0 != vocabs0.end(); ++i_vocabs0) {
      delete *i_vocabs0;
    }

    for (i_vocabs1 = vocabs1.begin(); i_vocabs1 != vocabs1.end(); ++i_vocabs1) {
      delete *i_vocabs1;
    }

    // return averaged MI's
    return (N0_all * e0 + N1_all * e1) / (N0_all + N1_all);
  }


  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  MMI() const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    double          e = 0.0;
    double          N_all = 0.0;

    NGramSubsets::const_iterator        i_subset;
    BasicVector<double>                 vocab_all(vocab_size);
    std::vector< BasicVector<double>* > vocabs(this->size());
    std::vector< BasicVector<double>* >::iterator i_vocabs;

    // initialize the vocabs
    for (i_vocabs = vocabs.begin(); i_vocabs != vocabs.end(); ++i_vocabs) {
      *i_vocabs = new BasicVector<double>(vocab_size);
    }

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i_vocabs=vocabs.begin();
         i_subset!=this->end(); 
         ++i_subset, ++i_vocabs) 
    {
      double          N = 0.0;

      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        (**i_vocabs)[*p_token] += counts;
        vocab_all[*p_token]    += counts;
        N                      += counts;
        N_all                  += counts;
      }
      
      // normalize
      // for (i=0; i<vocab_size; ++i) {
      //   (**i_vocabs)[i] /= N;
      // }
    }

    // normalize the overals
    for (i=0; i<vocab_size; ++i) {
      vocab_all[i] /= N_all;
    }

    // compute the MMI
    for (i_vocabs=vocabs.begin(); i_vocabs!=vocabs.end(); ++i_vocabs) {
      for (i = 0; i < vocab_size; ++i) {
        (**i_vocabs)[i] /= N_all;
        if ((**i_vocabs)[i] > 0.0 && vocab_all[i] > 0.0) {
          e += (**i_vocabs)[i] * log((**i_vocabs)[i]) 
             - (**i_vocabs)[i] * log(vocab_all[i]);
        }
      }
    }

    // delete vocabs
    for (i_vocabs = vocabs.begin(); i_vocabs != vocabs.end(); ++i_vocabs) {
      delete *i_vocabs;
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
    double          n_all = 0.0;
    double          e = 0.0;

    NGramSubsets::const_iterator        i_subset;
    std::vector< BasicVector<double>* > vocabs(this->size());
    std::vector< BasicVector<double>* >::iterator i_vocabs;

    // initialize the vocabs
    for (i_vocabs = vocabs.begin(); i_vocabs != vocabs.end(); ++i_vocabs) {
      *i_vocabs = new BasicVector<double>(vocab_size);
    }

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i_vocabs=vocabs.begin();
         i_subset!=this->end(); 
         ++i_subset, ++i_vocabs) 
    {
      double          N = 0.0;

      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        (**i_vocabs)[*p_token] += counts;
        N                      += counts;
        n_all                  += counts;
      }
      
      // normalize
      // for (i=0; i<vocab_size; ++i) {
      //   (**i_vocabs)[i] /= N;
      // }
    }

    // compute the MMI
    for (i_vocabs=vocabs.begin(); i_vocabs!=vocabs.end(); ++i_vocabs) {
      for (i = 0; i < vocab_size; ++i) {
        // normalize
        (**i_vocabs)[i] /= n_all;

        // add to entropy
        if ((**i_vocabs)[i] > 0.0) {
          e += (**i_vocabs)[i] * log((**i_vocabs)[i]);
        }
      }
    }

    // delete vocabs
    for (i_vocabs = vocabs.begin(); i_vocabs != vocabs.end(); ++i_vocabs) {
      delete *i_vocabs;
    }

    return -e;
  }


  //****************************************************************************
  //****************************************************************************
  double
  NGramSubsets::
  ParallelSplitEntropy(const BQuestion& rQuestion) const
  {
    size_t          vocab_size = this->front().mpPool->pTargetTable()->Size();
    size_t          i;
    double          e0 = 0.0;
    double          e1 = 0.0;
    double          N0_all = 0.0;
    double          N1_all = 0.0;

    NGramSubsets::const_iterator        i_subset;
    std::vector< BasicVector<double>* > vocabs0(this->size());
    std::vector< BasicVector<double>* > vocabs1(this->size());
    std::vector< BasicVector<double>* >::iterator i_vocabs0;
    std::vector< BasicVector<double>* >::iterator i_vocabs1;

    // initialize the vocabs
    for (i_vocabs0 = vocabs0.begin(); i_vocabs0 != vocabs0.end(); ++i_vocabs0) {
      *i_vocabs0 = new BasicVector<double>(vocab_size);
    }

    for (i_vocabs1 = vocabs1.begin(); i_vocabs1 != vocabs1.end(); ++i_vocabs1) {
      *i_vocabs1 = new BasicVector<double>(vocab_size);
    }

    // collect counts for each language and overal count for all langs
    for (i_subset=this->begin(), i_vocabs0=vocabs0.begin(), i_vocabs1=vocabs1.begin(); 
        i_subset!=this->end(); 
        ++i_subset, ++i_vocabs0, ++i_vocabs1) 
    {
      double          N0 = 0.0;
      double          N1 = 0.0;

      // collect vocabulary counts and total number of tokens
      NGramSubset::NGramContainer::const_iterator it;
      for (it=i_subset->mData.begin(); it!=i_subset->mData.end(); ++it) 
      {
        NGram::TokenType* p_token = &((**it)[0]);;
        NGram::ProbType   counts = (*it)->Counts();

        if (! rQuestion.Eval(**it)) {
          (**i_vocabs0)[*p_token] += counts;
          N0_all                  += counts;
          N0                      += counts;
        }
        else {
          (**i_vocabs1)[*p_token] += counts;
          N1_all                  += counts;
          N1                      += counts;
        }
      }

      // normalize the stats
      // for (i=0; i<vocab_size; ++i) {
      //   (**i_vocabs0)[i] /= N0;
      //   (**i_vocabs1)[i] /= N1;
      // }
    }

    // compute the entropy
    for (i_vocabs0=vocabs0.begin(), i_vocabs1=vocabs1.begin(); 
        i_vocabs0!=vocabs0.end(); 
        ++i_vocabs0, ++i_vocabs1) 
    {
      for (i = 0; i < vocab_size; ++i) {
        if ((**i_vocabs0)[i] > 0.0) {
          (**i_vocabs0)[i] /= N0_all;
          e0 += (**i_vocabs0)[i] * log((**i_vocabs0)[i]);
        }
        if ((**i_vocabs1)[i] > 0.0) {
          (**i_vocabs1)[i] /= N1_all;
          e1 += (**i_vocabs1)[i] * log((**i_vocabs1)[i]);
        }
      }
    }

    for (i_vocabs0 = vocabs0.begin(); i_vocabs0 != vocabs0.end(); ++i_vocabs0) {
      delete *i_vocabs0;
    }

    for (i_vocabs1 = vocabs1.begin(); i_vocabs1 != vocabs1.end(); ++i_vocabs1) {
      delete *i_vocabs1;
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
        e -= vocab_all[i] * log(vocab_all[i]);
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
      rData0.push_back(NGramSubset(it->Order()));
      rData1.push_back(NGramSubset(it->Order()));
      rData0.back().mpPool = it->mpPool;
      rData1.back().mpPool = it->mpPool;
      it->Split(rQuestion, rData0.back(), rData1.back());
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
      mass+=it->Mass();
    }
    return mass;
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
        e0 += vocab_all0[i] * log(vocab_all0[i]);
      }
      if (vocab_all1[i] > 0.0) {
        e1 += vocab_all1[i] * log(vocab_all1[i]);
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
