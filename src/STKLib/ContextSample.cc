#include "ContextSample.h"
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
  NGramPool::
  ~NGramPool()
  {
    NGramContainer::iterator i;

    // go through the data and delete the token arrays inside the NGrams and 
    // then the NGrams themself
    for (i=mData.begin(); i!=mData.end(); ++i) {
      delete [] (*i)->mpTokens;
      delete *i;
    }
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
          mData.reserve(mData.capacity() * 2); 
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
  void
  NGramSubset::
  Split(const BQuestion& rQuestion, NGramSubset& rData0, NGramSubset& rData1)
  {
    rData0.mData.reserve(mData.size());
    rData1.mData.reserve(mData.size());

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
  FLOAT
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
    // e = -sum(#token/N * log2(#token/N)) =
    //   = -1/N*sum(#token*log2(#token)) + log2(N)
    //   where sums are over all tokens
    double acc = 0;

    for (i_vocab=vocab.begin(); i_vocab!=vocab.end(); ++i_vocab) {
      acc += i_vocab->second * log2(i_vocab->second);
    }

    return -acc/N + log2(N);
  }


  //****************************************************************************
  //****************************************************************************
  FLOAT
  NGramSubset::
  SplitEntropy(const BQuestion& rQuestion) const
  {
    std::map<NGram::TokenType, NGram::ProbType> vocab0;
    std::map<NGram::TokenType, NGram::ProbType> vocab1;
    std::map<NGram::TokenType, NGram::ProbType>::iterator i_vocab;
    double  N0 = 0;
    double  N1 = 0;
    
    // collect vocabulary counts and total number of tokens
    NGramContainer::const_iterator i;
    for (i=mData.begin(); i!=mData.end(); ++i) {
      NGram::TokenType* p_token = &((**i)[0]);;
      NGram::ProbType   counts = (*i)->Counts();

      if (! rQuestion.Eval(**i)) {
        if (vocab0.find(*p_token) == vocab0.end()) {
          vocab0[*p_token] = counts;
        }
        else {
          vocab0[*p_token] += counts;
        }

        N0 += counts;
      }
      else {
        if (vocab1.find(*p_token) == vocab1.end()) {
          vocab1[*p_token] = counts;
        }
        else {
          vocab1[*p_token] += counts;
        }

        N1 += counts;
      }
    }
    
    // accumulate some stats for entropy computation
    // e = -sum(#token/N * log(#token/N)) =
    //   = -1/N*sum(#token*log(#token)) + log(N)
    //   where sums are over all tokens
    double acc0 = 0;
    double acc1 = 0;

    for (i_vocab=vocab0.begin(); i_vocab!=vocab0.end(); ++i_vocab) {
      acc0 += i_vocab->second * log(i_vocab->second);
    }
    for (i_vocab=vocab1.begin(); i_vocab!=vocab1.end(); ++i_vocab) {
      acc1 += i_vocab->second * log(i_vocab->second);
    }
    
    double p0 = N0/(N0+N1);
    double p1 = N1/(N0+N1);

    return p0*(-acc0/N0 + log(N0)) + p1*(-acc1/N1 + log(N1));
  }

} // namespace STK

// EOF
