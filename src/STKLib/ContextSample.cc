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

      for (size_t j=0; j<rWhat.Order(); ++j) {
        rOstr << r_ngram[j] << " ";
      }
      rOstr << r_ngram.Counts() << std::endl;
    }
  }


  //***************************************************************************/
  //***************************************************************************/
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
        std::vector<std::string>::iterator it;
        size_t token_index = 0;

        for (it=token_chunk.begin(); it!=token_chunk.end(); ++it, ++token_index) {
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



} // namespace STK

// EOF
