#include "common.h"
#include <algorithm>

#ifndef STK_MlfStream_tcc
#define STK_MlfStream_tcc

namespace STK
{
  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    BasicOMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    BasicOMlfStreamBuf(OStreamReference rOStream, size_t bufferSize)
    : mIsOpen(false), mOStream(rOStream)
    { }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    BasicOMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    ~BasicOMlfStreamBuf()
    {
      mOStream.flush();
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    int 
    BasicOMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    sync()
    {
      mOStream.flush();
      return 0;
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    typename _Traits::int_type 
    BasicOMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    overflow(typename _Traits::int_type c)
    {
      // we don't use buffer here... 
      if (mIsOpen) {
        if (_Traits::eof() == c) {
          return _Traits::not_eof(c);
        }
        // only pass the character to the stream
        mOStream.rdbuf()->sputc(c);

        // remember last char (in case we want to close)
        mLastChar = c;

        return c;
      }
      else {
        return _Traits::eof();
      }
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    void 
    BasicOMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    Close()
    {
      // if last character was not EOL, we need to insert it
      if (mLastChar != '\n') {
        mOStream.put('\n');
      }
      mOStream << ".\n";

      // flush the stream and declare the stream closed
      mOStream.flush();
      mIsOpen = false;
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    BasicOMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT> *
    BasicOMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    Open(const std::string& rFileName)
    {
      // retreive position
      std::streampos pos = mOStream.tellp();

      // write the initial "filename" in parantheses
      mOStream << '"' << rFileName << '"' << std::endl;
      mLastChar = '\n';

      // return NULL if we canot open
      if (!mOStream.good()) {
        return NULL;
      }

      // if ok, store the name position
      if (-1 != pos) {
        pos = mOStream.tellp();
        mLabels.Insert(rFileName, pos);
      }

      // set open flag and return this
      mIsOpen = true;
      return this;
    }


  //****************************************************************************
  //****************************************************************************
  // BasicIMlfStreamBuf section
  //
  //****************************************************************************
  //****************************************************************************

  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    BasicIMlfStreamBuf(IStreamReference rIStream, size_t bufferSize)
    : mIsOpen(false), mIsHashed(false), mIsEof(true), mState(IN_HEADER_STATE), 
      mIStream(rIStream), mLineBuffer()
    {
      // we reserve some place for the buffer...
      mLineBuffer.reserve(bufferSize);

      //StreamBufType::setg(mpBuffer, mpBuffer + bufferSize, mpBuffer + bufferSize);
      StreamBufType::setg(&(mLineBuffer.front()), &(mLineBuffer.back()), &(mLineBuffer.back()));
    }

  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    ~BasicIMlfStreamBuf()
    { 
    }

  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    void
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    Index()
    {
      // retreive position
      std::streampos orig_pos   = mIStream.tellg();
      int      orig_state = mState;

      // for streams like stdin, pos will by definition be -1, so we can only 
      // rely on sequential access and cannot hash it.
      if (-1 != orig_pos) {
        std::string aux_name;
        // we will constantly jump to next definition. the function automatically
        // hashes the stream if possible
        while (JumpToNextDefinition(aux_name)) 
        { }

        // move to the original position
        mIStream.clear();
        mIStream.seekg(orig_pos);
        mState = orig_state;
      }
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    bool
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    JumpToNextDefinition(std::string& rName)
    {
      if (!mIStream.good()) {
        return false;
      }

      // if we can, we will try to index the label
      std::streampos pos = mIStream.tellg();

      // we might be at a definition already, so first move one line further
      FillLineBuffer();

      // read lines till we get to definition again
      while (mIStream.good() && mState != IN_TITLE_STATE) {
        FillLineBuffer();
      }

      // decide what happened
      if (IN_TITLE_STATE == mState) {
        // if we can, we will try to index the label
        pos = mIStream.tellg();

        if (pos != static_cast<const std::streampos>(-1)) {
        // if (pos !=std::string::npos) {  // This line does not work under MinGW
          std::string line_buffer(mLineBuffer.begin(), mLineBuffer.end());
          STK::ParseHTKString(line_buffer, rName);
          mLabels.Insert(rName, pos);
        }

        return true;
      }
      else {
        // we have been hashing all the way through so we know that if this is 
        // is the EOF, we are done hashing this stream
        if (pos != static_cast<const std::streampos>(-1)) {
          mIsHashed = true;
        }

        // we are not in body state, so we just return false
        return false;
      }
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>*
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    Close()
    {
      if (!mIsOpen) {
        mIsEof = true;
        return NULL;
      }
      else {
        // if we try to close while in the body, we need to reach the end
        if (mState == IN_BODY_STATE) {
          while (mState == IN_BODY_STATE) {
            FillLineBuffer();
          }
        }

        // disable buffer mechanism
        StreamBufType::setg(&(mLineBuffer.front()), &(mLineBuffer.front()), 
            &(mLineBuffer.front()));

        mIsEof  = true;
        mIsOpen = false;

        return this;
      }
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>*
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    Open(const std::string& rFileName)
    {
      BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>* ret_val = NULL;

      // this behavior is compatible with ifstream
      if (mIsOpen) {
        Close();
        return NULL;
      }

      // retreive position
      std::streampos pos = mIStream.tellg();
      LabelRecord label_record;

      // for streams like stdin, pos will by definition be -1, so we can only 
      // rely on sequential access. At this place, we decide what to do
      if ((-1 != pos) && (mLabels.Find(rFileName, label_record))) {
        mIStream.seekg(label_record.mStreamPos);
        mState = IN_TITLE_STATE;

        // we don't want the other stream to be bad, so we transfer the 
        // flagbits to this stream
        if (!mIStream.good()) {
          mIStream.clear();
          mIsOpen = false;
          ret_val = NULL;
        }
        else {
          mIsOpen = true;
          mIsEof = false;
          ret_val = this;
        }
      }

      // we don't have sequential stream and we didn't find the label, but
      // we are hashed, so we can be sure, that we failed
      else if ((-1 != pos) && mIsHashed) {
        mIsOpen = false;
        ret_val = NULL;
      }

      // we either have sequential stream or didn't find anything, but we can 
      // still try to sequentially go and look for it
      else {
        bool        found = false;
        std::string aux_name;
        std::string aux_name2;

        while ((!found) && JumpToNextDefinition(aux_name)) {
          if (ProcessMask(rFileName, aux_name, aux_name2)) {
            mIsOpen = true;
            mIsEof  = false;
            found   = true;
            ret_val = this;
          }
        }

        if (!found) {
          mIsOpen = false;
          ret_val = NULL;
        }
      }

      return ret_val;
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    typename _Traits::int_type
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    underflow()
    {
      // we don't do anything if EOF
      if (mIsEof) {
        StreamBufType::setg(&(mLineBuffer.front()), &(mLineBuffer.front()), 
            &(mLineBuffer.front()));
        return _Traits::eof();
      }

      // read from buffer if we can
      if (StreamBufType::gptr() && (StreamBufType::gptr() < StreamBufType::egptr())) {
        return _Traits::not_eof(*StreamBufType::gptr());
      }

      // might happen that stream is in !good state
      if (!mIStream.good()) {
        mIsEof = true;
        StreamBufType::setg(&(mLineBuffer.front()), &(mLineBuffer.front()), 
            &(mLineBuffer.front()));
        return _Traits::eof();
      }

      // fill the line buffer and update my state
      FillLineBuffer();

      // if the whole line is just period or it's eof, declare EOF
      if (mState == OUT_OF_BODY_STATE) {
        mIsEof = true;
        StreamBufType::setg(&(mLineBuffer.front()), &(mLineBuffer.front()), 
            &(mLineBuffer.front()));
        return _Traits::eof();
      }

      // restore the buffer mechanism
      StreamBufType::setg(&(mLineBuffer.front()), &(mLineBuffer.front()),
          &(mLineBuffer.back()) + 1);

      return *StreamBufType::gptr();
    }


  //****************************************************************************
  //****************************************************************************
  template<
    typename _CharT, 
    typename _Traits,
    typename _CharTA,
    typename ByteT,
    typename ByteAT
  > 
    void
    BasicIMlfStreamBuf<_CharT, _Traits, _CharTA, ByteT, ByteAT>::
    FillLineBuffer()
    {
      // reset line buffer
      size_t capacity = mLineBuffer.capacity();
      mLineBuffer.clear();
      mLineBuffer.reserve(capacity);

      // read one line into buffer
      int c;
      while ((c = mIStream.get()) != '\n' && c != _Traits::eof()) {
        mLineBuffer.push_back(c);
      }

      // we want to be able to pass last eol symbol
      if (c == '\n') {
        mLineBuffer.push_back(c);
      }

      // we will decide where we are
      switch (mState) {
        case IN_HEADER_STATE:

        case OUT_OF_BODY_STATE:
          if (mLineBuffer[0] != '#') {
            mState = IN_TITLE_STATE;
          }
          break;

        case IN_TITLE_STATE:
          if (mLineBuffer[0] == '.' && (mLineBuffer.back() == '\n' || mIStream.eof())) {
            mState = OUT_OF_BODY_STATE;
          }
          else {
            mState = IN_BODY_STATE;
          }
          break;

        case IN_BODY_STATE:
          // period or EOF will end the file
          if (mLineBuffer[0] == '.' && (mLineBuffer.back() == '\n' || mIStream.eof())) {
            mState = OUT_OF_BODY_STATE;
          }
          if (mLineBuffer.size() == 0) {
            mState = OUT_OF_BODY_STATE;
          }
          break;
      }
    }
}; // namespace STK


#endif // STK_MlfStream_tcc
