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
      pos_type pos = mOStream.tellp();

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
    : mIsOpen(false), mIsHashed(false), mIStream(rIStream), mLineBuffer(), 
      mIsEof(true)
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
      pos_type orig_pos = mIStream.tellg();

      // for streams like stdin, pos will by definition be -1, so we can only 
      // rely on sequential access and cannot hash it.
      if (-1 != orig_pos) {
        std::basic_string<_CharT, _Traits> line_buffer;
        std::basic_string<_CharT, _Traits> name;
        bool        good_to_go = mIStream.good();
        pos_type    pos;

        // the reader will be a finite state automaton. We assume we are in the 
        // initial state, at the "#!MLF!# record
        int         state = 0;

        // continue while we get to the end of the line
        while (good_to_go) {
          if (std::getline(mIStream, line_buffer)) {
            switch (state) {
              // reading header
              case 0: 
                // TODO: Probably add option to require the header and do
                // something here

              // reading label title
              case 1:
                // skip out the comments, especially the #!MLF!#
                if (line_buffer[0] == '#') {
                  continue;
                }
                else {
                  pos = mIStream.tellg();
                  ::ParseHTKString(line_buffer, name);
                  mLabels.Insert(name, pos);

                  // switch state to reading labels
                  state = 2;
                }
                break;

              // reading label body
              case 2:
                // skip all lines
                if (line_buffer == ".") {
                  state = 1;
                }
                continue;
                break;

            } // switch
          } // if (std::getline(mIStream, line_buffer)) 
          else { 
            good_to_go = false;
          }
        } // while (good_to_go)

        // move to the original position
        mIStream.clear();
        mIStream.seekg(orig_pos);
        mIsHashed = true;
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
      mIsEof = true;

      if (!mIsOpen) {
        return NULL;
      }
      else {
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
      // this behavior is compatible with ifstream
      if (mIsOpen) {
        Close();
        return NULL;
      }

      // retreive position
      pos_type pos = mIStream.tellg();

      // for streams like stdin, pos will by definition be -1, so we can only 
      // rely on sequential access. At this place, we decide what to do
      if ((-1 == pos) || (!mIsHashed)) {
      }
      // this else = we DO have random access stream and it has been hashed
      else { 
        LabelRecord label_record;
        if (mLabels.Find(rFileName, label_record)) {
          mIStream.seekg(label_record.mStreamPos);

          // we don't want the other stream to be bad, so we transfer the 
          // flagbits to this stream
          if (!mIStream.good()) {
            mIStream.clear();
            mIsOpen = false;
            return NULL;
          }
          else {
            mIsOpen = true;
            mIsEof = false;
            return this;
          }
        }
        else {
          mIsOpen = false;
          return NULL;
        }
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

      // maybe we didn't read anything
      if (mLineBuffer.size() == 0) {
        mIsEof = true;
        StreamBufType::setg(&(mLineBuffer.front()), &(mLineBuffer.front()), 
            &(mLineBuffer.front()));
        return _Traits::eof();
      }

      // if the whole line is just period, declare EOF
      if (mLineBuffer[0] == '.' && (c == '\n' || mIStream.eof())) {
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

}; // namespace STK


#endif // STK_MlfStream_tcc
