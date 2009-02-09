#ifndef STK_StkStream_tcc
#define STK_StkStream_tcc

#include "common.h"
#include <cstring>
#include <iostream>
#pragma GCC system_header

namespace STK
{
  /*
  char * ExpandFilterCommand(const char *command, const char *filename)
  {

    char *out, *outend;
    const char *chrptr = command;
    int ndollars = 0;
    int fnlen = strlen(filename);

    while (*chrptr++) ndollars += (*chrptr ==  *gpFilterWldcrd);

    out = (char*) malloc(strlen(command) - ndollars + ndollars * fnlen + 1);
    //if (out == NULL) Error("Insufficient memory");

    outend = out;

    for (chrptr = command; *chrptr; chrptr++) {
      if (*chrptr ==  *gpFilterWldcrd) {
        strcpy(outend, filename);
        outend += fnlen;
      } else {
        *outend++ = *chrptr;
      }
    }
    *outend = '\0';
    return out;
  };
  */


  //******************************************************************************
  template<
    typename _CharT, 
    typename _Traits
  > 
  basic_stkbuf<_CharT, _Traits> *
  basic_stkbuf<_CharT, _Traits>::
  close(void)
  {
    // we only want to close an opened file
    if (this->is_open())
    {
      // we want to call the parent close() procedure
      std::basic_filebuf<_CharT, _Traits>::close();

      // and for different stream type we perform different closing
      if (mStreamType == basic_stkbuf::t_file)
      {
        fclose(mpFilePtr);
      }
      else if (mStreamType == basic_stkbuf::t_pipe)
      {
        pclose(mpFilePtr);
      }
      else if (mStreamType == basic_stkbuf::t_stdio)
      {

      }
      
      mpFilePtr     = NULL;
      mFilename    = "";
      mMode        = std::ios_base::openmode(0);
      mStreamType = basic_stkbuf::t_undef;
      return this;
    }
    else
      return 0;
  }


  template<
    typename _CharT, 
    typename _Traits
  > 
  void
  basic_stkbuf<_CharT, _Traits>::
  open_mode(std::ios_base::openmode __mode, int&, int&,  char* __c_mode)
  {
    bool __testb = __mode & std::ios_base::binary;
    bool __testi = __mode & std::ios_base::in;
    bool __testo = __mode & std::ios_base::out;
    bool __testt = __mode & std::ios_base::trunc;
    bool __testa = __mode & std::ios_base::app;

    if (!__testi && __testo && !__testt && !__testa)
      strcpy(__c_mode, "w");
    if (!__testi && __testo && !__testt && __testa)
      strcpy(__c_mode, "a");
    if (!__testi && __testo && __testt && !__testa)
      strcpy(__c_mode, "w");
    if (__testi && !__testo && !__testt && !__testa)
      strcpy(__c_mode, "r");
    if (__testi && __testo && !__testt && !__testa)
      strcpy(__c_mode, "r+");
    if (__testi && __testo && __testt && !__testa)
      strcpy(__c_mode, "w+");
    if (__testb)
      strcat(__c_mode, "b");
  }


  //******************************************************************************
  template<
    typename _CharT, 
    typename _Traits
  > 
  basic_stkbuf<_CharT, _Traits> *
  basic_stkbuf<_CharT, _Traits>::
  open(const char* pFName, std::ios::openmode m, const char* pFilter)
  {
    basic_stkbuf<_CharT, _Traits>* p_ret = NULL;

    if (NULL == pFName)
      return NULL;
      
    // we need to assure, that the stream is not open
    if (!this->is_open())
    {
      char mstr[4] = {'\0', '\0', '\0', '\0'};
      int __p_mode = 0;
      int __rw_mode = 0;

      // now we decide, what kind of file we open
      if (!strcmp(pFName,"-"))
      {
        if      ((m & std::ios::in) && !(m & std::ios::out))
        {
          mpFilePtr   = stdin;
          mMode       = std::ios::in;
          mFilename   = pFName;
          mStreamType = t_stdio;
          p_ret       = this;
        }
        else if ((m & std::ios::out) && !(m & std::ios::in))
        {
          mpFilePtr   = stdout;
          mMode       = std::ios::out;
          mFilename   = pFName;
          mStreamType = t_stdio;
          p_ret       = this;
        }
        else
          p_ret = NULL;
      }
      else if ( pFName[0] == '|' )
      {
        const char* command = pFName + 1;

        if      ((m & std::ios::in) && !(m & std::ios::out)) m = std::ios::in;
        else if ((m & std::ios::out) && !(m & std::ios::in)) m = std::ios::out;
        else return NULL;

        // we need to make some conversion
        // iostream -> stdio open mode string
        this->open_mode(m, __p_mode, __rw_mode, mstr);

        if ((mpFilePtr = popen(command, mstr)))
        {
          mFilename   = command;
          mMode       = m;
          mStreamType = t_pipe;
          p_ret       = this;
        }
        else
          p_ret = 0;
      }
      else
      {
        // maybe we have a filter specified
        if ( pFilter 
        && ('\0' != pFilter[0]))
        {
          char* command = STK::expandFilterCommand(pFilter, pFName);

          if      ((m & std::ios::in) && !(m & std::ios::out)) m = std::ios::in;
          else if ((m & std::ios::out) && !(m & std::ios::in)) m = std::ios::out;
          else return NULL;

          // we need to make some conversion
          // iostream -> stdio open mode string
          this->open_mode(m, __p_mode, __rw_mode, mstr);

          if ((mpFilePtr = popen(command, mstr)))
          {
            mFilename     = pFName;
            mMode         = m;
            mStreamType   = t_pipe;
            p_ret         = this;
          }
          else
            p_ret = 0;
        }
        else // if (!filter.empty())
        {
          // we need to make some conversion
          // iostream -> stdio open mode string
          this->open_mode(m, __p_mode, __rw_mode, mstr);

          if ((mpFilePtr = fopen(pFName, mstr)))
          {
            mFilename   = pFName;
            mMode       = m;
            mStreamType = t_file;
            p_ret       = this;
          }
          else {
            p_ret = NULL;
          }
        }
      }

      // here we perform what the stdio_filebuf would do
      if (p_ret) {
        superopen(mpFilePtr, m);
      }
    } //if (!isopen)

    return p_ret;
  }

  //******************************************************************************
  template<
    typename _CharT, 
    typename _Traits
  > 
  void
  basic_stkbuf<_CharT, _Traits>::
  superopen(std::__c_file* __f, std::ios_base::openmode __mode,
        size_t __size)
  {
    this->_M_file.sys_open(__f, __mode);
    if (this->is_open())
    {
      this->_M_mode = __mode;
      this->_M_buf_size = __size;
      this->_M_allocate_internal_buffer();
      this->_M_reading = false;
      this->_M_writing = false;
      this->_M_set_buffer(-1);
    }
  }
}

// STK_StkStream_tcc
#endif
