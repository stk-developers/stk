#include "labelreader.h"
#include "common.h"
#include <iostream>
#include <stdexcept>

using namespace std;
namespace STK
{
  //******************************************************************************
  MlfDefType MlfDefinition(const string & line)
  {
    // TODO
    // we are only looking for the -> and => symbols, but we don't really parse
    // line now, so next step would be to really parse the line
    // TODO

    size_t pos = 0;

    while (((pos = line.find_first_of("=-", pos)) != line.npos) ||
           (pos < line.size()))
    {
      if (line[pos++] == '>')
        return MLF_DEF_SUB_DIR_DEF;
    }
    return MLF_DEF_IMMEDIATE_TRANSCRIPTION;
  } // MlfDefType MlfDefinition(const string & line)

  //******************************************************************************
  LabelStore::
  ~LabelStore()
  {
    while (!this->mLabelList.empty())
    {
      delete this->mLabelList.back();
      this->mLabelList.pop_back();
    }
  }

  //******************************************************************************
  size_t
  LabelStore::
  DirDepth(const string & rPath)
  {
    size_t depth     = 0;
    size_t length    = rPath.length();
    const char * s   = rPath.c_str();

    for (size_t i = 0; i < length; i++)
    {
      if (*s == '/' || *s == '\\')
      {
        depth++;
      }
      s++;
    }
    return depth;
  }


  //******************************************************************************
  void
  LabelStore::
  Insert(const string &  rLabel,
         istream *       pStream,
         std::streampos  Pos)
  {
    LabelStream     ls;
    size_t          depth;
    LabelStream     tmp_ls;

    // we need to compute the depth of the label path if
    // wildcard is used
    // do we have a wildcard???
    if (rLabel[0] == '*')
    {
      depth = this->DirDepth(rLabel);
    }
    else
    {
      depth = MAX_LABEL_DEPTH;
    }

    // perhaps we want to store the depth of the path in the label for the wildcards
    // to work
    this->mDepths.insert(depth);

    // store the values
    ls.mpStream       = pStream;
    ls.mStreamPos     = Pos;
    ls.miLabelListLimit = mLabelList.end();


    if (mLabelList.begin() != mLabelList.end())
      ls.miLabelListLimit--;

    if (rLabel.find_first_of("*?%",1) == rLabel.npos)
    {
      if (!Find(rLabel, tmp_ls))
      {
        // look in the
        this->mLabelMap[rLabel] = ls;
      }
      else
      {
        // TODO
        // we want to print out some warning here, that user is adding
        // a label which will never be found because more general
        // definition was found in the list before
        // TODO
        std::cerr << "Warning... more general definition found when inserting " << rLabel << " ... ";
        std::cerr << "file: " << dynamic_cast<istkstream *>(tmp_ls.mpStream)->name() << " ";
        std::cerr << "label: " << MatchedPattern() ;
        std::cerr << std::endl;
      }
    }
    else
    {
      this->mLabelList.push_back(new std::pair<string,LabelStream>(rLabel, ls));
    }
  };

  //******************************************************************************
  void
  LabelStore::
  Insert(const string & rLabel, istream * pStream)
  {
    this->Insert(rLabel, pStream, pStream->tellg());
  };


  //******************************************************************************
  bool
  LabelStore::
  FindInHash(const string & rLabel, LabelStream & rLS)
  {
    bool run   = true;
    bool found = false;

    string str;

    // current depth within the str
    DepthType  current_depth    = MAX_LABEL_DEPTH;

    // current search position within the str
    size_t     prev             = rLabel.size() + 1;

    cerr << "mDepths.size = " << mDepths.size() << endl;
    // we will walk through the set depts bacwards so we begin at the end and move
    // to the front...
    std::set<DepthType>::reverse_iterator ri    (this->mDepths.end());
    std::set<DepthType>::reverse_iterator rlast (this->mDepths.begin());
    LabelHashType::iterator               lab;

    cout << "Searchin in Hash.." << endl;

    // we perform the search until we run to the end of the set or we find something
    while ((!found) && (ri != rlast))
    {
      // we don't need to do anything with the string if the depth is set to
      // max label depth since it contains no *
      if (*ri == MAX_LABEL_DEPTH)
      {
        found = ((lab=this->mLabelMap.find(rLabel)) != this->mLabelMap.end());
        if (found) str = rLabel;
      }
      // we will crop the string and put * in the begining and try to search
      else
      {
        // we know that we walk backwards in the depths, so we need to first find
        // the last / and
        if (current_depth == MAX_LABEL_DEPTH)
        {
          if (*ri > 0)
          {
            // we find the ri-th / from back
            for (DepthType i=1; (i <= *ri) && (prev != rLabel.npos); i++)
            {
              prev = rLabel.find_last_of("/\\", prev-1);
            }
          }
          else
          {
            prev = 0;
          }

          // check if finding succeeded (prev == str.npos => failure, see STL)
          if (prev != rLabel.npos)
          {
            // construct the new string beign sought for
            str.assign(rLabel, prev, rLabel.size());
            str = '*' + str;

            // now we try to find
            found = ((lab=this->mLabelMap.find(str)) != this->mLabelMap.end());

            // say, that current depth is *ri
            current_depth = *ri;
          }
          else
          {
            prev = rLabel.size() + 1;
          }
        }     // if (current_depth == MAX_LABEL_DEPTH)
        else
        {
          // now we know at which / we are from the back, so we search forward now
          // and we need to reach the ri-th /
          while (current_depth > *ri)
          {
            // we try to find next /
            if ((prev = rLabel.find_first_of("/\\", prev+1)) != rLabel.npos)
              current_depth--;
            else
              return false;
          }

          // construct the new string beign sought for
          str.assign(rLabel, prev, rLabel.size());
          str = '*' + str;

          // now we try to find
          found = ((lab=this->mLabelMap.find(str)) != this->mLabelMap.end());
        }
      }

      // move one element further (jump to next observed depth)
      ri++;
    } // while(run)

    // some debug info
    if (found)
    {
      rLS                   = lab->second;
      this->mMatchedPattern = str;
    }

    return found;
  }


  //******************************************************************************
  bool
  LabelStore::
  FindInList(const string & rLabel, LabelStream & rLS, bool limitSearch)
  {

    bool                      found = false;
    string                    str;
    LabelListType::iterator   lab   = mLabelList.begin();
    LabelListType::iterator   limit;

    cerr <<(int)limitSearch << "Searchin in List.. "  <<  "<<" << endl;

    if (limitSearch && (rLS.miLabelListLimit != NULL))
    {
      cout << "Limiting search... " << endl;
      limit = rLS.miLabelListLimit;
      limit++;
    }
    else
    {
      limit = this->mLabelList.end();
    }

    cerr <<(int)limitSearch << "Searchin in List.. "  <<  "<<" << endl;

    // we perform sequential search until we run to the end of the list or we find
    // something
    while ((!found) && (lab != limit))
    {
      if (ProcessMask(rLabel, (*lab)->first, str))
      {
        found               = true;
      }
      else
      {
        lab++;
      }
    } // while(run)

    // some debug info
    if (found)
    {
      rLS                       = (*lab)->second;
      this->mMatchedPattern     = (*lab)->first;
      this->mMatchedPatternMask = str;
    }
    return found;
  }


  //******************************************************************************
  bool
  LabelStore::
  Find(const string & rLabel, LabelStream & rLS)
  {
    // try to find the label in the Hash
    if (FindInHash(rLabel, rLS))
    {
      cerr << "Found in Hash in  mLabelStore.Find ... ("<< rLabel << ",...) " << endl;

      // we look in the list, but we limit the search.
      FindInList(rLabel, rLS, true);
      return true;
    } //if (this->mLabelStore.FindInHash(rLabel, label_stream))
    else
    {
      cerr << "Not found in Hash in  mLabelStore.Find ... ("<< rLabel << ",...) " << endl;
      // we didn't find it in the hash so we look in the list
      return FindInList(rLabel, rLS);
    }
  }


  //****************************************************************************
  LabelReader::
  ~LabelReader()
  {
    if (mpTempStream != NULL)
    {
      delete mpTempStream;
    }

    // go through the list and delete every stream object = close the streams
    while (!mStreamList.empty())
    {
      delete mStreamList.front();
      mStreamList.pop_front();
    }
  }


  //****************************************************************************
  bool
  LabelReader::
  FindInLabelStore(const string & rLabel, LabelStream & rLS)
  {
    // try to find the label in the Hash
    if (mLabelStore.Find(rLabel, rLS))
    {
      // we look in the list, but we limit the search.
      mpCurStream = rLS.mpStream;
      mpCurStream->seekg(rLS.mStreamPos);
      return true;
    } //if (this->mLabelStore.FindInHash(rLabel, label_stream))
    else
    {
      return false;
    }
  }


  //****************************************************************************
  bool
  LabelReader::
  FindInStreamList(const string & rLabel, LabelStream & rLS)
  {
    bool                       found = false;
    string                     line;
    string                     tmp_line;
    string                     mask_string;
    bool                       skip_next_reading = false;

    // we will use this enum to move in the MLF parsing automaton states
    LabelReadingState state = STATE_MLF_BEGIN;

    // read the MLFs and parse them
    while (!found && !IsHashed())
    {
      // read the line
      if (!skip_next_reading)
      {
        GetLine(line);
      }

      // next reading should not be skipped
      skip_next_reading = false;

      // finite automaton
      switch (state)
      {
        case STATE_MLF_BEGIN:
        {
          // we consider other than #!MLF!# to be MLFDef string
          if (line != "#!MLF!#")
          {
            // we want to skip next reading as the line is supposed
            // to be label definition
            skip_next_reading = true;
          }
          state             = STATE_MLF_DEF;
          break;
        } // case MLF_BEGIN:

        // MLF definition
        case STATE_MLF_DEF:
        {
          // decide what kind of definition we encountered
          // we only want to hash immediate definition
          if (MlfDefinition(line) == MLF_DEF_IMMEDIATE_TRANSCRIPTION)
          {
            // we store the string in internal mode, not HTK format
            ParseHTKString(line, tmp_line);
            // insert into the storage
            mLabelStore.Insert(tmp_line, mpCurStream);

            if (found = ProcessMask(rLabel, tmp_line, mask_string))
            {
              mMatchedPattern     = tmp_line;
              mMatchedPatternMask = mask_string;
            }

            // chante state
            state = STATE_MLF_BODY;
          }
          break;
        }

        case STATE_MLF_BODY:
        {
          // reading '.' means go to next MLF_DEF
          if (line == ".")
          {
            state = STATE_MLF_DEF;
          }
          break;
        }
      } // switch (state)
    }
    return found;
  }


  //****************************************************************************
  bool
  LabelReader::
  OpenLabelFile(const string & rLabel)
  {
    // we try to open the file
    mpTempStream = new istkstream;
    dynamic_cast<istkstream *>(mpTempStream)->open(rLabel);

    if (!mpTempStream->good())  // on error return false and clear the object
    {
      delete mpTempStream;
      mpTempStream = NULL;
      return false;
    }
    else                        // on success return true and set current
                                // stream
    {
      mpCurStream = mpTempStream;
      return true;
    }

  } // OpenLabelFile(const string & lFile)


  //****************************************************************************
  bool
  LabelReader::
  Open(const string & rLabel)
  {
    // helping strcuture when seeking in the LabelStore
    LabelStream       label_stream;
    bool              found = false;
    string            buffer;

    // it is possilble, that we look for a label that is currently
    // mpTempStream (only for istkstream)
    if ((mpTempStream != NULL) && (typeid(mpTempStream) == typeid(istkstream)))
    {
      // if we want to read from the same label, we only rewind
      if (dynamic_cast<istkstream *>(mpTempStream)->name() == rLabel)
      {
        // rewind the stream to the begining and read again
        mpTempStream->seekg(0, ios::beg);
        return true;
      }
      // else we delete the object and continue searching
      else
      {
        delete mpTempStream;
        mpTempStream = NULL;
      }
    }  // if (mpTempStream != NULL)

    // try to find the label in the Hash
    if (FindInLabelStore(rLabel, label_stream))
    {
      cerr << "Found in LabelStore ... " << endl;
      this->mpCurStream = label_stream.mpStream;
      this->mpCurStream->seekg(label_stream.mStreamPos);
      found = true;
    } //if (this->mLabelStore.FindInHash(rLabel, label_stream))
    else
    {
      // neither in hash nor in list => if all MLF's hashed, find file
      // else keep searching in MLF's and then find file
      if (IsHashed())
      {
        found = this->OpenLabelFile(rLabel);
      } // if (this->mIsHashed)
      else
      {
        if (!(found = FindInStreamList(rLabel, label_stream)))
        {
          found = this->OpenLabelFile(rLabel);
        }
      } // else (this->mIsHashed)
    } // else (FindInLabelStore(rLabel, label_stream))
    return found;
  }; // Open(const string & rLabel)


  //****************************************************************************
  bool
  LabelReader::
  RegisterMLF(const string & rFName)
  {
    istkstream *   new_stream = new istkstream;

    // open the MLF
    new_stream->open(rFName);

    // if OK, we add it to the MLF list
    if (new_stream->good())
    {
      // add to the list
      mStreamList.push_back(new_stream);
      return true;
    }    // if (new_stream.good())
    else // else we delete the object and trow an error
    {
      delete new_stream;
      throw std::logic_error("Cannot open " + rFName);
      return false;
    }
  }; //RegisterMLF


  //****************************************************************************
  bool
  LabelReader::
  GetLine(std::string & rLine)
  {
    bool ret;

    // we start reading
    if (this->mpCurStream == NULL)
    {
      if (!this->mStreamList.empty())
      {
        this->miCurStream = this->mStreamList.begin();
        this->mpCurStream = *(this->miCurStream);
      }
      else
      {
        return false;
      }
    }

    // if we encountered end of stream, we switch to the next stream in list
    if (this->mpCurStream->eof())
    {
      // move to the next stream
      (this->miCurStream)++;

      if (this->miCurStream == mStreamList.end())
      {
        mIsHashed = true;
        return false;
      }
      else
      {
        this->mpCurStream = *(this->miCurStream);
      }
    }

    // read the line here. we check for eof, EOL (end of label), etc
    if (std::getline(*(this->mpCurStream), rLine))
    {
      if (this->mpCurStream->eof() || (rLine == "."))
      {
        this->mStateFlags |= LabelReader::FLAG_EOL;

        if (this->mpCurStream->eof())
        {
          // move to the next stream
          (this->miCurStream)++;

          if (this->miCurStream == mStreamList.end())
          {
            mIsHashed = true;
          }
          else
          {
            mpCurStream = *(miCurStream);
            mpTempStream->seekg(0, ios::beg);
          }
        }
      }
      else
      {
        this->mStateFlags &= !(LabelReader::FLAG_EOL);
      }
      return true;
    }
    else
    {
      return false;
    }
  }


//****************************************************************************
  void
  LabelReader::
  HashAll()
  {
    // we will use this to iterate through the stream list
    StreamListType::iterator   cur_stream_it = this->mStreamList.begin();
    istream *                  cur_stream;
    std::streampos             cur_stream_pos;
    string                     line;
    string                     tmp_line;
    bool                       skip_next_reading = false;

    // we will use this enum to move in the MLF parsing automaton states
    enum
    {
      MLF_BEGIN,     // we begin parsing
      MLF_DEF,       // this is where label block is defined
      MLF_BODY,
    } state = MLF_BEGIN;

    // Go through the streamlist and parse the files
    while (cur_stream_it != this->mStreamList.end())
    {
      cur_stream     = *cur_stream_it;
      cur_stream_pos = cur_stream->tellg();

      // parse the file
      while (!cur_stream->eof())
      {
        // get current line position
        cur_stream_pos = cur_stream->tellg();

        // read the line
        if (!skip_next_reading)
        {
          std::getline(*cur_stream, line);
        }

        // next reading should not be skipped
        skip_next_reading = false;

        // finite automaton
        switch (state)
        {
          case MLF_BEGIN:
          {
            // we consider other than #!MLF!# to be MLFDef string
            if (line != "#!MLF!#")
            {
              // we want to skip next reading as the line is supposed
              // to be label definition
              skip_next_reading = true;
            }
            state             = MLF_DEF;
            break;
          } // case MLF_BEGIN:

          // MLF definition
          case MLF_DEF:
          {
            // decide what kind of definition we encountered
            // we only want to hash immediate definition
            if (MlfDefinition(line) == MLF_DEF_IMMEDIATE_TRANSCRIPTION)
            {
              // we store the string in internal mode, not HTK format
              ParseHTKString(line, tmp_line);
              // insert into the storage
              this->mLabelStore.Insert(tmp_line, cur_stream);
              // chante state
              state = MLF_BODY;
            }
            break;
          }

          case MLF_BODY:
          {
            // reading '.' means go to next MLF_DEF
            if (line == ".")
            {
              state = MLF_DEF;
            }
            break;
          }
        } // switch (state)
      }



      // move to the next stream
      cur_stream_it++;
    } // while (cur_stream_it != this->mStreamList.end())

    // if we got to this point, no error occured so say that we have Hashed all
    this->mIsHashed = true;
  }; // HashAll(..)


}; // namespace STK

