#include "MlfStream.h"
#include "common.h"


namespace STK
{
  //******************************************************************************
  LabelContainer::
  ~LabelContainer()
  {
    while (!this->mLabelList.empty())
    {
      delete this->mLabelList.back();
      this->mLabelList.pop_back();
    }
  }

  //******************************************************************************
  size_t
  LabelContainer::
  DirDepth(const std::string & rPath)
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
  LabelContainer::
  Insert(const std::string &  rLabel,
         std::streampos  Pos)
  {
    LabelRecord     ls;
    size_t          depth;
    LabelRecord     tmp_ls;

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
    ls.mStreamPos       = Pos;
    ls.miLabelListLimit = mLabelList.end();


    if (mLabelList.begin() != mLabelList.end()) {
      ls.miLabelListLimit--;
    }

    // if no wildcard chars, then we try to store in hash, otherwise store in 
    // list
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
        std::cerr << "label: " << MatchedPattern() ;
        std::cerr << std::endl;
      }
    }
    else
    {
      this->mLabelList.push_back(new std::pair<std::string,LabelRecord>(rLabel, ls));
    }
  };


  //******************************************************************************
  bool
  LabelContainer::
  FindInHash(const std::string & rLabel, LabelRecord & rLS)
  {
    //bool run   = true;
    bool found = false;

    std::string str;

    // current depth within the str
    DepthType  current_depth    = MAX_LABEL_DEPTH;

    // current search position within the str
    size_t     prev             = rLabel.size() + 1;

    // we will walk through the set depts bacwards so we begin at the end and move
    // to the front...
    std::set<DepthType>::reverse_iterator ri    (this->mDepths.end());
    std::set<DepthType>::reverse_iterator rlast (this->mDepths.begin());
    LabelHashType::iterator               lab;

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
    } // while (run)

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
  LabelContainer::
  FindInList(const std::string & rLabel, LabelRecord & rLS, bool limitSearch)
  {

    bool                      found = false;
    std::string                    str;
    LabelListType::iterator   lab   = mLabelList.begin();
    LabelListType::iterator   limit;

    if (limitSearch && (rLS.miLabelListLimit != mLabelList.end()))
    {
      limit = rLS.miLabelListLimit;
      limit++;
    }
    else
    {
      limit = this->mLabelList.end();
    }

    // we perform sequential search until we run to the end of the list or we find
    // something
    while ((!found) && (lab != limit))
    {
      if (ProcessMask(rLabel, (*lab)->first, str))
      {
        found = true;
      }
      else
      {
        lab++;
      }
    } // while (run)

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
  LabelContainer::
  Find(const std::string & rLabel, LabelRecord & rLS)
  {
    // try to find the label in the Hash
    if (FindInHash(rLabel, rLS))
    {
      // we look in the list, but we limit the search.
      FindInList(rLabel, rLS, true);
      return true;
    } //if (this->mLabelContainer.FindInHash(rLabel, label_stream))
    else
    {
      // we didn't find it in the hash so we look in the list
      return FindInList(rLabel, rLS);
    }
  }

}; // namespace STK

