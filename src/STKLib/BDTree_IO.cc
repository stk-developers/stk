/***************************************************************************
 *   copyright            : (C) 2004 by Lukas Burget,UPGM,FIT,VUT,Brno     *
 *   email                : burget@fit.vutbr.cz                            *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "BDTree.h"
#include <cstring>
#include <stdexcept>

namespace STK
{
  const char* BDTreeHeader::INTRO = "#STKBTREE "; 


  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Read(std::istream& rStream, BDTreeHeader& rHeader)
  {
    DOUBLE_64 aux_double;
    INT_32    i, j;
    INT_32    size;

    // read training data mass for this node (leaf)
    rStream.read(reinterpret_cast<char*>(&aux_double), sizeof(aux_double));
    mN = aux_double;

    // read number of records for this distribution
    rStream.read(reinterpret_cast<char*>(&size), sizeof(size));

    // In case we use binary format for LM
    if(rHeader.mFileVersion == 1)
    {
      DOUBLE_64 z = 0;
      mVec.clear();
      mVec.reserve(rHeader.mVocabSize);

      //bool current_zero = true;
      // read records
      for (i = 0; i < size; ++i) 
      {
        rStream.read(reinterpret_cast<char*>(&aux_double), sizeof(aux_double));

        if(aux_double > 0 && aux_double <= 1)
          mVec.push_back(aux_double);
        // singular null is kept as null (not to mix with prob 1)
        else if(aux_double == 0)
            mVec.push_back(z);
        else
	{
          for(j = 0; j < aux_double; j++)
            mVec.push_back(z);
	}
      }
    }
    else
    {
      assert(size == rHeader.mVocabSize);

      // reserve place
      mVec.clear();
      mVec.reserve(size);

      // read all records
      for (i = 0; i < size; ++i) {
        rStream.read(reinterpret_cast<char*>(&aux_double), sizeof(aux_double));
        mVec.push_back(aux_double);
      }
    }
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  VecDistribution::
  Write(std::ostream& rStream, BDTreeHeader& rHeader)
  {
    VecDistribution::Container::const_iterator i;
    DOUBLE_64 aux_double;
    INT_32    aux_int;

    // write training data mass for this node (leaf)
    aux_double = mN;
    rStream.write(reinterpret_cast<char*>(&aux_double), sizeof(aux_double));


   // In case we use binary format for LM
    if(rHeader.mFileVersion == 1)
    {
      // First count the size we have
      // For the distribution, if we have several zeroes in the row - we encode this as an integer: 0 0 0 -> 3. Probabilities stay
      DOUBLE_64 prev_i = 1;
      int counter = 0; 
      for (i = mVec.begin(); i != mVec.end(); ++i) 
      {
        aux_double = *i;
        if(aux_double != 0 || (aux_double == 0 && prev_i != 0))
          counter++;

        prev_i = aux_double;
      }
      // write number of records in this distribution
      aux_int = counter;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      // write set records
      counter = 0;
      prev_i = 0;
      for (i = mVec.begin(); i != mVec.end(); ++i) 
      {
        aux_double = *i;

        if(aux_double > 0)
	{
          if(prev_i == 0 && i != mVec.begin())
	  {
            // we leave singular null as null
            if(counter == 1)
              prev_i = (DOUBLE_64)0;
            else
              prev_i = (DOUBLE_64)counter;

            rStream.write(reinterpret_cast<char*>(&prev_i), sizeof(prev_i));
            counter = 0;
          }

          rStream.write(reinterpret_cast<char*>(&aux_double), sizeof(aux_double));
	}
        else
          counter++;

        if(i == mVec.end()-1 && aux_double == 0)
	{
          // we leave singular null as null
          if(counter == 1)
            aux_double = (DOUBLE_64)0;
          else
            aux_double = (DOUBLE_64)counter;

          rStream.write(reinterpret_cast<char*>(&aux_double), sizeof(aux_double));
	}

        prev_i = aux_double;
      }
    }
    else
    {
      // write number of records in this distribution
      aux_int = mVec.size();
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      for (i = mVec.begin(); i != mVec.end(); ++i) {
        aux_double = *i;
        rStream.write(reinterpret_cast<char*>(&aux_double), sizeof(aux_double));
      }
    }
  }



  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  void
  BDTreeHeader::
  Write(std::ostream& rStream)
  {
    //char    aux_char;
    INT_32  aux_int;
    std::streampos stream_pos;

    stream_pos = rStream.tellp();

    rStream.write(BDTreeHeader::INTRO, 10);
    rStream.put(mBinary ? 'b' : 'a');

    mHeaderSize = (6*sizeof(INT_32));

    if (mBinary) {
      aux_int = mHeaderSize;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      
      aux_int = mFileVersion;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      
      aux_int = mOrder;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      
      aux_int = mVocabSize;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      aux_int = mExtra0;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      aux_int = mExtra1;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
    }
    else {
      throw std::runtime_error("ASCII format not supported yet");
    }
  }
  //***************************************************************************/
  // Write in binary format for LVCSR. In that case Target and Predictor vocabs 
  // may differ, so we write both sizes
  //***************************************************************************/

  void
  BDTreeHeader::
  Write_bin1(std::ostream& rStream)
  {
    //char    aux_char;
    INT_32  aux_int;
    std::streampos stream_pos;

    stream_pos = rStream.tellp();

    rStream.write(BDTreeHeader::INTRO, 10);
    rStream.put(mBinary ? 'b' : 'a');

    mHeaderSize = 7*sizeof(INT_32);

    if (mBinary) {
      aux_int = mHeaderSize;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      
      aux_int = mFileVersion;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      
      aux_int = mOrder;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      
      aux_int = mVocabSize;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      aux_int = mPredictorVocabSize;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      aux_int = mExtra0;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      aux_int = mExtra1;
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
    }
    else {
      throw std::runtime_error("ASCII format not supported yet");
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTreeHeader::
  Read(std::istream& rStream)
  {
    char    aux_char;
    INT_32  aux_int;
    char    p_intro[strlen(BDTreeHeader::INTRO)];

    // read intro and check it
    rStream.read(p_intro, strlen(BDTreeHeader::INTRO));
    if (strcmp(p_intro, BDTreeHeader::INTRO)) {
    }
    
    // binary or ascii???
    rStream.get(aux_char);
    switch(aux_char) {
      case 'a': mBinary = false; break;
      case 'b': mBinary = true; break;
      default : throw std::runtime_error("Format not supported");
    }

    if (mBinary) {
      rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      mHeaderSize  = aux_int;

      rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      mFileVersion = aux_int;

      rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      mOrder       = aux_int;

      rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      mVocabSize   = aux_int;

      if(mFileVersion == 1)
      {
        rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
        mPredictorVocabSize   = aux_int;
      }
      else {
        mPredictorVocabSize   = mVocabSize;
      }

      rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      mExtra0      = aux_int;

      rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      mExtra1      = aux_int;
    }
    else { // if (Binary) 
      // TODO: Implement this
      throw std::runtime_error("ASCII format not supported yet");
    }
  }


  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  LoadFile(const std::string& rName)
  {
    //OStkStream model_stream(model_name.c_str());
    //if (!model_stream.good()) {
    //  throw runtime_error(string("Error opening output file ") + 
    //      model_name);
    //}

    //BDTreeHeader file_header;
    //file_header.mBinary = true;
    //file_header.mOrder  = 

    //// output header
    //file_header.Write(model_stream);
    //// output model
    //this->Write(model_stream, file_header);

    //model_stream.close();
  }

  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  SaveFile(const std::string& rName)
  {
  }
  
  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  Read(std::istream& rStream, BDTreeHeader& rHeader)
  {
    INT_32 flag = 0;

    // read main flag
    rStream.read(reinterpret_cast<char*>(&flag), sizeof(flag));
    
    if (flag & 1) {
      mpQuestion  = new BSetQuestion(rStream, rHeader);
      mpTree0     = new BDTree(rStream, rHeader);
      mpTree1     = new BDTree(rStream, rHeader);
      mpDist      = NULL;
      mpBackoffDist = NULL;
    }
    else {
      mpDist      = new VecDistribution(rStream, rHeader);
      mpQuestion  = NULL;
      mpTree0     = NULL;
      mpTree1     = NULL;
      mpBackoffDist = NULL;
    }
  } // Read(std::istream& rStream, BDTreeHeader& rHeader)


  //***************************************************************************/
  //***************************************************************************/
  void
  BDTree::
  Write(std::ostream& rStream, BDTreeHeader& rHeader)
  {
    INT_32 flag = 0;

    if (!IsLeaf()) {
      flag |= 1;
    }

    // write main flag
    rStream.write(reinterpret_cast<char*>(&flag), sizeof(flag));

    if (IsLeaf()) {
      mpDist->Write(rStream, rHeader);
    }
    else {
      mpQuestion->Write(rStream, rHeader);
      mpTree0->Write(rStream, rHeader);
      mpTree1->Write(rStream, rHeader);
    }
  } // Write(std::ostream& rStream, BDTreeHeader& rHeader)


  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  //***************************************************************************/
  void
  BSetQuestion::
  Read(std::istream& rStream, BDTreeHeader& rHeader)
  {
    INT_32 aux_int;
    INT_32 i, j, k;
    INT_32 size;

    // which predictor
    rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
    mPred = aux_int;

    // number of records for set
    rStream.read(reinterpret_cast<char*>(&size), sizeof(size));

    // In case we use binary format for LM
    if(rHeader.mFileVersion == 1)
    {
      mSet.clear();
      mSet.reserve(rHeader.mPredictorVocabSize);

      bool current_zero = true;
      // read records
      for (i = 0; i < size; ++i) 
      {
        rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
        // there can be a leading "0" - we should just skip it
        if(i == 0 && aux_int == 0)
          ;
        else
	{
          if(current_zero)
            k = 0;
          else
            k = 1;

          for(j = 0; j < aux_int; j++)
            mSet.push_back(k);
	}

        if(current_zero)
          current_zero = false;
        else
          current_zero = true;
      }
    }
    else
    {
      assert(size == rHeader.mPredictorVocabSize);

      mSet.clear();
      mSet.reserve(size);

      // read records
      for (i = 0; i < size; ++i) {
        rStream.read(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
        mSet.push_back(aux_int);
      }
      }
  }


  //***************************************************************************/
  //***************************************************************************/
  void

  BSetQuestion::
  Write(std::ostream& rStream, BDTreeHeader& rHeader)
  {
    BSetQuestion::SetType::const_iterator i;
    INT_32 aux_int;

    // write predictor
    aux_int = mPred;
    rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

    // In case we use binary format for LM
    if(rHeader.mFileVersion == 1)
    {
      // First count the size we have
      INT_32 prev_i = 0;
      INT_32 counter = 1; // we always start from zero - even if there is none, we'll write "0"
      for (i = mSet.begin(); i != mSet.end(); ++i) 
      {
        aux_int = *i ? 1 : 0;
        if(prev_i != aux_int)
	{
          counter++;
          prev_i = aux_int;
	}
      }
      aux_int = counter;
      // write set size
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      // write set records
      counter = 0;
      prev_i = 0;
      for (i = mSet.begin(); i != mSet.end(); ++i) 
      {
        aux_int = *i ? 1 : 0;
        if(prev_i == aux_int)
          counter++;
        else
	{
          prev_i = counter;
          rStream.write(reinterpret_cast<char*>(&prev_i), sizeof(prev_i));
          counter = 1;
          prev_i = aux_int;
	}

        if(i == mSet.end()-1)
	{
          aux_int = counter;
          rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
	}
      }

    }
    else
    {
      aux_int = mSet.size();
      // write set size
      rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));

      // write set records
      for (i = mSet.begin(); i != mSet.end(); ++i) {
        aux_int = *i ? 1 : 0;
        rStream.write(reinterpret_cast<char*>(&aux_int), sizeof(aux_int));
      }
      }
  }
} //namespace STK
