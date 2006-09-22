/***************************************************************************
 *   copyright           : (C) 2004-2005 by Lukas Burget,UPGM,FIT,VUT,Brno *
 *   email               : burget@fit.vutbr.cz                             *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef STK_SigP_h
#define STK_SigP_h

#include "common.h"
#include "Matrix.h"

namespace STK
{
  /**
    * @brief Generates matrix of @nRepetitions DCT transforms
    * @param rMat reference to matrix
    * @param basesNum number of bases for each transform
    * @param vctLen vector length
    * @param nRepetitions number of repettition 
    * @param includeC0 include zeroth base component (C0)
    */
  void GenerateDCTMatrix(Matrix<FLOAT> &rMat, const size_t nBasis, const size_t vctLen, const size_t nRepetitions, bool includeC0);

  /**
    * @brief Generates constant matrix
    * @param rMat reference to matrix
    * @param nRows number of rows
    * @param nCols number of columns
    * @param constant floating point constant value
    */
  void GenerateConstantMatrix(Matrix<FLOAT> &rMat, const size_t nRows, const size_t nCols, const FLOAT constant);

  /**
    * @brief Generates diagonal matrix
    * @param rMat reference to matrix
    * @param nRows number of rows
    * @param nCols number of columns
    * @param val value on diagonal
    */
  void GenerateDiagMatrix(Matrix<FLOAT> &rMat, const size_t nRows, const size_t nCols, const FLOAT val);  
  
  /**
    * @brief Generates random matrix
    * @param rMat reference to matrix
    * @param nRows number of rows
    * @param nCols number of columns
    * @param minVal minimal value
    * @param maxVal maximal value
    * @param seed for random generator
    */
  void GenerateRandomMatrix(Matrix<FLOAT> &rMat, const size_t nRows, const size_t nCols, const FLOAT minVal, const FLOAT maxVal, const unsigned int seed);

  /**
    * @brief Generates vector of @nRepetitions Hamming windows
    * @param rVct reference to vector
    * @param partLen length of one Hamming window
    * @param nRepetitions number of repettition 
    */  
  void GenerateHammingWindow(BasicVector<FLOAT> &rVct, const size_t partLen, const size_t nRepetitions);  

  /**
    * @brief Generates random window
    * @param rVct reference to vector
    * @param vctLen length of the vector
    * @param minVal minimal value
    * @param maxVal maximal value
    * @param seed for random generator
    */  
  void GenerateRandomWindow(BasicVector<FLOAT> &rVct, const size_t vctLen, const FLOAT minVal, const FLOAT maxVal, unsigned int seed);

  /**
    * @brief Generates constant window
    * @param rVct reference to vector
    * @param vctLen vector length
    * @param constant constant floating point value
    */  
  void GenerateConstantWindow(BasicVector<FLOAT> &rVct, const size_t vctLen, const FLOAT constant);

  /**
    * @brief Generates vector of @nRepetitions LinSpace windows
    * @param rVct reference to vector
    * @param windowLen length of one window
    * @param startValue start value (first value of window)
    * @param endValue end value (last value of window)
    * @param nRepetitions number of repettition 
    */  
  void GenerateLinSpaceWindow(BasicVector<FLOAT> &rVct, const size_t windowLen, const FLOAT startVal, const FLOAT endVal, const size_t nRepetitions);

  /**
    * @brief Generates vector of @nRepetitions Triangular windows
    * @param rVct reference to vector
    * @param windowLen length of one window
    * @param nRepetitions number of repettition 
    */  
  void GenerateTriangWindow(BasicVector<FLOAT> &rVct, const size_t windowLen, const size_t nRepetitions);

}; // namespace STK

#endif
