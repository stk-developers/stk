//
// C++ Interface: %{MODULE}
//
// Description: 
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef STK_Features_h
#define STK_Features_h

//*****************************************************************************
//*****************************************************************************
// Specific includes
//
#include "common.h"
#include "Matrix.h"
#include "stkstream.h"

//*****************************************************************************
//*****************************************************************************
// Standard includes
//
#include <list>
#include <string>


//*****************************************************************************
//*****************************************************************************
// Code
//

namespace STK
{
  /** *************************************************************************
   * @brief 
   */
  class FeatureRepository
  {
  private:
    /// List (queue) of input feature files
    std::list<FileListElem>     mInputQueue;
    
    /// current stream
    IStkStream                  mStream;
      
  public:
    /**
     * @brief Adds a single feature file to the repository
     * @param rFileName file to read features from
     */
    void
    AddFile(const std::string & rFileName);
    
    /**
     * @brief Adds a list of feature files to the repository
     * @param rFileName feature list file to read from
     */
    void
    AddFileList(const std::string & rFileName, const std::string & rFilter = "");
  
    /**
     * @brief Returns current feature file logical name
     * @return 
     */
    const std::string &
    LogicalName() const;
    
    /**
     * @brief Opens next feature file
     */
    void
    OpenNext();
    
    /**
     * @brief Reads full feature matrix from a feature file
     * @param pMatrix matrix to be created and filled with read data
     * @return number of successfully read feature vectors
     */
    int
    ReadFullMatrix(Matrix<FLOAT>* pMatrix, const std::string & rFeatureName);    
    
    /**
     * @brief Reads feature vectors from a feature file
     * @param pMatrix matrix to be (only!) filled with read data. 
     * @return number of successfully read feature vectors
     * 
     * The function tries to fill @c pMatrix with feature vectors comming from
     * the current stream. If there are less vectors left in the stream, 
     * they are used and true number of successfuly read vectors is returned.
     */
    int
    ReadPartialMatrix(Matrix<FLOAT>* pMatrix);    
    
    /**
     * @brief Returns true if there are no feature files left on input
     */
    const bool
    EndOfList() const { return mInputQueue.empty(); }
    
    
  }; // class FeatureStream
}; //namespace STK

#endif // STK_Features_h
