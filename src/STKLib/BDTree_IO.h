#ifndef BDTree_IO_h
#define BDTree_IO_h

#include <iostream>

namespace STK
{
  class BDTreeHeader;

  /** 
   * @brief BDTree IO header encapsulation
   */
  struct BDTreeHeader
  {
    static const char* INTRO ;

    int mHeaderSize;    ///< header size
    int mFileVersion;   ///< file version
    int mOrder;         ///< n-gram order
    int mVocabSize;     ///< target vocabulary size
    int mPredictorVocabSize; // predictor vocabulary size
    int mBinary;        ///< true if binary, false if ASCII format
    int mExtra0;
    int mExtra1;

    BDTreeHeader()
    : mHeaderSize(sizeof(BDTreeHeader)), 
      mFileVersion(0), 
      mOrder(0),
      mVocabSize(0),
      mPredictorVocabSize(0),
      mBinary(1),
      mExtra0(0),
      mExtra1(0)
    { }

    /** 
     * @brief Reads header from open stream
     * @param rStream std::istream
     */
    void
    Read(std::istream& rStream);

    /** 
     * @brief Writes header to stream
     * @param rStream std::ostream
     */
    void
    Write(std::ostream& rStream);

    /** 
     * @brief Writes header to stream with additional Predictor vocab size parameter, used for LVCSR binary format
     * @param rStream std::ostream
     */
    void
    Write_bin1(std::ostream& rStream);
  }; //struct BDTreeHeader


} // namespace STK

#endif

