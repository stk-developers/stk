#include "Models.h"
#include "Decoder.h"

namespace STK
{

  class ZNorm
  {
    private:
      ModelSet* mpModelSet;
      Hmm*      mpUbm;
      Hmm*      mpSpkModel;

      FLOAT   mUbmScore;   ///< UBM score
      FLOAT   mNAccum;     ///< Accumlates occurances
      FLOAT   mMAccum;     ///< First order statistics
      FLOAT   mSAccum;     ///< Second order statistics


    public:
      /**
       * @brief The constructor
       */
      ZNorm() : mpModelSet(NULL), mpUbm(NULL), mpSpkModel(NULL)
      { Init(); }

      /**
       * @brief The destructor
       */
      ~ZNorm() 
      {}
      

      /**
       * @brief Resets the accumulators for new Z-Norm computation
       */
      void
      Init()
      { 
        mUbmScore = 0.0;
        mNAccum   = 0.0;
        mMAccum   = 0.0;
        mSAccum   = 0.0;
      }


      void
      SetModelSet(ModelSet* pModelSet)
      { mpModelSet = pModelSet; }
      
      
      void
      SetUbm(Hmm* pUbm)
      { mpUbm = pUbm; }

      
      Hmm*
      pUbm()
      { return mpUbm; }


      void
      SetSpkModel(Hmm* pSpkModel)
      { mpSpkModel = pSpkModel; }

      
      Hmm*
      pSpkModel()
      { return mpSpkModel; }


      /** ************************************************************************
       *  ************************************************************************ 
       * @brief Computes the feature matrix score and adds it to the accumulators
       *
       * @param rFMatrix Feature matrix to process
       */
      void
      ParseFeatureMatrix(const Matrix<FLOAT>& rFMatrix);


      /** ************************************************************************
       *  ************************************************************************ 
       * @brief Gives access to the mUbmScore
       *
       * @return Returns a refference to the mUbmScore variable
       */
      FLOAT&
      rUbm()
      {
        return mUbmScore;
      }
      

      /** ************************************************************************
       *  ************************************************************************ 
       * @brief Adds a score to the accumulators
       *
       * @param score Score that will be added to the accumulators
       */
      void 
      AddScore(FLOAT score) 
      {
        mMAccum += score;
        mSAccum += score * score;
        mNAccum += 1; 
      }


      /** ************************************************************************
       *  ************************************************************************ 
       * @brief Finalizes the computation of Z-Norm mean and variance
       */
      void
      Finalize()
      {
        mMAccum = mMAccum / mNAccum;
        mSAccum = sqrt(mSAccum / mNAccum - mMAccum * mMAccum);
      }


      /** ************************************************************************
       *  ************************************************************************ 
       * @brief Returns the Z-Norm standard deviation
       *
       * @return The Z-Norm standard deviation
       *
       * It has to be provided that Finalize() has been called EXACTLY once
       * otherwise second order statisctic accumulator is returned
       */
      const FLOAT&
      StdDev() const
      {
        return mSAccum;
      }


      /** ************************************************************************
       *  ************************************************************************ 
       * @brief Returns the Z-Norm mean
       *
       * @return The Z-Norm mean
       *
       * It has to be provided that Finalize() has been called EXACTLY once
       * otherwise first order statisctic accumulator is returned
       */
      const FLOAT&
      Mean() const
      {
        return mMAccum;
      }
  };
}; // namespace STK

