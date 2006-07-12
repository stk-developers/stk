#include "Models.h"
#include "Viterbi.h"

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

      void
      SetSpkModel(Hmm* pSpkModel)
      { mpSpkModel = pSpkModel; }


      /** ************************************************************************
       *  ************************************************************************ 
       * @brief Computes the feature matrix score and adds it to the accumulators
       *
       * @param rFMatrix Feature matrix to process
       */
      void
      ParseFeatureMatrix(const Matrix<FLOAT>& rFMatrix)
      {
        State*  ubm_state = mpUbm->mpState[0];
        State*  spk_state = mpSpkModel->mpState[0];
        FLOAT   ubm_score = 0.0;
        FLOAT   spk_score = 0.0;
        int     time      = - mpModelSet->mTotalDelay;
        size_t  i;


        mpModelSet->ResetXformInstances();

        // compute UBM score
        for (size_t t = 0; t < rFMatrix.Rows(); t++)
        {
          // per feature vector score accums
          FLOAT     tmp_ubm_score(LOG_0);
          FLOAT     tmp_spk_score(LOG_0);

          mpModelSet->UpdateStacks(rFMatrix[t], ++time, FORWARD);

          if (time <= 0) 
            continue;
          

          // compute the UBM frame-score (log-add all mixture scores)
          for (i = 0; i < ubm_state->mNMixtures; i++)
          {
            Mixture*  ubm_mix     = ubm_state->mpMixture[i].mpEstimates;
            FLOAT     ubm_weight  = ubm_state->mpMixture[i].mWeight;
            FLOAT*    ubm_l_obs   = XformPass(ubm_mix->mpInputXform, rFMatrix[t], time, FORWARD);
            assert(ubm_l_obs != NULL);
            FLOAT     ubm_glike   = DiagCGaussianDensity(ubm_mix, ubm_l_obs, NULL);
            tmp_ubm_score = LogAdd(tmp_ubm_score, ubm_glike + ubm_weight);
          }

          ubm_score += tmp_ubm_score;
          
          
          // compute SpkModel frame score
          for (i = 0; i < spk_state->mNMixtures; i++)
          {
            FLOAT     spk_weight  = spk_state->mpMixture[i].mWeight;
            Mixture*  spk_mix     = spk_state->mpMixture[i].mpEstimates; 
            FLOAT*    spk_l_obs   = XformPass(spk_mix->mpInputXform, rFMatrix[t], time, FORWARD);
            assert(spk_l_obs != NULL);
            FLOAT     spk_glike   = DiagCGaussianDensity(spk_mix, spk_l_obs, NULL);
            tmp_spk_score = LogAdd(tmp_spk_score, spk_glike + spk_weight);
          }

          spk_score += tmp_spk_score;
        } // for(t = 0; t < header.mNSamples; t++)  


        // update the accumulators
        FLOAT score = spk_score - ubm_score;
        mMAccum += score;
        mSAccum += score * score;
        mNAccum += 1; 
      }



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

