
#include "FeaCatAdapter.h"

#include "common.h"
#include "Models.h"

namespace STK {
 
  
  FeaCatAdapter::
  FeaCatAdapter() 
    : mpModelSet(NULL), mpXformInstance(NULL), mNframes(0), mTime(0)
  { }


  FeaCatAdapter::
  ~FeaCatAdapter()
  { delete mpModelSet; delete mpXformInstance; }


  void
  FeaCatAdapter::
  Init(const char* mmfFilename)
  {
    mpModelSet = new ModelSet;
    mpModelSet->Init(MODEL_SET_WITH_ACCUM);
    mpModelSet->mUseNewMatrix = true;
    mpModelSet->ParseMmf(mmfFilename, NULL);
    mpModelSet->ExpandPredefXforms();

    if(NULL != mpModelSet->mpInputXform) {
      mpXformInstance = mpModelSet->mpInputXform;
    } else {
      Error((std::string("No InputXform in the file:")+mmfFilename).c_str());
    }
  }


  void
  FeaCatAdapter::
  XReset()
  {
    mpModelSet->ResetXformInstances();
    mNframes = 0; 
    mTime = -mpXformInstance->mTotalDelay;
  }


  const float* 
  FeaCatAdapter::
  XPass(float* in_vec)
  {
#if DOUBLEPRECISION==0
    const float* out;
    mTime++;
    mNframes++;
    mpModelSet->UpdateStacks(in_vec, mTime, FORWARD);

    //It's terribly slow, why??? 
    out = XformPass(mpXformInstance, in_vec, mTime, FORWARD);

    if(mTime <= 0) { return NULL; }

    return out;
#else
    Error("FeaCatAdapter in double precision not supported!!!");
    return NULL;
#endif
  }


  size_t
  FeaCatAdapter::
  InputDim()
  {
    return 0;
  }


  size_t
  FeaCatAdapter::
  OutputDim()
  {
    return mpXformInstance->OutSize();
  }


  size_t
  FeaCatAdapter::
  Delay()
  {
    return mpXformInstance->mTotalDelay;
  }


}

