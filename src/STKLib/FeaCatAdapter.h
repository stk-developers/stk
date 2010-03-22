#ifndef _TNET_FEA_CAT_ADAPTER_H_
#define _TNET_FEA_CAT_ADAPTER_H_

#include <cstring>

namespace STK {
  
  class ModelSet;
  class XformInstance;

  class FeaCatAdapter {
    public:
      FeaCatAdapter();
      ~FeaCatAdapter();

      void Init(const char* mmfFilename);
      void XReset();
      const float* XPass(float* in_vec);

      size_t InputDim();
      size_t OutputDim();
      size_t Delay();

      size_t Nframes()
      { return mNframes; }
      int Time()
      { return mTime; }

    private:
      ModelSet* mpModelSet;
      XformInstance* mpXformInstance;

      size_t mNframes;
      int mTime;

  };

}


#endif

