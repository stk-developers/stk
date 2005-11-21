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

#ifndef VITERBI_H
#define VITERBI_H

#define bordel_staff

#include "hmms.h"
#include "labels.h"
#include "dict.h"
#include "net.h"

//#define DEBUG_MSGS
//#define TRACE_TOKENS

typedef struct _Network Network;
typedef struct _SubNet SubNet;
typedef struct _Cache Cache;
typedef struct _WLR WLR;


struct _Cache {
  FLOAT value;
  long  time;
};

struct _SubNet {
  int   nnodes;
  Node  *first;
  Node  *last;
  Node  *nodes;
};

typedef enum {
  NO_ALIGNMENT    = 0,
  WORD_ALIGNMENT  = 1,
  MODEL_ALIGNMENT = 2,
  STATE_ALIGNMENT = 4,
  FRAME_ALIGNMENT = 8
} Alignment;

typedef enum {
  AT_ML=0,
  AT_MPE,
  AT_MFE,
  AT_MCE
} AccumType;

struct _Network {
//Subnet part
//  int   nnodes;
  Node  *first;
  Node  *last;
//  Node  *nodes;

  Node  *activeModels;
  Node  *activeNodes;
  int   activeTokens;

  int   nNetStates;
  Token *auxTokens;
  Cache *outPCache;
  Cache *mixPCache;

  long  time;
  Token *bestToken;
  Node  *bestNode;
  struct my_hsearch_data wordHash;
// Passing parameters
  State *threshState;
  FLOAT beamThresh;
  FLOAT wordThresh;

  FLOAT wPenalty;
  FLOAT mPenalty;
  FLOAT pronScale;
  FLOAT lmScale;
  FLOAT tranScale;
  FLOAT outpScale;
  FLOAT ocpScale;
  FLOAT pruningThresh;
  enum {SP_All, SP_TrueOnly} SearchPaths;

  FLOAT (*OutputProbability)(State *state, FLOAT *observation, Network *network);
  int (*PassTokenInNetwork)(Token *from, Token *to, FLOAT addLogLike);
  int (*PassTokenInModel)(Token *from, Token *to, FLOAT addLogLike);
  PropagDir propagDir;
  int alignment;
  int collectAlphaBeta;
//  int mmi_den_pass;
  AccumType accumType;
  HMMSet *hmmSet;
  HMMSet *hmmSetToUpdate;
};

#define IS_ACTIVE(token) ((token).like > LOG_MIN)

struct _Token {
  double like;
  WLR   *wlr;
  FloatInLog accuracy;
#ifdef bordel_staff
  WLR   *twlr;
  FLOAT bestlike;
#endif
};

struct AlphaBeta {
   FLOAT alpha;
   FLOAT beta;
   FloatInLog alphaAccuracy;
   FloatInLog betaAccuracy;
};

struct _FWBWR {
 struct _FWBWR *next;
 int time;
 struct AlphaBeta state[1];
};

struct _WLR {
  Node  *node;
  int   state_idx;
  FLOAT like;
  long  time;
  WLR   *next;
  int   refs;
#ifdef DEBUG_MSGS
  WLR   *tmpNext;
  BOOL  freed;
#endif
};

#ifdef __cplusplus
  extern "C" {
#endif

void LoadRecognitionNetwork(
  char *netFileName,
  HMMSet *hmms,
  Network *net);

void ReadRecognitionNetwork(
  FILE *fp,
  HMMSet *hmms,
  Network *network,
  LabelFormat labelFormat,
  long sampPeriod,
  const char *file_name,
  const char *in_MLF);

void InitNetwork(Network *net, Node *first, HMMSet *hmms, HMMSet *hmmsToUptade);
void ReleaseNetwork(Network *network);
void ViterbiInit(Network *network);
void ViterbiStep(Network *network, FLOAT *observation);
FLOAT ViterbiDone(Network *network, Label **labels);
Label *GetLabels(Token *token);
void KillToken(Token *token);

FLOAT *StateOccupationProbability(
  Network *network,
  FLOAT *observationMx,
  HMMSet *hmmset,
  int nFrames,
  FLOAT **outProbOrMahDist,
  int getMahalDist);
  
FLOAT MCEReest(
  Network *net,
  FLOAT *obsMx,
  FLOAT *obsMx2,
  int nFrames,
  FLOAT weight,
  FLOAT sigSlope);

FLOAT BaumWelchReest(
  Network *net,
  FLOAT *obsMx,
  FLOAT *obsMx2,
  int nFrames,
  FLOAT weight);

FLOAT ViterbiReest(
  Network *net,
  FLOAT *observationMx,
  FLOAT *observationMx2,
  int nFrames,
  FLOAT weight);

FLOAT DiagCGaussianMixtureDensity(State *state, FLOAT *obs, Network *network);
FLOAT FromObservationAtStateId(State *state, FLOAT *obs, Network *network);
int PassTokenMax(Token *from, Token *to, FLOAT addLogLike);
int PassTokenSum(Token *from, Token *to, FLOAT addLogLike);
int PassTokenSumUnlogLikes(Token *from, Token *to, FLOAT addLogLike);

WLR *TimePruning(Network *network, int frame_delay);

#ifdef __cplusplus
}
#endif

#endif // VITERBI_H
