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

#include "Models.h"
#include "labels.h"
#include "dict.h"
#include "net.h"

//#define DEBUG_MSGS
//#define TRACE_TOKENS

//typedef struct _Network Network;
//typedef struct _SubNet SubNet;
//typedef struct _Cache Cache;
//typedef struct _WLR WLR;
class Cache;
class WLR;
class SubNet;
class Network;

typedef enum 
{
  SP_ALL, 
  SP_TRUE_ONLY
} SearchPathsType;


class Cache 
{
public:
  FLOAT value;
  long  time;
};

class SubNet 
{
public:
  int   nnodes;
  Node  *first;
  Node  *last;
  Node  *nodes;
};

typedef enum 
{
  NO_ALIGNMENT    = 0,
  WORD_ALIGNMENT  = 1,
  MODEL_ALIGNMENT = 2,
  STATE_ALIGNMENT = 4,
  FRAME_ALIGNMENT = 8
} Alignment;

typedef enum 
{
  AT_ML=0,
  AT_MPE,
  AT_MFE,
  AT_MCE
} AccumType;

class Network 
{
public:
//Subnet part
//  int   nnodes;
  Node  *       first;
  Node  *       last;
//  Node  *nodes;

  Node  *       activeModels;
  Node  *       activeNodes;
  int           activeTokens;

  int           nNetStates;
  Token *       auxTokens;
  Cache *       outPCache;
  Cache *       mixPCache;

  long          mTime;
  Token *       bestToken;
  Node  *       bestNode;
  struct        my_hsearch_data wordHash;
// Passing parameters
  State *       threshState;
  FLOAT         beamThresh;
  FLOAT         wordThresh;

  FLOAT         wPenalty;
  FLOAT         mPenalty;
  FLOAT         pronScale;
  FLOAT         lmScale;
  FLOAT         tranScale;
  FLOAT         outpScale;
  FLOAT         ocpScale;
  FLOAT         pruningThresh;
  SearchPathsType SearchPaths;

  FLOAT (*OutputProbability)(State *state, FLOAT *observation, Network *network);
  int (*PassTokenInNetwork)(Token *from, Token *to, FLOAT addLogLike);
  int (*PassTokenInModel)(Token *from, Token *to, FLOAT addLogLike);
  PropagDir propagDir;
  int alignment;
  int collectAlphaBeta;
//  int mmi_den_pass;
  AccumType accumType;
  ModelSet *hmmSet;
  ModelSet *hmmSetToUpdate;
};

#define IS_ACTIVE(token) ((token).like > LOG_MIN)

class Token 
{
public:
  double like;
  WLR   *wlr;
  FloatInLog accuracy;
#ifdef bordel_staff
  WLR   *twlr;
  FLOAT bestlike;
#endif
};

struct AlphaBeta 
{
  FLOAT alpha;
  FLOAT beta;
  FloatInLog alphaAccuracy;
  FloatInLog betaAccuracy;
};


class FWBWR 
{
public:
  FWBWR *     mpNext;
  int         mTime;
  AlphaBeta   mpState[1];
};


class WLR 
{
public:
  Node  *     node;
  int         state_idx;
  FLOAT       like;
  long        mTime;
  WLR   *     mpNext;
  int         refs;
#ifdef DEBUG_MSGS
  WLR   *     tmpNext;
  BOOL        freed;
#endif
};


void LoadRecognitionNetwork(
  char *netFileName,
  ModelSet *hmms,
  Network *net);

void ReadRecognitionNetwork(
  FILE *fp,
  ModelSet *hmms,
  Network *network,
  LabelFormat labelFormat,
  long sampPeriod,
  const char *file_name,
  const char *in_MLF);

void InitNetwork(Network *net, Node *first, ModelSet *hmms, ModelSet *hmmsToUptade);
void ReleaseNetwork(Network *network);
void ViterbiInit(Network *network);
void ViterbiStep(Network *network, FLOAT *observation);
FLOAT ViterbiDone(Network *network, Label **labels);
Label *GetLabels(Token *token);
void KillToken(Token *token);

FLOAT *StateOccupationProbability(
  Network *network,
  FLOAT *observationMx,
  ModelSet *hmmset,
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


#endif // VITERBI_H
