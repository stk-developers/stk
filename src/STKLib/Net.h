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

#ifndef STK_Net_h
#define STK_Net_h

#include "Models.h"
#include "labels.h"
#include "dict.h"
#include "common.h"

namespace STK
{
  class Node;
  class Link;
  class Token;
  class FWBWR;
  
  
  enum NodeType 
  {
    NT_Undef  = 0x00,
    NT        = 0x01,
    NT_Model  = 0x02,
    NT_Phone  = 0x04,
    NT_Subnet = 0x08,
    NT_Tee    = 0x10,
    NT_Sticky = 0x20,
    NT_True   = 0x40
  };
  
  
  class Link 
  {
  public:
    Node *      mpNode;
    FLOAT       mLike;
  };
  
  
  class Node
  {
  public:
  //  union {
      char   *mpName;
  //    struct {
        Hmm  *hmm;
        Hmm  *hmmToUpdate;
  //    };
      Pronun *pronun;
  //    SubNet *subnet;
  //  };
  
    int           aux;
    int           mType;
    Node  *       mpNext;
    Node  *       backnext;
    int           nlinks;
    int           nbacklinks;
    Link  *       links;
    Link  *       backlinks;
  
    //time range when model can be active - apply only for model type
    long long start;
    long long stop;
    FLOAT phoneAccuracy;
  
  #ifndef EXPANDNET_ONLY
    Token *tokens;
    Token *exitToken;
  
    //id of first emiting state - apply only for model type
    int   estate_id;
  
    FWBWR *alphaBetaList;
    FWBWR *alphaBetaListReverse;
  
    Node  *nextActiveModel;
    Node  *prevActiveModel;
  
    Node  *nextActiveNode;
    Node  *prevActiveNode;
    int   isActive;
    int   isActiveNode;
  #ifndef NDEBUG
    int   aux2;
  #endif
  #endif
  };
  
  //typedef struct _STKNetworkOutputFormat STKNetworkOutputFormat;
  
  class STKNetworkOutputFormat 
  {
  public:
    unsigned no_LM_likes    : 1;
    unsigned no_times       : 1;
    unsigned start_times    : 1;
    unsigned no_word_nodes  : 1;
    unsigned no_model_nodes : 1;
    unsigned no_pronun_vars : 1;
    unsigned no_defaults    : 1;
    unsigned all_field_names: 1;
    unsigned arc_defs_to_end: 1;
    unsigned arc_defs_with_J: 1;
    unsigned base62_labels  : 1;
    unsigned aprox_accuracy : 1;
  
    //Have no effect yet
    unsigned no_acc_likes   : 1;
    unsigned strip_triphones: 1;
    unsigned lin_node_seqs  : 1;
  };
  
  enum NotInDictAction 
  {
    WORD_NOT_IN_DIC_UNSET = 0,
    WORD_NOT_IN_DIC_ERROR = 1,
    WORD_NOT_IN_DIC_WARN  = 2,
    PRON_NOT_IN_DIC_ERROR = 4
  };
  
  //typedef struct _ExpansionOptions ExpansionOptions;
  class ExpansionOptions 
  {
  public:
    unsigned no_optimization    : 1;
    unsigned no_word_expansion  : 1;
    unsigned respect_pronun_var : 1;
    unsigned remove_words_nodes : 1;
    unsigned CD_phone_expansion : 1;
    unsigned strict_timing      : 1;
    unsigned trace_flag;
  };
  
  
  Node *MakeNetworkFromLabels(Label *labels, enum NodeType node_type);
  
  void ExpandWordNetworkByDictionary(
    Node *first,
    struct my_hsearch_data *dict,
    int keep_word_nodes,
    int multiple_pronun);
  
  void ExpandMonophoneNetworkToTriphones(
    Node *first,
    struct my_hsearch_data *nonCDphones,
    struct my_hsearch_data *CDphones);
  
  void LatticeLocalOptimization(
    Node *first,
    int strictTiming,
    int trace_flag);
  
  Node *DiscardUnwantedInfoInNetwork(
    Node *first,
    STKNetworkOutputFormat format);
  
  void WriteSTKNetwork(
    FILE *flp,
    Node *node,
    STKNetworkOutputFormat format,
    long sampPeriod,
    const char *label_file,
    const char *out_MNF);
  
  void WriteSTKNetworkInOldFormat(
    FILE *flp,
    Node *node,
    LabelFormat labelFormat,
    long sampPeriod,
    const char *label_file,
    const char *out_MNF);
  
  void FreeNetwork(Node *node);
  
  Node *ReadSTKNetwork(
    FILE *lfp,
    struct my_hsearch_data *word_hash,
    struct my_hsearch_data *phone_hash,
    int notInDict,
    LabelFormat labelFormat,
    long sampPeriod,
    const char *file_name,
    const char *in_MLF);
  
  Node *ReadSTKNetworkInOldFormat(
    FILE *lfp,
    struct my_hsearch_data *word_hash,
    struct my_hsearch_data *phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char *file_name,
    const char *in_MLF);
  
  Node *ReadHTKLattice(
    FILE *lfp,
    struct my_hsearch_data *word_hash,
    struct my_hsearch_data *phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char *file_name);
  
  void ComputeAproximatePhoneAccuracy(
    Node *first,
    int type);
  
  void SelfLinksToNullNodes(Node *first);
  int RemoveRedundantNullNodes(Node *first);
  
  void NetworkExpansionsAndOptimizations(
    Node *node,
    ExpansionOptions expOptions,
    STKNetworkOutputFormat out_net_fmt,
    struct my_hsearch_data *dictHash,
    struct my_hsearch_data *nonCDphHash,
    struct my_hsearch_data *triphHash);

}; // namespace STK

#endif // STK_Net_h
