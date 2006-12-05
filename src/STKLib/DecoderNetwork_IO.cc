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

#include "DecoderNetwork.h"

#include <sstream>

#define INIT_NODE_HASH_SIZE 1000

// CODE
//
namespace STK
{
  
  //###########################################################################
  //###########################################################################
  // GENERAL FUNCTIONS
  //###########################################################################
  //###########################################################################
  
  
  //***************************************************************************
  //***************************************************************************
  DecoderNetwork::NodeType*
  find_or_create_node(struct MyHSearchData *node_hash, const char *node_id, DecoderNetwork::NodeType** last)
  // Auxiliary function used by ReadHTKLattice_new. (Optionally initialize
  // uninitialized node_hash) Search for node record at key node_id. If found,
  // pointer to this record is returned, otherwise new node record is allocated
  // with type set to NT_UNDEF and entered to has at key node_id and pointer
  // to this new record is returned.
  {
    DecoderNetwork::NodeType*   p_node;
    ENTRY   e = {0}; //{0} is just to make compiler happy
    ENTRY*  ep;
  
    if (node_hash->mTabSize == 0 && !my_hcreate_r(INIT_NODE_HASH_SIZE, node_hash))
      Error("Insufficient memory");
    
    e.key = const_cast<char *>(node_id);
    my_hsearch_r(e, FIND, &ep, node_hash);
  
    if (ep != NULL) 
      return (DecoderNetwork::NodeType*) ep->data;  
  
    p_node = (DecoderNetwork::NodeType*) calloc(1, sizeof(DecoderNetwork::NodeType));
    
    if (p_node == NULL) 
      Error("Insufficient memory");
  
    p_node->rNLinks()     = 0;
    p_node->rpLinks()     = NULL;
    p_node->rNBackLinks() = 0;
    p_node->rpBackLinks() = NULL;
    p_node->mType       = NT_WORD;
    p_node->SetStart(UNDEF_TIME);
    p_node->SetStop(UNDEF_TIME);
    p_node->mpPronun    = NULL;
    p_node->mpBackNext  = *last;
    p_node->mPhoneAccuracy = 1.0;
    
    *last  = p_node;
    e.key  = strdup(node_id);
    e.data = p_node;
  
    if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, node_hash))
      Error("Insufficient memory");
    
    return p_node;
  } // find_or_create_node(...)

  
  //###########################################################################
  //###########################################################################
  // OUTPUT FUNCTIONS
  //###########################################################################
  //###########################################################################
  
  //***************************************************************************
  //***************************************************************************
  int fprintBase62(FILE *fp, int v)
  {
    int i = 0;
    char str[16];
    char *tab="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    if (v == 0) {
      fputc(tab[0], fp);
      return 1;
    }
    while (v) {
      str[i++] = tab[v % 62];
      v /= 62;
    }
    v=i;
    while (i--) fputc(str[i], fp);
    return v;
  }

  
  
    
  //###########################################################################
  //###########################################################################
  // INPUT FUNCTIONS
  //###########################################################################
  //###########################################################################
  
  //***************************************************************************
  //***************************************************************************
  int 
  getInteger(char *str, char **endPtr, const char *file_name, int line_no)
  {
    long l = strtoul(str, endPtr, 10);
  
    if (str == *endPtr || (**endPtr && !isspace(**endPtr))) {
      Error("Invalid integral value (%s:%d)", file_name, line_no);
    }
    while (isspace(**endPtr)) ++*endPtr;
  
    return l;
  }
  
  //***************************************************************************
  //***************************************************************************
  float 
  getFloat(char *str, char **endPtr, const char *file_name, int line_no)
  {
    double d = strtod(str, endPtr);
  
    if (str == *endPtr || (**endPtr && !isspace(**endPtr))) {
      Error("Invalid float value (%s:%d)", file_name, line_no);
    }
    while (isspace(**endPtr)) ++*endPtr;
  
    return d;
  }
  
  //***************************************************************************
  //***************************************************************************
  static int 
  getNodeNumber(int nnodes, char *str, char **endPtr, const char *file_name,
      int line_no)
  {
    long node_id = getInteger(str, endPtr, file_name, line_no);
  
    if (node_id < 0 || node_id >= nnodes) {
      Error("Node number out of range (%s:%d)", file_name, line_no);
    }
    return node_id;
  }
    
  //***************************************************************************
  //***************************************************************************
  int 
  RemoveCommentLines(FILE *fp)
  {
    int lines = 0, ch;
  
    for (;;) {
      while (isspace(ch=fgetc(fp))) {
        if (ch == '\n') {
          lines ++;
        }
      }
      if (ch != '#') {
        ungetc(ch, fp);
        return lines;
      }
      lines++;
      while ((ch=fgetc(fp)) != '\n' && ch != EOF);
    }
  }
  
  
} // namespace STK

