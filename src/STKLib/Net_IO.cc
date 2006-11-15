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


// PROJECT INCLUDES
//
#include "common.h"
#include "labels.h"
#include "Net.h"


#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <malloc.h>
#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <sstream>
#include <math.h>

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
  Node<NODE_REGULAR, LINK_REGULAR>*
  find_or_create_node(struct MyHSearchData *node_hash, const char *node_id, Node<NODE_REGULAR, LINK_REGULAR> **last)
  // Auxiliary function used by ReadHTKLattice_new. (Optionally initialize
  // uninitialized node_hash) Search for node record at key node_id. If found,
  // pointer to this record is returned, otherwise new node record is allocated
  // with type set to NT_UNDEF and entered to has at key node_id and pointer
  // to this new record is returned.
  {
    Node<NODE_REGULAR, LINK_REGULAR>*   p_node;
    ENTRY   e = {0}; //{0} is just to make compiler happy
    ENTRY*  ep;
  
    if (node_hash->mTabSize == 0 && !my_hcreate_r(INIT_NODE_HASH_SIZE, node_hash))
      Error("Insufficient memory");
    
    e.key = const_cast<char *>(node_id);
    my_hsearch_r(e, FIND, &ep, node_hash);
  
    if (ep != NULL) 
      return (Node<NODE_REGULAR, LINK_REGULAR> *) ep->data;  
  
    p_node = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
    
    if (p_node == NULL) 
      Error("Insufficient memory");
  
    p_node->mNLinks     = 0;
    p_node->mpLinks     = NULL;
    p_node->mNBackLinks = 0;
    p_node->mpBackLinks = NULL;
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

  



  //***************************************************************************
  //***************************************************************************
  /*
  void 
  WriteSTKNetwork(
    FILE*                   pFp,
    Node<NODE_REGULAR, LINK_REGULAR>*                   pFirst,
    STKNetworkOutputFormat  format,
    long                    sampPeriod,
    const char*             net_file,
    const char*             out_MNF)
  {
    int     n;
    int     l=0;
    Node<NODE_REGULAR, LINK_REGULAR>*   node;
  
    for (n = 0, node = pFirst; node != NULL; node = node->mpNext, n++)  
    {
      node->mAux = n;
      l += node->mNLinks;
    }
    
    fprintf(pFp,"N=%d L=%d\n", n, l);
    for (node = pFirst; node != NULL; node = node->mpNext)
    {
      int j;
  
      if (format.mAllFieldNames) fputs("I=", pFp);
      if (format.mBase62Labels) fprintBase62(pFp, node->mAux);
      else                      fprintf(pFp,"%d", node->mAux);
  
      if (!format.mNoTimes && node->Stop() != UNDEF_TIME) {
        fputs(" t=", pFp);
  
        if (node->Start() != UNDEF_TIME && format.mStartTimes) {
          fprintf(pFp,"%g,", node->Start() * 1.0e-7 * sampPeriod);
        }
        fprintf(  pFp,"%g",  node->Stop()  * 1.0e-7 * sampPeriod);
      }

      if (!(node->mType & NT_WORD && node->mpPronun == NULL)
        || !format.mNoDefaults) 
      {
        putc(' ', pFp);
        putc(node->mType & NT_WORD   ? 'W' :
             node->mType & NT_SUBNET ? 'S' :
                                       'M', pFp); // NT_MODEL, NT_PHONE
        putc('=', pFp);
        fprintHTKstr(pFp, node->mType & NT_MODEL   ? node->mpHmm->mpMacro->mpName   :
                          node->mType & NT_WORD    ? (!node->mpPronun ? "!NULL" :
                                                    node->mpPronun->mpWord->mpName) :
                                                    node->mpName); // NT_PHONE (NT_SUBNET)
      }

      if (!format.mNoPronunVars && node->mType & NT_WORD
      && node->mpPronun != NULL && node->mpPronun->mpWord->npronuns > 1
      && (node->mpPronun->variant_no > 1 || !format.mNoDefaults))
      {
        fprintf(pFp," v=%d", node->mpPronun->variant_no);
      }

      if (node->mType & NT_TRUE || node->mType & NT_STICKY) 
      {
        fputs(" f=", pFp);
        if (node->mType & NT_TRUE)   putc('T', pFp);
        if (node->mType & NT_STICKY) putc('K', pFp);
      }

      if (node->mType & NT_PHONE && node->mPhoneAccuracy != 1.0) {
        fprintf(pFp," p="FLOAT_FMT, node->mPhoneAccuracy);
      }

      if (!format.mArcDefsToEnd) 
      {
        if (format.mAllFieldNames) fprintf(pFp," J=%d", node->mNLinks);
  
        for (j = 0; j < node->mNLinks; j ++) 
        {
          putc(' ', pFp);
          if (format.mAllFieldNames) fputs("E=", pFp);
          if (format.mBase62Labels) fprintBase62(pFp, node->mpLinks[j].mpNode->mAux);
          else                     fprintf(pFp,"%d", node->mpLinks[j].mpNode->mAux);

          // compute the language likelihood
          FLOAT aux_like = node->mpLinks[j].mAcousticLike;
          FLOAT l_like   = node->mpLinks[j].mLmLike;
          //double eps(l_like * EPSILON);

          //printf("eps = %f \n", eps);

          //if (_ABS(eps) <= EPSILON) l_like = 0.0;
          //if (_ABS(eps) <= EPSILON + EPSILON) l_like = 0.0;
           
          if ((!close_enough(node->mpLinks[j].mLmLike, 0.0, 10)) 
          &&  (!format.mNoLMLikes))
          {
            fprintf(pFp," l="FLOAT_FMT, node->mpLinks[j].mLmLike);
          }

          if ((node->mpLinks[j].mAcousticLike != 0.0) 
          && !(format.mNoAcousticLikes))
          {
            fprintf(pFp," a="FLOAT_FMT, node->mpLinks[j].mAcousticLike);
          }
        }
      }

      fputs("\n", pFp);

      if (ferror(pFp)) {
        Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
      }
    }
  
    if (format.mArcDefsToEnd) 
    {
      l = 0;
      for (node = pFirst; node != NULL; node = node->mpNext) 
      {
        int j;
  
        for (j = 0; j < node->mNLinks; j ++) 
        {
          if (format.mAllFieldNames) 
            fprintf(pFp, format.mArcDefsWithJ ? "J=%d S=" : "I=", l++);

          if (format.mBase62Labels) fprintBase62(pFp, node->mAux);
          else                     fprintf(pFp,"%d", node->mAux);
          putc(' ', pFp); // space = ' ';
          if (format.mAllFieldNames) fputs("E=", pFp);

          if (format.mBase62Labels) fprintBase62(pFp, node->mpLinks[j].mpNode->mAux);
          else                      fprintf(pFp,"%d", node->mpLinks[j].mpNode->mAux);

          // output language probability
          FLOAT aux_like = node->mpLinks[j].mAcousticLike;
          FLOAT l_like = node->mpLinks[j].mLmLike - aux_like;
          //double eps(l_like * EPSILON);

          //printf("eps = %f \n", eps);

          //if (_ABS(eps) <= EPSILON) l_like = 0.0;
          //if (_ABS(eps) <= EPSILON + EPSILON) l_like = 0.0;
           
          if ((!close_enough(node->mpLinks[j].mLmLike, 0.0, 10)) 
          && (!format.mNoLMLikes))
          {
            fprintf(pFp," l="FLOAT_FMT, node->mpLinks[j].mLmLike);
          }

          // output acoustic probability
          if ((node->mpLinks[j].mAcousticLike != 0.0) 
          && !(format.mNoAcousticLikes))
          {
            fprintf(pFp," a="FLOAT_FMT, node->mpLinks[j].mAcousticLike);
          }

          fputs("\n", pFp);

          if (ferror(pFp)) 
            Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
          
        }
      }
    }
  }
  */
  
  //***************************************************************************
  //***************************************************************************
  void WriteSTKNetworkInOldFormat(
    FILE        *lfp,
    Node<NODE_REGULAR, LINK_REGULAR>       *first,
    LabelFormat labelFormat,
    long        sampPeriod,
    const char  *net_file,
    const char  *out_MNF)
  {
    int   i;
    Node<NODE_REGULAR, LINK_REGULAR>* node;
  
    for (i = 0, node = first; node != NULL; node = node->mpNext, i++)  {
      node->mAux = i;
    }
    fprintf(lfp,"NUMNODES: %d\n", i);
    for (i = 0, node = first; node != NULL; node = node->mpNext, i++) {
      int j,
      type = node->mType & NT_MODEL       ? 'M'  :
            node->mType & NT_PHONE       ? 'M'  :
            node->mType & NT_SUBNET      ? 'S'  :
            node->mType & NT_WORD        ?
              (node->mpPronun == NULL     ?
                (node->mType & NT_STICKY ? 'F'  :
                                          'N') :
                (node->mType & NT_STICKY ? 'K'  :
                                          'W')):
                                          '?';
      if (!(node->mType & NT_TRUE)) {
        type = tolower(type);
      }
      fprintf(lfp,"%d\t%c %s",
              i, type,
              node->mType & NT_MODEL   ? node->mpHmm->mpMacro->mpName :
              node->mType & NT_PHONE   ? node->mpName :
              node->mType & NT_SUBNET  ? node->mpName :
              node->mType & NT_WORD    ?
                (node->mpPronun == NULL ? "-" :
                                        node->mpPronun->mpWord->mpName):
                                        "?");
      if (node->mType & NT_WORD && node->mpPronun) {
        if (node->mpPronun->mpWord->mpName != node->mpPronun->outSymbol) {
          fprintf(lfp," [%s]", node->mpPronun->outSymbol);
        }
        if (node->mpPronun->prob != 0.0 || node->mpPronun->mpWord->npronuns > 1) {
          fprintf(lfp," {%d "FLOAT_FMT"}",
                  node->mpPronun->variant_no,
                  node->mpPronun->prob);
        }
      }
      if (node->mType & NT_PHONE && node->mPhoneAccuracy != 1.0) {
        fprintf(lfp," {"FLOAT_FMT"}", node->mPhoneAccuracy);
      }
      if (!(labelFormat.TIMES_OFF) &&
        node->Start() != UNDEF_TIME && node->Stop() != UNDEF_TIME) {
        int ctm = labelFormat.CENTRE_TM;
        fprintf   (lfp," (");
        fprintf_ll(lfp, sampPeriod * (2 * node->Start() + ctm) / 2 - labelFormat.left_extent);
        fprintf   (lfp," ");
        fprintf_ll(lfp, sampPeriod * (2 * node->Stop() - ctm)  / 2 + labelFormat.right_extent);
        fprintf   (lfp,")");
      }
      fprintf(lfp,"\t%d", node->mNLinks);
      for (j = 0; j < node->mNLinks; j ++) {
        fprintf(lfp," %d",node->mpLinks[j].pNode()->mAux);
        if (node->mpLinks[j].LmLike() != 0.0) {
          fprintf(lfp," {"FLOAT_FMT"}", node->mpLinks[j].LmLike());
        }
      }
      fputs("\n", lfp);
      if (ferror(lfp)) {
        Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
      }
    }
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
  
  //***************************************************************************
  //***************************************************************************
  Node<NODE_REGULAR, LINK_REGULAR> *ReadSTKNetworkInOldFormat(
    FILE*                     lfp,
    struct MyHSearchData*     word_hash,
    struct MyHSearchData*     phone_hash,
    LabelFormat               labelFormat,
    long                      sampPeriod,
    const char*               file_name,
    const char*               in_MLF)
  {
    size_t    numOfNodes;
    size_t    i;
    size_t    j;
    size_t    ch;
    
    int       nodeId = 0;
    int       linkId;
    int       numOfLinks;
    int       pronunVar;
    long long start;
    long long stop;
    char      nodeType;
    char      wordOrModelName[1024] = {'\0'};
    double    linkLike;
    double    pronunProb;
    Node<NODE_REGULAR, LINK_REGULAR>*     node;
    Node<NODE_REGULAR, LINK_REGULAR>**    nodes;
  
    std::stringstream ss; 
    
    RemoveCommentLines(lfp);
  
    if (fscanf(lfp," %1023[^0-9]", wordOrModelName) == 1) {
      for ( i =0; i < strlen(wordOrModelName); i++) {
        wordOrModelName[i] = toupper(wordOrModelName[i]);
      }
      while (--i>=0 && (wordOrModelName[i] == '='||isspace(wordOrModelName[i]))) {
        wordOrModelName[i] = '\0';
      }
    }

    int t1, t2;    
    if ((strcmp(wordOrModelName, "NUMNODES:") &&
        strcmp(wordOrModelName, "NUMBEROFNODES")) ||        //Obsolete NumerOfNodes
      fscanf(lfp," %d NumberOfArcs=%d", &t1, &t2)<1){       //Obsolete NumerOfArcs
      Error("Syntax error in file %s\nKeyword NumNodes: is missing", file_name);
    }
    numOfNodes = t1;
    
    if ((nodes = (Node<NODE_REGULAR, LINK_REGULAR> **) calloc(numOfNodes, sizeof(Node<NODE_REGULAR, LINK_REGULAR> *))) == NULL) {
      Error("Insufficient memory");
    }
    for (i=0; i < numOfNodes; i++) {
      if ((nodes[i] = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL) {
        Error("Insufficient memory");
      }
      nodes[i]->mType = NT_UNDEF;
    }
    for (i=0; i < numOfNodes; i++) {
      RemoveCommentLines(lfp);
  
      switch (fscanf(lfp, "%d %c %1023s", &nodeId, &nodeType, wordOrModelName)) 
      {
        case  3:
          break; //OK
  
        case -1:
          for (j=0; j < numOfNodes; j++) 
          {
            if (nodes[j] == NULL) {
              Error("Node %d is not defined in file %s", (int) j, file_name);
            }
          }
          
        default:
          Error("Invalid syntax in definition of node %d in file %s",
                nodeId, file_name);
      }
      
      if (static_cast<size_t>(nodeId) >= numOfNodes) {
        Error("Invalid definition of node %d in file %s.\n"
              "Node Id is bigger than number of nodes", nodeId, file_name);
      }
      node = nodes[nodeId];
  
      if (node->mType != NT_UNDEF)
        Error("Redefinition of node %d in file %s", nodeId, file_name);
      
      if (toupper(nodeType) != nodeType) {
        nodeType = toupper(nodeType);
        node->mType = 0;
      } else {
        node->mType = NT_TRUE;
      }
      if (nodeType != 'M' && nodeType != 'W' && nodeType != 'N' &&
        nodeType != 'S' && nodeType != 'K' && nodeType != 'F') {
        Error("Invalid definition of node %d in file %s.\n"
              "Supported values for node type are: M - model, W - word, N - null, S - subnet, K - keyword, F - filler",
              nodeId, file_name);
      }
      node->SetStart(UNDEF_TIME); 
      node->SetStop (UNDEF_TIME);
      node->mPhoneAccuracy = 1.0;
  //    node->mAux = *totalNumOfNodes;
  //    ++*totalNumOfNodes;
  
      if (nodeType == 'S') {
        FILE *snfp;
        Node<NODE_REGULAR, LINK_REGULAR> *subnetFirst;
  
  //      --*totalNumOfNodes; // Subnet node doesn't count
  
        if ((snfp = fopen(wordOrModelName, "rt")) == NULL) {
          Error("Cannot open network file: %s", wordOrModelName);
        }
        subnetFirst = ReadSTKNetworkInOldFormat(
                        snfp, word_hash, phone_hash, labelFormat,
                        sampPeriod, wordOrModelName, NULL);
        subnetFirst->mNBackLinks = node->mNBackLinks;
        *node = *subnetFirst;
        free(subnetFirst);
        fclose(snfp);
      } else if (nodeType == 'M') {
        ENTRY e, *ep;
        node->mType |= NT_PHONE;
        e.key  = wordOrModelName;
        e.data = NULL;
        my_hsearch_r(e, FIND, &ep, phone_hash);
  
        if (ep == NULL) {
          e.key  = strdup(wordOrModelName);
          e.data = e.key;
  
          if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, phone_hash)) {
            Error("Insufficient memory");
          }
          ep->data = e.data;
        }
        node->mpName = (char *) ep->data;
        fscanf(lfp, " {%lf}", &pronunProb);
        // We are not interested in PhoneAccuracy
      } else {
        node->mType |= NT_WORD;
  
        if (nodeType == 'K' || nodeType == 'F') {
        node->mType |= NT_STICKY;
        nodeType = nodeType == 'K' ? 'W' : 'N';
        }
        if (nodeType == 'W') {
          ENTRY e = {0}; //{0} is just to make compiler happy
          ENTRY *ep;
          Word *word = NULL;
  
          e.key  = wordOrModelName;
          my_hsearch_r(e, FIND, &ep, word_hash);
  
          if (ep == NULL) {
            Error("Invalid definition of node %d in file %s.\n"
                  "Unknown word '%s'", nodeId, file_name, e.key);
          }
          word = (Word *) ep->data;
  
          while (isspace(ch = fgetc(lfp)));
          if (ch != '[') {
            ungetc(ch, lfp);
          } else {
            if (fscanf(lfp, "%1023[^]]", wordOrModelName) != 1) {
  //           pronun.outSymbol = NULL;
            } else{
  //           pronun.outSymbol = wordOrModelName;
            }
            if (fgetc(lfp) != ']'){
              Error("Invalid definition of node %d in file %s.\n"
                    "Missing ']' after output symbol definition", nodeId, file_name);
            }
          }
          if (fscanf(lfp, "{%d %lf}", &pronunVar, &pronunProb) != 2) {
            pronunProb = 0.0;
            pronunVar = 0;
          } else {
            pronunVar--;
          }
          if (word->npronuns <= pronunVar) {
            Error("Invalid definition of node %d in file %s.\n"
                  "Word %s does not have pronunciation varian %d",
                  nodeId, file_name, word->mpName, pronunVar+1);
          }
          node->mpPronun = word ? word->pronuns[pronunVar] : NULL;
        } else {
          node->mpPronun = NULL;
        }
      }
      if (nodeType != 'S') {
        if (fscanf(lfp, " (%lld %lld)", &start, &stop)==2 && !(labelFormat.TIMES_OFF)) {
          long center_shift = labelFormat.CENTRE_TM ? sampPeriod / 2 : 0;
          node->SetStart((start - center_shift - labelFormat.left_extent)  / sampPeriod);
          node->SetStop((stop  + center_shift + labelFormat.right_extent) / sampPeriod);
        }
      }
      if (fscanf(lfp, "%d ", &numOfLinks) != 1) {
        Error("Invalid definition of node %d in file %s.\n"
              "Number of links is expected", nodeId, file_name);
      }
      if (nodeType == 'S') { // Add links to the final node of the subnetwork
        while (node->mpNext != NULL) node = node->mpNext;
      }
      if (numOfLinks) {
        if ((node->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(numOfLinks * sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) {
          Error("Insufficient memory");
        }
      } else {
        if (nodeType == 'M') {
          Error("Invalid definition of node %d in file %s.\n"
                "Model node must have at least one link", nodeId, file_name);
        }
        node->mpLinks = NULL;
      }
      node->mNLinks = numOfLinks;
  
      for (j=0; j < static_cast<size_t>(numOfLinks); j++) {
        if (fscanf(lfp, "%d ", &linkId) != 1) {
          Error("Invalid definition of node %d in file %s.\n"
                "Link<NODE_REGULAR, LINK_REGULAR> Id is expected in link list", nodeId, file_name);
        }
        if (static_cast<size_t>(linkId) >= numOfNodes) {
          Error("Invalid definition of node %d in file %s.\n"
                "Link Id is bigger than number of nodes", nodeId, file_name);
        }
        if (fscanf(lfp, "{%lf} ", &linkLike) != 1) {
          linkLike = 0.0;
        }
        node->mpLinks[j].SetNode(nodes[linkId]);
        node->mpLinks[j].SetLmLike(linkLike);
        nodes[linkId]->mNBackLinks++;
      }
    }
    for (i = 1; i < numOfNodes-1; i++) {
      if (nodes[i]->mNLinks == 0) {
        if (nodes[numOfNodes-1]->mNLinks == 0) {
          Error("Network contains multiple nodes with no successors (%s)",
                file_name);
        }
        node = nodes[numOfNodes-1];
        nodes[numOfNodes-1] = nodes[i];
        nodes[i] = node;
      }
      if (nodes[i]->mNBackLinks == 0) {
        if (nodes[0]->mNBackLinks == 0) {
          Error("Network contains multiple nodes with no predecessor (%s)",
                file_name);
        }
        node = nodes[0];
        nodes[0] = nodes[i];
        nodes[i] = node;
        i--;
        continue; // Check this node again. Could be the first one
      }
    }
    if (nodes[0]->mNBackLinks != 0 || nodes[numOfNodes-1]->mNLinks != 0) {
      Error("Network contain no start node or no final node (%s)", file_name);
    }
    if (!(nodes[0]           ->mType & NT_WORD) || nodes[0]           ->mpPronun != NULL ||
      !(nodes[numOfNodes-1]->mType & NT_WORD) || nodes[numOfNodes-1]->mpPronun != NULL) {
      Error("Start node and final node must be Null nodes (%s)", file_name);
    }
    for (i = 0; i < numOfNodes-1; i++) {
      nodes[i]->mpNext = nodes[i+1];
    }
  
    // create back links
    for (i = 0; i < numOfNodes; i++) {
      if (!nodes[i]->mpBackLinks) // Could be allready alocated for subnetwork
        nodes[i]->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(nodes[i]->mNBackLinks * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
      if (nodes[i]->mpBackLinks == NULL) Error("Insufficient memory");
      nodes[i]->mNBackLinks = 0;
    }
    for (i = 0; i < numOfNodes; i++) {
      for (j=0; j < static_cast<size_t>(nodes[i]->mNLinks); j++) {
        Node<NODE_REGULAR, LINK_REGULAR> *forwNode = nodes[i]->mpLinks[j].pNode();

        forwNode->mpBackLinks[forwNode->mNBackLinks].SetNode(nodes[i]);
        forwNode->mpBackLinks[forwNode->mNBackLinks].SetLmLike(nodes[i]->mpLinks[j].LmLike());
        forwNode->mNBackLinks++;
      }
    }
    node = nodes[0];
    free(nodes);
  
    if (in_MLF) {
      char *chptr;
      do {
          char line[1024];
          if (fgets(line, sizeof(line), lfp) == NULL) {
            Error("Missing '.' at the end of network '%s' in NMF '%s'",
                  file_name, in_MLF);
          }
          chptr = line + strspn(line, " \n\t");
      } while (!*chptr);
  
      if ((chptr[0] != '.' || (chptr[1] != '\0' && !isspace(chptr[1])))) {
        Error("Missing '.' at the end of network '%s' in NMF '%s'",
              file_name, in_MLF);
      }
    }
    return node;
  }
  //***************************************************************************
  //***************************************************************************
  


  //***************************************************************************
  //***************************************************************************
  Node<NODE_REGULAR, LINK_REGULAR>*
  ReadSTKNetwork(
    FILE*                     lfp,
    struct MyHSearchData *    word_hash,
    struct MyHSearchData *    phone_hash,
    int                       notInDict,
    LabelFormat               labelFormat,
    long                      sampPeriod,
    const char *              file_name,
    const char *              in_MLF,
    bool                      compactRepresentation)
  {
    Node<NODE_REGULAR, LINK_REGULAR>*            node;
    Node<NODE_REGULAR, LINK_REGULAR>*            enode = NULL;
    Node<NODE_REGULAR, LINK_REGULAR>*            first = NULL;
    Node<NODE_REGULAR, LINK_REGULAR>*            last  = NULL;
    Node<NODE_REGULAR, LINK_REGULAR>*            fnode = NULL;
    Node<NODE_REGULAR, LINK_REGULAR>*            lnode;
    char*            line;
    int              line_no   =  0;
    int              nnodes    =  0;
    char*            chptr;
    char*            valptr;
    char*            phn_marks = NULL;
    Word*            word     = NULL;
    int              i;
    int              pron_var  = 1;
    enum {LINE_START, AFTER_J, HEADER_DEF, ARC_DEF, NODE_DEF} state;
    MyHSearchData node_hash = {0};
    struct ReadlineData   rld       = {0};
    
    NodeBasic<NODE_REGULAR, LINK_REGULAR>* first_basic;
  
    for (;;) 
    {
      do 
      {
        if ((chptr = line = readline(lfp, &rld)) == NULL) 
          break;
          
        if (chptr[0] == '.' && (chptr[1] == '\0' || isspace(chptr[1]))) {
          chptr = NULL;
          break;
        }
        line_no++;
        while (isspace(*chptr)) chptr++;
      } while (!*chptr || *chptr == '#');
  
      if (chptr == NULL) break; // End of file
  
      state = LINE_START;
      node = NULL;

      while (*chptr) 
      {
        for (valptr=chptr; isalnum(*valptr); valptr++);
  
        if (*valptr == '=') {
          *valptr = '\0';
          valptr++;
        } else if (!*valptr || isspace(*valptr)) { // label definition (field without '=' )
          valptr = chptr;
          chptr = "";
        } else {
          Error("Invalid character '%c' (%s:%d, char %d)",
                *valptr, file_name, line_no, valptr-line+1);
        }

        if (state == LINE_START && !strcmp(chptr, "J")) {
          getInteger(valptr, &chptr, file_name, line_no);
          state = AFTER_J;
          continue;
        }

        if (state == AFTER_J) {
          if (*chptr && strcmp(chptr,"START") && strcmp(chptr,"S")) {
            Error("Term 'J=' must be followed by term 'S=' (%s:%d)", file_name, line_no);
          }
          state = LINE_START;
          chptr="";
        }
        
        if (state == LINE_START) {
          if (!*chptr || !strcmp(chptr, "I")) {
            if (compactRepresentation)
            {
               node = reinterpret_cast<Node<NODE_REGULAR, LINK_REGULAR> *>(&first_basic[getInteger(valptr, &chptr, file_name, line_no)]);
               if (node->mType == NT_UNDEF)
                 node->mType = NT_WORD;
            }
            else
            {
              if (getHTKstr(valptr, &chptr)) {
                Error("%s (%s:%d)", chptr, file_name, line_no);
              }
              node = find_or_create_node(&node_hash, valptr, &last);
            }
            word     = NULL;
            pron_var = 1;
            state    = NODE_DEF;
            continue;
          }
          state = HEADER_DEF;
        }

        if (state == HEADER_DEF) { // label definition
          if (!strcmp(chptr, "S") || !strcmp(chptr, "SUBLAT")) {
            Error("%s not supported (%s:%d)", chptr, file_name, line_no);
          } else if (!strcmp(chptr, "N") || !strcmp(chptr, "NODES")) {
            nnodes  = getInteger(valptr, &chptr, file_name, line_no);
            if (compactRepresentation)
            {
              if (first != NULL)
                Error("Redefinition of N= (NODES=) is not allowed in CSTK format (%s:%d)", file_name, line_no);
                
              first = reinterpret_cast<Node<NODE_REGULAR, LINK_REGULAR> *>(first_basic = new NodeBasic<NODE_REGULAR, LINK_REGULAR>[nnodes]);
            }
            else if (node_hash.mTabSize == 0 && !my_hcreate_r(nnodes, &node_hash)) 
            {
              Error("Insufficient memory");
            }
          } else { // Skip unknown header term
            if (getHTKstr(valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
          }
          continue;
        }

        // node definition .....................................................
        if (state == NODE_DEF) 
        {
          if ((!strcmp(chptr, "time") || !strcmp(chptr, "t")) 
          && !(labelFormat.TIMES_OFF) && !compactRepresentation) {
            char *colonptr=valptr;
            while (*colonptr && !isspace(*colonptr) && *colonptr != ',') colonptr++;
  
            if (*colonptr == ',') {
              if (colonptr != valptr) {
                *colonptr = ' ';
                node->SetStart(100 * (long long) (0.5 + 1e5 *
                              getFloat(valptr, &chptr, file_name, line_no)));
              }
              valptr = colonptr+1;
            }
            node->SetStop(100 * (long long) (0.5 + 1e5 *
                        getFloat(valptr, &chptr, file_name, line_no)));
          } 
          else if (!strcmp(chptr, "var") || !strcmp(chptr, "v")) 
          {
            pron_var = getInteger(valptr, &chptr, file_name, line_no);
            if (pron_var < 1) {
              Error("Invalid pronunciation variant (%s:%d)", file_name, line_no);
            }
          } 
          else if (!strcmp(chptr, "p") && !compactRepresentation) 
          {
            node->mPhoneAccuracy = getFloat(valptr, &chptr, file_name, line_no);
          } 
          else if (!strcmp(chptr, "flag") || !strcmp(chptr, "f")) 
          {
            if (getHTKstr(valptr, &chptr)) 
            {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            
            for (; *valptr; valptr++) {
              switch (toupper(*valptr)) {
                case 'K':
                case 'F':  node->mType |= NT_STICKY; break;
                case 'T':  node->mType |= NT_TRUE;   break;
                default:
                  Error("Invalid flag '%c' (%s:%d)", *valptr, file_name, line_no);
              }
            }
          } else if (!strcmp(chptr, "L")) {
            Error("Sub-lattice nodes are not yet supported (%s:%d)",
                  *valptr, file_name, line_no);
          } else if (!strcmp(chptr, "WORD") || !strcmp(chptr, "W")) {
            ENTRY e = {0}; //{0} is just to make compiler happy
            ENTRY *ep;
            if (getHTKstr(e.key = valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            if (!strcmp(e.key, "!NULL")) {
              word = NULL;
            } else {
              my_hsearch_r(e, FIND, &ep, word_hash);
  
              if (ep != NULL) {
                word = static_cast<Word *>(ep->data);
              } else {
                if (notInDict & WORD_NOT_IN_DIC_ERROR) {
                  Error("Word '%s' not in dictionary (%s:%d)", e.key, file_name, line_no);
                } else if (notInDict & WORD_NOT_IN_DIC_WARN) {
                  Warning("Word '%s' not in dictionary (%s:%d)", e.key, file_name, line_no);
                }
  
                e.key  = strdup(e.key);
                word = (Word *) malloc(sizeof(Word));
                if (e.key == NULL || word  == NULL) {
                  Error("Insufficient memory");
                }
                word->mpName = e.key;
                word->npronuns = 0;
                word->npronunsInDict = 0;
                word->pronuns  = NULL;
                e.data = word;
  
                if (!my_hsearch_r(e, ENTER, &ep, word_hash)) {
                  Error("Insufficient memory");
                }
              }
            }
            node->mType &= ~(NT_MODEL | NT_PHONE);
            node->mType |= NT_WORD;
          } else if (!strcmp(chptr, "MODEL") || !strcmp(chptr, "M")) {
            ENTRY e = {0}; //{0} is just to make compiler happy
            ENTRY *ep;
  
            if (getHTKstr(e.key = valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            my_hsearch_r(e, FIND, &ep, phone_hash);
  
            if (ep == NULL) {
              e.key  = strdup(valptr);
              e.data = e.key;
  
              if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, phone_hash)) {
                Error("Insufficient memory");
              }
              ep->data = e.data;
            }
            node->mpName = (char *) ep->data;
            node->mType &= ~NT_WORD;
            node->mType |= NT_PHONE;
          } else if (*chptr=='\0' || !strcmp(chptr,"END") || !strcmp(chptr,"E")) {
            state = ARC_DEF;
          } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
            Error("%s (%s:%d)", chptr, file_name, line_no);
          }          
          if (state == ARC_DEF || *chptr == '\0') {
            // Node definition is over. For NT_WORD, select right pronun according to
            // word and pron_var; and continue with parsing the arc definition below
            if (node->mType & NT_WORD && word != NULL) {
              if (word->npronuns < pron_var) {
                // Word does not have so many pronuns; add new empty pronuns...
                if (notInDict & PRON_NOT_IN_DIC_ERROR && word->npronuns != 0) {
                  Error("Word '%s' does not have pronunciation variant %d (%s:%d)",
                        word->mpName, pron_var, file_name, line_no);
                }
  
                word->pronuns = (Pronun **) realloc(word->pronuns,
                                                    pron_var * sizeof(Pronun *));
                if (word->pronuns == NULL) Error("Insufficient memory");
  
                for (i = word->npronuns; i < pron_var; i++) {
                  word->pronuns[i] = (Pronun *) malloc(sizeof(Pronun));
                  if (word->pronuns[i] == NULL) Error("Insufficient memory");
  
                  word->pronuns[i]->mpWord     = word;
                  word->pronuns[i]->outSymbol  = word->mpName;
                  word->pronuns[i]->nmodels    = 0;
                  word->pronuns[i]->model      = NULL;
                  word->pronuns[i]->variant_no = i+1;
                  word->pronuns[i]->prob       = 0.0;
                }
                word->npronuns = pron_var;
              }
              node->mpPronun = word->pronuns[pron_var-1];
            }
            if (state == ARC_DEF) {
              // Count number of link definitions on the rest of the line and 
              // prealocate memory for links
              char *pCh;
              int nl = 1;
              
              if (skipHTKstr(valptr, &pCh)) 
                Error("%s (%s:%d)", pCh, file_name, line_no);
              
              while(*pCh != '\0') {
                if (!strncmp("END=", pCh, 4) || !strncmp("E=", pCh, 2)) nl++;
              
              
                while(isalnum(*pCh)) pCh++;
                
                if (*pCh == '=')
                {
                  if (skipHTKstr(pCh, &pCh)) 
                    Error("%s (%s:%d)", pCh, file_name, line_no);      
                    
                  continue;
                }
                
                nl++;
                if(!isspace(*pCh) && *pCh != '\0')
                  Error("Invalid character '%c' (%s:%d, char %d)",
                        *pCh, file_name, line_no, pCh-line+1);
                
                while(isspace(*pCh)) pCh++;
              }
              
              node->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) 
                realloc(node->mpLinks, (node->mNLinks + nl) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));

              // initialize the new links
              for (size_t new_i = node->mNLinks; new_i < node->mNLinks + nl; new_i++)
              {
                node->mpLinks[new_i].SetLmLike(0.0);
                node->mpLinks[new_i].SetAcousticLike(0.0);
              }

              if (node->mpLinks == NULL) Error("Insufficient memory");              
            }
          }
        }

        // arc definition ......................................................
        if (state == ARC_DEF) 
        {
          if (!*chptr || !strcmp(chptr, "END") || !strcmp(chptr, "E")) 
          {
            if (compactRepresentation)
            {
               enode = reinterpret_cast<Node<NODE_REGULAR, LINK_REGULAR> *>(&first_basic[getInteger(valptr, 
                     &chptr, file_name, line_no)]);
               if (enode->mType == NT_UNDEF)
                 enode->mType = NT_WORD;
            }
            else
            {
              if (getHTKstr(valptr, &chptr)) {
                Error("%s (%s:%d)", chptr, file_name, line_no);
              }
              enode = find_or_create_node(&node_hash, valptr, &last);
            }
  
            int nl = ++node->mNLinks;
            
            // Links are counted and node->mpLinks is properly realocated 
            // at the end of the node definition above
            node->mpLinks[nl-1].SetNode(enode);
            node->mpLinks[nl-1].SetLmLike(0.0);
            node->mpLinks[nl-1].SetAcousticLike(0.0);
  
            if (!compactRepresentation)
              ++enode->mNBackLinks;  
          } 
          
          else if (!strcmp(chptr, "language") || !strcmp(chptr, "l")) 
          {
            FLOAT lm_like = getFloat(valptr, &chptr, file_name, line_no);
            
            // Set LM score to link starting in node. This link can possibly
            // lead to a phone node already inserted (div=) between 'node' and'enode'
            node->mpLinks[node->mNLinks-1].SetLmLike(lm_like);
          } 
          
          else if (!strcmp(chptr, "acoustic") || !strcmp(chptr, "a")) 
          {
            // Set acoustic score to link starting in node. This link can possibly
            // lead to a phone node already inserted (div=) between 'node' and'enode'
            FLOAT acoustic_like = getFloat(valptr, &chptr, file_name, line_no);
            node->mpLinks[node->mNLinks-1].SetAcousticLike(acoustic_like);
          } 

          else if (!strcmp(chptr, "div") || !strcmp(chptr, "d")) 
          {
            
            if (compactRepresentation)
              Error("d= or div= is not allowed in CSTK format (%s:%d)", file_name, line_no);
            
            ENTRY e, *ep;
            char  name[1024];
            float time;
            int   n;
            Node<NODE_REGULAR, LINK_REGULAR>*  last = node;
            FLOAT lm_like  = node->mpLinks[node->mNLinks-1].LmLike();
  
            if (node->mpLinks[node->mNLinks-1].pNode() != enode) {
              Error("Redefinition of  (%s:%d)", chptr, file_name, line_no);
            }

            if (getHTKstr(phn_marks=valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }

            time = -FLT_MAX;

            while (sscanf(phn_marks, ":%[^,:]%n,%f%n", name, &n, &time, &n) > 0) 
            {
              Node<NODE_REGULAR, LINK_REGULAR> *tnode;
              phn_marks+=n;
  
              if ((tnode          = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL
              || (tnode->mpLinks  = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(
                                        sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) 
              {
                Error("Insufficient memory");
              }
  
              // initialize the new links
              tnode->mpLinks->SetLmLike(0.0);
              tnode->mpLinks->SetAcousticLike(0.0);

              //Use special type to mark nodes inserted by d=..., they will need
              //special treatment. Later, they will become ordinary NT_PHONE nodes
              tnode->mType      = NT_PHONE | NT_MODEL;
              tnode->mNLinks    = tnode->mNBackLinks = 1;
              tnode->mpBackNext  = enode->mpBackNext;
              enode->mpBackNext  = tnode;
              tnode->mPhoneAccuracy = 1.0;
              
              //Store phone durations now. Will be replaced by absolute times below.
              tnode->SetStart(time != -FLT_MAX 
                              ?  100 * (long long) (0.5 + 1e5 * time) 
                              : UNDEF_TIME);

              tnode->SetStop (UNDEF_TIME);

              e.key  = name;
              e.data = NULL;
              my_hsearch_r(e, FIND, &ep, phone_hash);
  
              if (ep == NULL) {
                e.key  = strdup(name);
                e.data = e.key;
  
                if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, phone_hash)) {
                  Error("Insufficient memory");
                }
                ep->data = e.data;
              }
              tnode->mpName = (char *) ep->data;
              last->mpLinks[last->mNLinks-1].SetNode(tnode);
              last->mpLinks[last->mNLinks-1].SetLmLike(0.0);
              last->mpLinks[last->mNLinks-1].SetAcousticLike(0.0);
              last = tnode;
            }

            if (strcmp(phn_marks,":")) {
              Error("Invalid specification of phone marks (d=) (%s:%d)",
                    file_name, line_no);
            }

            last->mpLinks[last->mNLinks-1].SetNode(enode);
            last->mpLinks[last->mNLinks-1].SetLmLike(lm_like);
          } 
          
          else if (getHTKstr(valptr, &chptr)) 
          { // Skip unknown term
            Error("%s (%s:%d)", chptr, file_name, line_no);
          }
        }
      }
    }
    
    if (compactRepresentation)
    {
      if (nnodes == 0)
        Error("No node defined in the network file (%s)", file_name);
        
      for(i = 0; i < nnodes; i++)
      {
        if (first_basic[i].mType == NT_UNDEF)
          Error("Node %d not defined in network file (%s)", i, file_name);
        
        for (int j = 0; j < first_basic[i].mNLinks; j++)
          if (first_basic[i].mpLinks[j].pNode() == first_basic)
            Error("Node 0 must be the initial network node (%s)", file_name);
      }
      
      if (first_basic[nnodes-1].mNLinks != 0)
        Error("Node with the highest id (%d) must be the final network node (%s)",
              nnodes-1, file_name);
    }
    else
    {
      if (last == NULL)
        Error("No node defined in the network file (%s)", file_name);
    
      my_hdestroy_r(&node_hash, 1);
      lnode = last;
      first = last = NULL;

      if (lnode) 
        lnode->mpNext = NULL;
    
      for (node = lnode; node != NULL; node = node->mpBackNext)
      {
        // pre-allocate space for back-links
        node->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) 
          malloc(node->mNBackLinks * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));

          // initialize the new links
          for (size_t new_i = 0; new_i < node->mNBackLinks; new_i++)
          {
            node->mpBackLinks[new_i].SetLmLike(0.0);
            node->mpBackLinks[new_i].SetAcousticLike(0.0);
          }

        if (node->mpBackLinks == NULL) 
          Error("Insufficient memory");
        
        if (node->mpBackNext) 
          node->mpBackNext->mpNext = node;
    
        if (node->mNLinks == 0) 
        {
          if (last)
            Error("Network has multiple nodes with no successors (%s)", file_name);
          last = node;
        }
        
        if (node->mNBackLinks == 0) 
        {
          if (first)
            Error("Network has multiple nodes with no predecessor (%s)", file_name);
          first = node;
        }
      
        node->mNBackLinks = 0;
      }

      if (!first || !last) {
        Error("Network contain no start node or no final node (%s)", file_name);
      }

      for (node = lnode; node != NULL; node = node->mpBackNext)
      {
        int i;
        for (i = 0; i < node->mNLinks; i++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* p_forwnode = node->mpLinks[i].pNode();
          // 
          p_forwnode->mpBackLinks[p_forwnode->mNBackLinks].SetNode(node);
          p_forwnode->mpBackLinks[p_forwnode->mNBackLinks].SetLmLike(node->mpLinks[i].LmLike());
          p_forwnode->mNBackLinks++;
        }
      }
    
      for (node = lnode; node != NULL; fnode = node, node = node->mpBackNext) 
      {
        //If only stop time is specified, set start time to lowest predecessor stop time
        if (node->Start() == UNDEF_TIME && node->Stop() != UNDEF_TIME) {
          int i;
          for (i = 0; i < node->mNBackLinks; i++) 
          {
            Node<NODE_REGULAR, LINK_REGULAR>* p_backnode = node->mpBackLinks[i].pNode();
            // skip nodes inserted by d=...
            while (p_backnode->mType == (NT_PHONE | NT_MODEL)) 
            {
              assert(p_backnode->mNBackLinks == 1);
              p_backnode = p_backnode->mpBackLinks[0].pNode();
            }

            if (p_backnode->Stop() != UNDEF_TIME) {
              node->SetStart(node->Start() == UNDEF_TIME
                             ? p_backnode->Stop()
                             : LOWER_OF(p_backnode->Stop(), node->Start()));
            }
          }
          if (node->Start() == UNDEF_TIME) 
            node->SetStart(0);
        }

        //For model nodes defined by d=... (NT_PHONE | NT_MODEL), node->Start() contains
        //only phone durations. Absolute times must be computed derived starting from
        //the end time of the node to which arc with d=... definition points.
        if (node->mType == (NT_PHONE | NT_MODEL)) 
        {
          assert(node->mNLinks == 1);

          node->SetStop (node->mpLinks[0].pNode()->mType == (NT_PHONE | NT_MODEL)
                           && node->Start() != UNDEF_TIME
                         ? node->mpLinks[0].pNode()->Start() 
                         : node->mpLinks[0].pNode()->Stop());

          node->SetStart(node->Start() != UNDEF_TIME && node->Stop() != UNDEF_TIME
                         ? node->Stop() - node->Start() 
                         : node->mpLinks[0].pNode()->Start());
        }
      }
      
      if (first != fnode) {
        if (first->mpNext)     first->mpNext->mpBackNext = first->mpBackNext;
        if (first->mpBackNext) first->mpBackNext->mpNext = first->mpNext;
        if (first == lnode)  lnode = first->mpBackNext;
    
        first->mpBackNext = NULL;
        fnode->mpBackNext = first;
        first->mpNext = fnode;
      }

      if (last != lnode) {
        if (last->mpNext)     last->mpNext->mpBackNext = last->mpBackNext;
        if (last->mpBackNext) last->mpBackNext->mpNext = last->mpNext;
        last->mpNext = NULL;
        lnode->mpNext = last;
        last->mpBackNext = lnode;
      }
      
      for (node = first; node != NULL; node = node->mpNext) {
        if (node->mType == (NT_PHONE | NT_MODEL)) {
          node->mType = NT_PHONE;
        }
        if (node->Start() != UNDEF_TIME) {
          node->SetStart((node->Start() - labelFormat.left_extent) / sampPeriod);
        }
        if (node->Stop()  != UNDEF_TIME) {
          node->SetStop((node->Stop() + labelFormat.right_extent) / sampPeriod);
        }
      }

      if (first->mpPronun != NULL) {
        node = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
        if (node == NULL) Error("Insufficient memory");
        node->mpNext       = first;
        node->mpBackNext   = NULL;
        first->mpBackNext  = node;
        node->mType       = NT_WORD;
        node->mpPronun     = NULL;
        node->SetStart(UNDEF_TIME);
        node->SetStop(UNDEF_TIME);
        node->mNBackLinks = 0;
        node->mpBackLinks  = NULL;
        node->mNLinks     = 1;
        node->mpLinks      = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (node->mpLinks == NULL) Error("Insufficient memory");
        node->mpLinks[0].SetLmLike(0.0);
        node->mpLinks[0].SetAcousticLike(0.0);
        node->mpLinks[0].SetNode(first);
        first->mNBackLinks = 1;
        first->mpBackLinks  = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (first->mpBackLinks == NULL) Error("Insufficient memory");
        first->mpBackLinks[0].SetLmLike(0.0);
        first->mpBackLinks[0].SetAcousticLike(0.0);
        first->mpBackLinks[0].SetNode(node);
        first = node;
      }
      if (last->mpPronun != NULL) {
        node = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
        if (node == NULL) Error("Insufficient memory");
        last->mpNext      = node;
        node->mpNext      = NULL;
        node->mpBackNext  = last;
        node->mType      = NT_WORD;
        node->mpPronun    = NULL;
        node->SetStart(UNDEF_TIME);
        node->SetStop(UNDEF_TIME);
        node->mNLinks    = 0;
        node->mpLinks     = NULL;
        last->mNLinks = 1;
        last->mpLinks  = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (last->mpLinks == NULL) Error("Insufficient memory");
        last->mpLinks[0].SetLmLike(0.0);
        last->mpLinks[0].SetAcousticLike(0.0);
        last->mpLinks[0].SetNode(node);
        node->mNBackLinks = 1;
        node->mpBackLinks  = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (node->mpBackLinks == NULL) Error("Insufficient memory");
        node->mpBackLinks[0].SetLmLike(0.0);
        node->mpBackLinks[0].SetAcousticLike(0.0);
        node->mpBackLinks[0].SetNode(last);
      }
    }

    return first;
  } 
  // Node *ReadSTKNetwork(...)
  //***************************************************************************


    
  //***************************************************************************
  //***************************************************************************
  Node<NODE_REGULAR, LINK_REGULAR>* 
  ReadHTKLattice(
    FILE *            lfp,
    MyHSearchData *   word_hash,
    MyHSearchData *   phone_hash,
    LabelFormat       labelFormat,
    long              sampPeriod,
    const char *      file_name)
  {
    Node<NODE_REGULAR, LINK_REGULAR> **nodes = NULL, *node,
        *first = NULL, *last = NULL,
        *fnode, *lnode = NULL;
    char line[1024];
    int nnodes  = -1;
    int nnread  =  0;
    int node_id =  0;
    int line_no =  0;
    char *chptr, *valptr, *phn_marks = NULL;
    int       arc_start = 0;
    int       arc_end   = 0;
    FLOAT     arc_like  = 0.0;
    long long node_time = UNDEF_TIME;
    int       node_var  = 0;
    Word     *node_word = NULL;
    enum {LINE_START, ARC_DEF, NODE_DEF} state;
  
    for (;;) {
      do {
        if (fgets(line, sizeof(line), lfp) == NULL) {
          chptr = NULL;
          break;
        }
        line_no++;
        if (strlen(line) == sizeof(line) - 1) {
          Error("Line is too long (%s:%d)", file_name, line_no);
        }
        chptr = line + strspn(line, " \n\t");
      } while (!*chptr || *chptr == '#');
  
      if (chptr == NULL) break; // End of file
  
      state = LINE_START;
      while (*chptr) {
        if ((valptr = strchr(chptr, '=')) == NULL) {
          Error("'=' expected (%s:%d)", file_name, line_no);
        }
        valptr++;
  
        if (state == LINE_START) {
  
          if (chptr[0] == 'S') {
            Error("%s not supported (%s:%d)", chptr, file_name, line_no);
          }
          if (chptr[0] == 'N') {
            if (nnodes > -1) {
              Error("Redefinition of number of nodes (%s:%d)",file_name,line_no);
            }
            nnodes = strtoull(valptr, NULL, 10);
            nodes  = (Node<NODE_REGULAR, LINK_REGULAR> **) calloc(nnodes, sizeof(Node<NODE_REGULAR, LINK_REGULAR> *));
            if (nodes == NULL) Error("Insufficient memory");
            break;
          }
          if (chptr[0] == 'I') {
            state = NODE_DEF;
            node_id = getNodeNumber(nnodes, valptr, &chptr, file_name, line_no);
  
            if (nodes[node_id] != NULL) {
              Error("Redefinition of node %d (%s:%d)",
              node_id, file_name, line_no);
            }
            nodes[node_id] = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
            if (nodes[node_id] == NULL) Error("Insufficient memory");
  
            if (last == NULL) {
              first = nodes[node_id];
            } else {
              last->mpNext = nodes[node_id];
            }
            last = nodes[node_id];
            nodes[node_id]->mpNext = NULL;
  
            node_word = NULL;
            node_var  = 0;
            node_time = UNDEF_TIME;
            nnread++;
          } else if (chptr[0] == 'J') {
            if (nnodes < nnread) {
              for (node_id = 0; nodes[node_id] != NULL; node_id++);
              Error("Definition of node %d is missing (%s:%d)",
                    node_id, file_name, line_no);
            }
            state = ARC_DEF;
            getInteger(valptr, &chptr, file_name, line_no);
            arc_like  = 0.0;
            arc_start = -1;
            arc_end   = -1;
            phn_marks = NULL;
          } else {
            break; // Ignore line with unknown initial term
          }
        } else if (state == NODE_DEF) {
          if (chptr[0] == 't' && !(labelFormat.TIMES_OFF)) {
            node_time = 100 * (long long)
                        (0.5 + 1e5*getFloat(valptr, &chptr, file_name, line_no));
          } else if (chptr[0] == 'v') {
            node_var = getInteger(valptr, &chptr, file_name, line_no) - 1;
          } else if (chptr[0] == 'W') {
            ENTRY e = {0}; //{0} is just to make compiler happy
            ENTRY *ep;
            if (getHTKstr(e.key = valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
  
            if (strcmp(e.key, "!NULL")) {
              my_hsearch_r(e, FIND, &ep, word_hash);
  
              if (ep == NULL) {
                Error("Unknown word '%s' (%s:%d)", e.key, file_name, line_no);
              }
              node_word = (Word *) ep->data;
            }
          } else {
            if (getHTKstr(valptr, &chptr)) { // Skip unknown term
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
          }
        } else {
          assert(state == ARC_DEF);
          if (chptr[0] == 'S') {
            arc_start = getNodeNumber(nnodes, valptr, &chptr, file_name,line_no);
          } else if (chptr[0] == 'E') {
            arc_end   = getNodeNumber(nnodes, valptr, &chptr, file_name,line_no);
          } else if (chptr[0] == 'l') {
            arc_like  = getFloat(valptr, &chptr, file_name, line_no);
          } else if (chptr[0] == 'd') {
            if (getHTKstr(phn_marks=valptr, &chptr)) { // Skip unknown term
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
          } else {
            if (getHTKstr(valptr, &chptr)) { // Skip unknown term
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
          }
        }
      }
      if (state == NODE_DEF) {
        nodes[node_id]->mType       = NT_WORD;
        nodes[node_id]->mNLinks     = 0;
        nodes[node_id]->mNBackLinks = 0;
        nodes[node_id]->mpLinks      = NULL;
        nodes[node_id]->mpBackLinks  = NULL;
        if (node_word && node_word->npronuns <= node_var) {
          Error("Word %s does not have pronunciation varian %d (%s:%d)",
                node_word->mpName, node_var+1, file_name, line_no);
        }
        nodes[node_id]->mpPronun     = node_word ? node_word->pronuns[node_var]
                                              : NULL;
        nodes[node_id]->SetStart(UNDEF_TIME);
        nodes[node_id]->SetStop(node_time + labelFormat.right_extent);
      } else if (state == ARC_DEF) {
        if (arc_start == -1 || arc_end == -1) {
          Error("Start node or end node not defined (%s:%d)", file_name, line_no);
        }
  
        int linkId = nodes[arc_start]->mNLinks++;
        nodes[arc_start]->mpLinks =
          (Link<NODE_REGULAR, LINK_REGULAR> *) realloc(nodes[arc_start]->mpLinks, (linkId+1) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
  
        if (nodes[arc_start]->mpLinks == NULL) Error("Insufficient memory");
  
        last = nodes[arc_start];
  
        if (phn_marks) {
          ENTRY e, *ep;
          char  name[1024];
          float time;
          int   n;
  
          while (sscanf(phn_marks, ":%[^,],%f%n", name, &time, &n) > 1) {
            phn_marks+=n;
  
            if ((node            = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
              (node->mpLinks     = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
              (node->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) {
              Error("Insufficient memory");
            }
            node->mType = NT_PHONE;
            node->mNLinks = node->mNBackLinks = 1;
            node->mpNext   = last->mpNext;
            last->mpNext = node;
            node->mPhoneAccuracy = 1.0;
  
            if (!(labelFormat.TIMES_OFF)) {
              node->SetStart(last->Stop() - labelFormat.left_extent - labelFormat.right_extent);
              node->SetStop( last->Stop() + 100 * (long long) (0.5 + 1e5 * time)
                                      + labelFormat.right_extent);
            } else {
              node->SetStart(UNDEF_TIME);
              node->SetStop(UNDEF_TIME);
            }

            e.key  = name;
            e.data = NULL;
            my_hsearch_r(e, FIND, &ep, phone_hash);
  
            if (ep == NULL) {
              e.key  = strdup(name);
              e.data = e.key;
  
              if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, phone_hash)) {
                Error("Insufficient memory");
              }
              ep->data = e.data;
            }
            node->mpName = (char *) ep->data;
            last->mpLinks[last->mNLinks-1].SetNode(node);
            last->mpLinks[last->mNLinks-1].SetLmLike(0.0);
            node->mpBackLinks[0].SetNode(last);
            node->mpBackLinks[0].SetLmLike(0.0);
            last = node;
          }
          if (strcmp(phn_marks,":")) {
            Error("Invalid specification of phone marks (d=) (%s:%d)",
                  file_name, line_no);
          }
        }
        linkId = nodes[arc_end]->mNBackLinks++;
        nodes[arc_end]->mpBackLinks =
          (Link<NODE_REGULAR, LINK_REGULAR> *) realloc(nodes[arc_end]->mpBackLinks, (linkId+1) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
  
        if (nodes[arc_end]->mpBackLinks == NULL) Error("Insufficient memory");
  
        last->mpLinks[last->mNLinks-1].SetNode(nodes[arc_end]);
        last->mpLinks[last->mNLinks-1].SetLmLike(arc_like);
        nodes[arc_end]->mpBackLinks[linkId].SetNode(last);
        nodes[arc_end]->mpBackLinks[linkId].SetLmLike(arc_like);
  
        if (nodes[arc_start]->Stop() != UNDEF_TIME) 
        {
          if (nodes[arc_end]->Start() == UNDEF_TIME) 
          {
            nodes[arc_end]->SetStart( nodes[arc_start]->Stop() 
                - labelFormat.left_extent - labelFormat.right_extent);
          } 
          else 
          {
            nodes[arc_end]->SetStart( LOWER_OF(nodes[arc_end]->Start(), 
                  nodes[arc_start]->Stop() - labelFormat.left_extent 
                  - labelFormat.right_extent));
          }
        }
      }
    }
    free(nodes);
    fnode = first;
    first = last = NULL;
    if (fnode) fnode->mpBackNext = NULL;
    for (node = fnode; node != NULL; lnode = node, node = node->mpNext) {
      if (node->mpNext) node->mpNext->mpBackNext = node;
  
      if (node->mNLinks == 0) {
        if (last)
          Error("Network has multiple nodes with no successors (%s)", file_name);
        last = node;
      }
      if (node->mNBackLinks == 0) {
        if (first)
          Error("Network has multiple nodes with no predecessor (%s)", file_name);
        first = node;
      }
    }
    if (!first || !last) {
      Error("Network contain no start node or no final node (%s)", file_name);
    }
    if (first != fnode) {
      if (first->mpNext) first->mpNext->mpBackNext = first->mpBackNext;
      first->mpBackNext->mpNext = first->mpNext;
      first->mpBackNext = NULL;
      fnode->mpBackNext = first;
      first->mpNext = fnode;
    }
    if (last != lnode) {
      if (last->mpBackNext) last->mpBackNext->mpNext = last->mpNext;
      last->mpNext->mpBackNext = last->mpBackNext;
      last->mpNext = NULL;
      lnode->mpNext = last;
      last->mpBackNext = lnode;
    }
    if (first->mpPronun != NULL) {
      node = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
      if (node == NULL) Error("Insufficient memory");
      node->mpNext       = first;
      node->mpBackNext   = NULL;
      first->mpBackNext  = node;
      node->mType       = NT_WORD;
      node->mpPronun     = NULL;
      node->SetStart(UNDEF_TIME);
      node->SetStop(UNDEF_TIME);
      node->mNBackLinks = 0;
      node->mpBackLinks  = NULL;
      node->mNLinks     = 1;
      node->mpLinks      = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
      if (node->mpLinks == NULL) Error("Insufficient memory");
      node->mpLinks[0].SetLmLike(0.0);
      node->mpLinks[0].SetNode(first);
  
      first->mNBackLinks = 1;
      first->mpBackLinks  = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
      if (first->mpBackLinks == NULL) Error("Insufficient memory");
      first->mpBackLinks[0].SetLmLike(0.0);
      first->mpBackLinks[0].SetNode(node);
      first = node;
    }
    if (last->mpPronun != NULL) {
      node = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
      if (node == NULL) Error("Insufficient memory");
      last->mpNext      = node;
      node->mpNext      = NULL;
      node->mpBackNext  = last;
  
      node->mType      = NT_WORD;
      node->mpPronun    = NULL;
      node->SetStart(UNDEF_TIME);
      node->SetStop(UNDEF_TIME);
      node->mNLinks    = 0;
      node->mpLinks     = NULL;
  
      last->mNLinks = 1;
      last->mpLinks  = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
      if (last->mpLinks == NULL) Error("Insufficient memory");
      last->mpLinks[0].SetLmLike(0.0);
      last->mpLinks[0].SetNode(node);
  
      node->mNBackLinks = 1;
      node->mpBackLinks  = (Link<NODE_REGULAR, LINK_REGULAR>*) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
      if (node->mpBackLinks == NULL) Error("Insufficient memory");
      node->mpBackLinks[0].SetLmLike(0.0);
      node->mpBackLinks[0].SetNode(last);
    }
    return first;
  }
  
  
  
} // namespace STK



