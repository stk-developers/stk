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
//#define VERSION "0.2 "__TIME__" "__DATE__


// PROJECT INCLUDES
//
#include "common.h"
#include "Models.h"
#include "labels.h"
#include "dict.h"
#include "Net.h"

#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <malloc.h>
#include <assert.h>
#include <ctype.h>
#include <stdarg.h>



#define SIGNIFICANT_PROB_DIFFERENCE (0.01)


// CODE
//

namespace STK
{
  //***************************************************************************
  //***************************************************************************

  void 
  FreeNetwork(Node<NODE_REGULAR, LINK_REGULAR> * pNode, bool compactRepresentation) 
  {
    if (!compactRepresentation)
    {
      Node<NODE_REGULAR, LINK_REGULAR> *  tnode;
      while (pNode) 
      {
        tnode = pNode->mpNext;
        free(pNode->mpLinks);
        free(pNode->mpBackLinks);
        free(pNode);
        pNode = tnode;
      }
    }
    else
    {
      NodeBasic<NODE_REGULAR, LINK_REGULAR> *  tnode;
      for(tnode =  reinterpret_cast<NodeBasic<NODE_REGULAR, LINK_REGULAR>*>(pNode); tnode->mNLinks != 0; tnode++) 
        free(tnode->mpLinks);
      
      free(pNode);
    }
  }




  //***************************************************************************
  //***************************************************************************
  Node<NODE_REGULAR, LINK_REGULAR> *
  MakeNetworkFromLabels(Label * pLabels, NodeKind nodeKind)
  {
    Label * lp;
    Node<NODE_REGULAR, LINK_REGULAR>  * first;
    Node<NODE_REGULAR, LINK_REGULAR>  * last = NULL;
    Node<NODE_REGULAR, LINK_REGULAR>  * node;
  
    if ((first             = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
        (last              = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
        (first->mpLinks    = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL    ||
        (last->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) 
    {
      Error("Insufficient memory");
    }

    first->mpName       = last->mpName      = NULL;     
    first->mType        = last->mType       = NT_WORD;
    first->mpPronun     = last->mpPronun    = NULL;
    first->mNLinks      = last->mNBackLinks = 1;
    first->mNBackLinks  = last->mNLinks     = 0;
    first->mpBackLinks  = last->mpLinks     = NULL;
    first->SetStart(UNDEF_TIME);
    first->SetStop(UNDEF_TIME);
    
    last->SetStart(UNDEF_TIME);
    last->SetStop(UNDEF_TIME);

  //  first->mpTokens        = last->mpTokens        = NULL;
  //  first->mpExitToken     = last->mpExitToken     = NULL;
  
    node = first;
    
    for (lp = pLabels; lp != NULL; lp = lp->mpNext) 
    {
      Node<NODE_REGULAR, LINK_REGULAR> *  tnode;
  
      if ((tnode              = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
          (tnode->mpLinks     = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL    ||
          (tnode->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) 
      {
        Error("Insufficient memory");
      }
      
      tnode->mpName = NULL;
      
      node->mpLinks[0].SetNode(tnode);
      node->mpLinks[0].SetLmLike(0.0);
      node->mpLinks[0].SetAcousticLike(0.0);
      
      switch (nodeKind) 
      {
        case NT_WORD:  tnode->mpPronun = ((Word *) lp->mpData)->pronuns[0]; break;
        case NT_MODEL: tnode->mpHmm    =   (Hmm *) lp->mpData;              break;
        case NT_PHONE: tnode->mpName   =  (char *) lp->mpData;              break;
        default:       Error("Fatal: Invalid node type");
      }
      
      tnode->mType       = nodeKind;
      tnode->mNLinks     = 1;
      tnode->mNBackLinks = 1;
      tnode->SetStart(lp->mStart);
      tnode->SetStop(lp->mStop);
      tnode->mpBackLinks[0].SetNode(node);
      tnode->mpBackLinks[0].SetLmLike(0.0);
      tnode->mpBackLinks[0].SetAcousticLike(0.0);
      node->mpNext = tnode;
      node = tnode;
    }
    
    node->mpNext = last;
    node->mpLinks[0].SetNode(last);
    node->mpLinks[0].SetLmLike(0.0);
    node->mpLinks[0].SetAcousticLike(0.0);
    last->mpBackLinks[0].SetNode(node);
    last->mpBackLinks[0].SetLmLike(0.0);
    last->mpBackLinks[0].SetAcousticLike(0.0);
    last->mpNext = NULL;
    return first;
  }
  
  
  void ExpandWordNetworkByDictionary(
    Node<NODE_REGULAR, LINK_REGULAR> *        pFirstNode,
    MyHSearchData *dict,
    int keep_word_nodes,
    int multiple_pronun)
  {
    Node<NODE_REGULAR, LINK_REGULAR>* node; 
    Node<NODE_REGULAR, LINK_REGULAR>* prev = NULL;
    int   i;
    int   j;
  
    Pronun    singlePronun;
    Pronun*   singlePronunPtr;
    Word      singlePronunWrd;

    singlePronunWrd.npronuns = 1;
    singlePronunWrd.pronuns  = &singlePronunPtr;
    singlePronunPtr = &singlePronun;
  
    assert(pFirstNode != NULL || pFirstNode->mType & NT_WORD || pFirstNode->mpPronun == NULL);
  
    for (node = pFirstNode; node != NULL; prev = node, node = node->mpNext) 
    {
      if (!(node->mType & NT_WORD)) continue;
  
      if (node->mpPronun == NULL) continue;
      Word *word = node->mpPronun->mpWord;
  
      //Do not expand non-dictionary words, which where added by ReadSTKNetwork
      if (word->npronunsInDict == 0) continue;
  
      if (!multiple_pronun) 
      {
        singlePronunWrd.mpName = node->mpPronun->mpWord->mpName;
        word = &singlePronunWrd;
        *word->pronuns = node->mpPronun;
      }
  
      // Remove links to current node form backlinked nodes and realloc
      // link arrays of backlinked nodes to hold word->npronuns more backlinks
      for (i = 0; i < node->mNBackLinks; i++) 
      {
        Node<NODE_REGULAR, LINK_REGULAR>* bakcnode = node->mpBackLinks[i].pNode();

        for (j=0; j<bakcnode->mNLinks && bakcnode->mpLinks[j].pNode() != node; j++)
        {}

        assert(j < bakcnode->mNLinks); // Otherwise link to 'node' is missing
                                      // from which backlink exists
        bakcnode->mpLinks[j] = bakcnode->mpLinks[bakcnode->mNLinks-1];
  
        bakcnode->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *)
          realloc(bakcnode->mpLinks,
                (bakcnode->mNLinks - 1 + word->npronuns) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (bakcnode->mpLinks == NULL) Error("Insufficient memory");
        bakcnode->mNLinks--;// += word->npronuns-1;
      }
  
      // Remove backlinks to current node form linked nodes and realloc
      // backlink arrays of linked nodes to hold word->npronuns more backlinks
      for (i=0; i < node->mNLinks; i++) 
      {
        Node<NODE_REGULAR, LINK_REGULAR>* forwnode = node->mpLinks[i].pNode();
        for (j=0;j<forwnode->mNBackLinks&&forwnode->mpBackLinks[j].pNode()!=node;j++);
        assert(j < forwnode->mNBackLinks);
        // Otherwise link to 'node' is missing from which backlink exists
        forwnode->mpBackLinks[j] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
  
        forwnode->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *)
          realloc(forwnode->mpBackLinks,
                (forwnode->mNBackLinks - 1 + word->npronuns) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
        forwnode->mNBackLinks--;
      }

      for (i = 0; i < word->npronuns; i++) 
      {
        Pronun *pronun = word->pronuns[i];
        Node<NODE_REGULAR, LINK_REGULAR> *pronun_first = NULL, *pronun_prev = NULL, *tnode;
  
        for (j = 0; j < pronun->nmodels; j++) 
        {
          tnode = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
          
          if (tnode == NULL) 
            Error("Insufficient memory");
  
          tnode->mType       = NT_PHONE | (node->mType & NT_TRUE);
          tnode->mpName      = pronun->model[j].mpName;
          tnode->SetStart(node->Start());
          tnode->SetStop (node->Stop());
          tnode->mPhoneAccuracy = 1.0;
  
          if (j == 0) 
          {
            pronun_first = tnode;
          } 
          else 
          {
            if ((pronun_prev->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
                (tnode->mpBackLinks   = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) 
            {
              Error("Insufficient memory");
            }
            
            tnode->mNBackLinks              = 1;
            tnode->mpBackLinks[0].SetNode(pronun_prev);
            tnode->mpBackLinks[0].SetLmLike(0.0);
            tnode->mpBackLinks[0].SetAcousticLike(0.0);
            pronun_prev->mNLinks            = 1;
            pronun_prev->mpLinks[0].SetNode(tnode);
            pronun_prev->mpLinks[0].SetLmLike(0.0);
            pronun_prev->mpLinks[0].SetAcousticLike(0.0);
            pronun_prev->mpNext             = tnode;
          }
          pronun_prev = tnode;
        }
        
        if (keep_word_nodes || j == 0) 
        {
          tnode = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>));
          
          if (tnode == NULL) Error("Insufficient memory");
  
          tnode->mpName   = NULL;
          tnode->mType    = NT_WORD | (node->mType & NT_TRUE);
          tnode->mpPronun = keep_word_nodes ? word->pronuns[i] : NULL;
          tnode->SetStart(node->Start());
          tnode->SetStop (node->Stop());
  
          if (j == 0) 
          {
            pronun_first = tnode;
          } 
          else 
          {
            if ((pronun_prev->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
              (tnode->mpBackLinks   = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) {
              Error("Insufficient memory");
            }

            tnode->mNBackLinks          = 1;
            tnode->mpBackLinks[0].SetNode(pronun_prev);
            tnode->mpBackLinks[0].SetLmLike(0.0);
            tnode->mpBackLinks[0].SetAcousticLike(0.0);
            
            pronun_prev->mNLinks        = 1;
            pronun_prev->mpLinks[0].SetNode(tnode);
            pronun_prev->mpLinks[0].SetLmLike(0.0);
            pronun_prev->mpLinks[0].SetAcousticLike(0.0);

            pronun_prev->mpNext         = tnode;
          }
          pronun_prev = tnode;
        }
        if ((pronun_prev->mpLinks =
              (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>) * node->mNLinks))==NULL ||
          (pronun_first->mpBackLinks =
              (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>) * node->mNBackLinks)) == NULL) {
          Error("Insufficient memory");
        }
        pronun_prev->mNLinks      = node->mNLinks;
        pronun_first->mNBackLinks = node->mNBackLinks;
  
        for (j = 0; j < node->mNBackLinks; j++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* backnode = node->mpBackLinks[j].pNode();
          backnode->mpLinks[backnode->mNLinks].SetNode(pronun_first);
          backnode->mpLinks[backnode->mNLinks].SetLmLike(node->mpBackLinks[j].LmLike());
          backnode->mpLinks[backnode->mNLinks].SetAcousticLike(node->mpBackLinks[j].AcousticLike());
          backnode->mNLinks++;
          pronun_first->mpBackLinks[j] = node->mpBackLinks[j];
        }

        for (j=0; j < node->mNLinks; j++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* forwnode = node->mpLinks[j].pNode();
          forwnode->mpBackLinks[forwnode->mNBackLinks].SetNode(pronun_prev);
          forwnode->mpBackLinks[forwnode->mNBackLinks].SetLmLike(node->mpLinks[j].LmLike());
          forwnode->mpBackLinks[forwnode->mNBackLinks].SetAcousticLike(node->mpLinks[j].AcousticLike());
          forwnode->mNBackLinks++;
          pronun_prev->mpLinks[j] = node->mpLinks[j];
        }
        if (prev != NULL) prev->mpNext = pronun_first;
        prev = pronun_prev;
      }
      prev->mpNext = node->mpNext;
      free(node->mpLinks);
      free(node->mpBackLinks);
      free(node);
      node = prev;
    }
  }

  
  //****************************************************************************
  //****************************************************************************
  Node<NODE_REGULAR, LINK_REGULAR>* 
  DiscardUnwantedInfoInNetwork(Node<NODE_REGULAR, LINK_REGULAR> *pFirstNode, STKNetworkOutputFormat format)
  // The function discard the information in network records that is not to be
  // saved to the output. This should allow for more effective network
  // optimization, which will be run after calling this function and before
  // saving network to file.
  {
    Node<NODE_REGULAR, LINK_REGULAR>* node;
    int   i;
  
    for (node = pFirstNode; node != NULL; node = node->mpNext)  
    {
      if (format.mNoLMLikes) 
      {
        for (i=0; i < node->mNLinks;     i++) 
          node->mpLinks    [i].SetLmLike(0.0);
        
        for (i=0; i < node->mNBackLinks; i++) 
          node->mpBackLinks[i].SetLmLike(0.0);
      }

      if (format.mNoAcousticLikes) 
      {
        for (i=0; i < node->mNLinks;     i++) 
          node->mpLinks    [i].SetAcousticLike(0.0);

        for (i=0; i < node->mNBackLinks; i++) 
          node->mpBackLinks[i].SetAcousticLike(0.0);
      }

      if (format.mNoTimes) {
        node->SetStart(UNDEF_TIME);
        node->SetStop(UNDEF_TIME);
      }
      if (format.mNoWordNodes && node->mType & NT_WORD) {
        node->mpPronun = NULL;
      }
      if (format.mNoModelNodes && (node->mType&NT_MODEL || node->mType&NT_PHONE)) {
        node->mType = NT_WORD;
        node->mpPronun = NULL;
      }
      if (format.mNoPronunVars && node->mType & NT_WORD && node->mpPronun != NULL) {
        node->mpPronun = node->mpPronun->mpWord->pronuns[0];
      }
    }
    return pFirstNode;
  }

  

  //***************************************************************************
  //***************************************************************************
  static int 
  LatticeLocalOptimization_ForwardPass(Node<NODE_REGULAR, LINK_REGULAR>* pFirstNode, int strictTiming)
  {
    int i, j, k, l, m, rep;
    Node<NODE_REGULAR, LINK_REGULAR> *node, *tnode;
    int node_removed = 0;
    FLOAT tlike;

    for (node = pFirstNode; node != NULL; node = node->mpNext) 
    {
/**/      for (i = 0; i < node->mNLinks; i++) 
      {
      //for (tnode = inode; tnode != NULL; tnode = (tnode == inode ? jnode : NULL)) {
      
        tnode = node->mpLinks[i].pNode();
        
        if (tnode->mNLinks == 0) 
          continue;
  
        // Weight pushing
        tlike = tnode->mpBackLinks[0].LmLike();
        
        for (l=1; l <  tnode->mNBackLinks; l++) 
          tlike = HIGHER_OF(tlike, tnode->mpBackLinks[l].LmLike());
        
        for (l=0; l < tnode->mNBackLinks; l++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* backnode = tnode->mpBackLinks[l].pNode();
          tnode->mpBackLinks[l].AddLmLike(-tlike);
          
          for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].pNode()!=tnode; k++)
          {}
          
          assert(k < backnode->mNLinks);
          backnode->mpLinks[k].AddLmLike(-tlike);
#ifndef NDEBUG
          for (k++; k<backnode->mNLinks && backnode->mpLinks[k].pNode()!=tnode; k++)
          {}
#endif
          assert(k == backnode->mNLinks);
        }
        
        for (l=0; l < tnode->mNLinks; l++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* forwnode = tnode->mpLinks[l].pNode();
          tnode->mpLinks[l].AddLmLike(tlike);
          
          for (k=0; k<forwnode->mNBackLinks && forwnode->mpBackLinks[k].pNode()!=tnode;k++)
          {}
          
          assert(k < forwnode->mNBackLinks);
          forwnode->mpBackLinks[k].AddLmLike(tlike);
#ifndef NDEBUG
          for (k++; k<forwnode->mNBackLinks && forwnode->mpBackLinks[k].pNode()!=tnode;k++)
          {}
#endif
          assert(k == forwnode->mNBackLinks);
        }
      }
/**/      
  //dnet(pFirstNode, 1, node);
  
      // For current node 'node', check for each possible pair of its successors
      // ('inode' and 'jnode') whether the pair may be merged to single node.
      for (i = 0; i < node->mNLinks-1; i++) 
      {
        for (j = i+1; j < node->mNLinks; j++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR> *inode = node->mpLinks[i].pNode();
          Node<NODE_REGULAR, LINK_REGULAR> *jnode = node->mpLinks[j].pNode();

          // Final node may be never merged.
          if (inode->mNLinks == 0 || jnode->mNLinks == 0) 
            continue;


          // Two nodes ('inode' and 'jnode') may be mergeg if they are of the 
          // same type, name, ... with the same predecessors and with the same
          // weights on the links from predecesors.
          if ((inode->mType & ~NT_TRUE) != (jnode->mType & ~NT_TRUE)
          || ( inode->mType & NT_PHONE && inode->mpName   != jnode->mpName)
          || ( inode->mType & NT_WORD  && inode->mpPronun != jnode->mpPronun)

  //          &&  (inode->mpPronun == NULL ||
  //             jnode->mpPronun == NULL ||
  //             inode->mpPronun->mpWord       != jnode->mpPronun->mpWord ||
  //             inode->mpPronun->outSymbol  != jnode->mpPronun->outSymbol ||
  //             inode->mpPronun->variant_no != jnode->mpPronun->variant_no ||
  //             inode->mpPronun->prob       != jnode->mpPronun->prob)
          || (inode->mNBackLinks != jnode->mNBackLinks)) 
          {
            continue;
          }
          
          if (strictTiming && (inode->Start() != jnode->Start()
                          ||  inode->Stop()  != jnode->Stop())) 
          {
            continue;
          }

//Weights on the links from predecesors does not have to be exactely the same, but the must not
//differ more than by SIGNIFICANT_PROB_DIFFERENCE
          for (l=0; l < inode->mNBackLinks; l++) 
          {
            if (inode->mpBackLinks[l].pNode() != jnode->mpBackLinks[l].pNode()) break;
            FLOAT ldiff =  inode->mpBackLinks[l].LmLike() - jnode->mpBackLinks[l].LmLike();
            if (ldiff < -SIGNIFICANT_PROB_DIFFERENCE ||
              ldiff >  SIGNIFICANT_PROB_DIFFERENCE) 
            {
              break;
            }
          }
          
          if (l < inode->mNBackLinks) 
            continue;
  
  /*        if (memcmp(inode->mpBackLinks, jnode->mpBackLinks,
                    inode->mNBackLinks * sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) {
            continue;
          }*/
            // inode and jnode are the same nodes with the same predeccessors
            // Remove jnode and add its links to inode
  
          assert(inode->mNLinks && jnode->mNLinks);
  
//TraceLog("Removing node: %s", inode->mType & NT_PHONE ? inode->mpName : 
//                              inode->mType & NT_WORD  ? (inode->mpPronun ? inode->mpPronun->mpWord->mpName : "!NULL") : "UNKNOWN_TYPE");


            // Remove links to jnode form predeccessors
          for (l=0; l < jnode->mNBackLinks; l++) 
          {
            Node<NODE_REGULAR, LINK_REGULAR> *backnode = jnode->mpBackLinks[l].pNode();
            for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].pNode()!=jnode; k++);
            assert(k < backnode->mNLinks);
            // Otherwise link to 'node' is missing from which backlink exists
            memmove(backnode->mpLinks+k, backnode->mpLinks+k+1,
                    (backnode->mNLinks-k-1) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
            backnode->mNLinks--;
          }
  
          // Merge jnode's links with inode links
  
          //Count jnode's links not present among inode's links
          rep = l = k = 0;
          while (k < jnode->mNLinks) 
          {
            Link<NODE_REGULAR, LINK_REGULAR> *ill = inode->mpLinks+l;
            Link<NODE_REGULAR, LINK_REGULAR> *jlk = jnode->mpLinks+k;
            if (l == inode->mNLinks || ill->pNode() > jlk->pNode())
            {
              // k-th link of jnode will be included among inode's links.
              // Redirect corresponding baclink to inode
              for (m = 0; m < jlk->pNode()->mNBackLinks
                          && jlk->pNode()->mpBackLinks[m].pNode() != jnode; m++)
              {}
              
              assert(m < jlk->pNode()->mNBackLinks);
              jlk->pNode()->mpBackLinks[m].SetNode(node);
              qsort(jlk->pNode()->mpBackLinks, jlk->pNode()->mNBackLinks,
                    sizeof(Link<NODE_REGULAR, LINK_REGULAR>), lnkcmp);
              k++;
            } 
            else  if (ill->pNode() == jlk->pNode()) 
            {
              // l-th link of inode and k-th link of jnode points to
              // the same node. Link<NODE_REGULAR, LINK_REGULAR> from jnode is redundant.
              // Remove backlinks to jnode form jnode's succesors
              for (m = 0; m < jlk->pNode()->mNBackLinks
                        && jlk->pNode()->mpBackLinks[m].pNode() != jnode; m++);
              {}
              
              assert(m < jlk->pNode()->mNBackLinks);
              memmove(jlk->pNode()->mpBackLinks+m, jlk->pNode()->mpBackLinks+m+1,
                      (jlk->pNode()->mNBackLinks-m-1) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
              jlk->pNode()->mNBackLinks--;
  
              ill->SetLmLike(HIGHER_OF(ill->LmLike(), jlk->LmLike()));
              jlk->SetNode(NULL); // Mark link to be removed
              rep++; k++, l++;
            } 
            else 
            {
              l++;
            }
          }
          
          l = inode->mNLinks;
          inode->mNLinks += jnode->mNLinks-rep;
          inode->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) realloc(inode->mpLinks,
                                          inode->mNLinks * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
          
          if (inode->mpLinks == NULL) 
            Error("Insufficient memory");
  
          for (k = 0; k < jnode->mNLinks; k++) 
          {
            if (jnode->mpLinks[k].pNode() != NULL) 
              inode->mpLinks[l++] = jnode->mpLinks[k];
          }
          
          qsort(inode->mpLinks, inode->mNLinks, sizeof(Link<NODE_REGULAR, LINK_REGULAR>), lnkcmp);
  
          inode->SetStart(inode->Start() == UNDEF_TIME || jnode->Start() == UNDEF_TIME
                          ? UNDEF_TIME 
                          : LOWER_OF(inode->Start(), jnode->Start()));
  
          inode->SetStop (inode->Stop() == UNDEF_TIME || jnode->Stop() == UNDEF_TIME
                          ? UNDEF_TIME 
                          : HIGHER_OF(inode->Stop(), jnode->Stop()));
  
          if (inode->mAux > jnode->mAux) 
          {
          // Make sure that topological order of new inode's links
          // (inherited from jnode) is higher than inode's order.
          // In the 'next' list, move inode to jnode's lower position
            inode->mpBackNext->mpNext = inode->mpNext;
            inode->mpNext->mpBackNext = inode->mpBackNext;
            inode->mpNext           = jnode->mpNext;
            inode->mpBackNext       = jnode->mpBackNext;
            inode->mpBackNext->mpNext = inode;
            inode->mpNext->mpBackNext = inode;
            inode->mAux = jnode->mAux;
          } 
          else
          {
            jnode->mpNext->mpBackNext = jnode->mpBackNext;
            jnode->mpBackNext->mpNext = jnode->mpNext;
          }
          
          inode->mType |= jnode->mType & NT_TRUE;
          free(jnode->mpLinks);
          free(jnode->mpBackLinks);
          free(jnode);
          --j; // Process j-th node again
              // there is new shifted node on this index
  
          node_removed = 1;
        }
      }
    }
    return node_removed;
  }
  
  //***************************************************************************
  //***************************************************************************
  Node<NODE_REGULAR, LINK_REGULAR>*
  ReverseNetwork(Node<NODE_REGULAR, LINK_REGULAR> *pFirstNode)
  {
    Node<NODE_REGULAR, LINK_REGULAR>*  node;
    Node<NODE_REGULAR, LINK_REGULAR>*  last = NULL;
    
    for (node = pFirstNode; node != NULL; node = node->mpBackNext) 
    {
      Link<NODE_REGULAR, LINK_REGULAR> *  links       = node->mpLinks;
      int     nlinks      = node->mNLinks;
      node->mpLinks       = node->mpBackLinks;
      node->mNLinks       = node->mNBackLinks;
      node->mpBackLinks   = links;
      node->mNBackLinks   = nlinks;
      Node<NODE_REGULAR, LINK_REGULAR> *  next        = node->mpNext;
      node->mpNext        = node->mpBackNext;
      node->mpBackNext    = next;
      node->mAux          = -node->mAux;
      last                = node;
    }
    return last;
  }
  
  //***************************************************************************
  //***************************************************************************
  static int LatticeLocalOptimization_BackwardPass(Node<NODE_REGULAR, LINK_REGULAR> *pFirstNode, int strictTiming)
  {
    int node_removed;
    Node<NODE_REGULAR, LINK_REGULAR> *last = ReverseNetwork(pFirstNode);
    node_removed = LatticeLocalOptimization_ForwardPass(last, strictTiming);
  //  if (!node_removed) dnet(last, 0);
    ReverseNetwork(last);
  //  if (!node_removed) dnet(pFirstNode, 0);
    return node_removed;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int nbacklinkscmp(const void *a, const void *b) 
  {
    return (*(Node<NODE_REGULAR, LINK_REGULAR> **) a)->mNBackLinks - (*(Node<NODE_REGULAR, LINK_REGULAR> **) b)->mNBackLinks;
  }
  
  //***************************************************************************
  //***************************************************************************
  void LatticeLocalOptimization(Node<NODE_REGULAR, LINK_REGULAR> *pFirstNode, int strictTiming, int trace_flag)
  {
    Node<NODE_REGULAR, LINK_REGULAR> *    node;
    Node<NODE_REGULAR, LINK_REGULAR> *    lastnode;
    int       i;
    int       j;
    int       unreachable = 0;
    
    // For each node, sort links by pointer value to allow
    // for easy comparison whether two nodes have the same set of links
    for (node = pFirstNode; node != NULL; node = node->mpNext)  
    {
      node->mAux = 0;
      node->mpBackNext = node->mpNext;
      qsort(node->mpLinks, node->mNLinks, sizeof(Link<NODE_REGULAR, LINK_REGULAR>), lnkcmp);
      qsort(node->mpBackLinks, node->mNBackLinks, sizeof(Link<NODE_REGULAR, LINK_REGULAR>), lnkcmp);
    }
  
  
    // Sort nodes in topological order
    // printf("Sorting nodes...\n");
    pFirstNode->mAux = 1;
    for (lastnode = node = pFirstNode; node != NULL; node = node->mpNext) 
    {
      for (i=0; i < node->mNLinks; i++) 
      {
        Node<NODE_REGULAR, LINK_REGULAR>* lnknode = node->mpLinks[i].pNode();
        
        if (lnknode->mAux == 0) 
        {
          for (j=0; j<lnknode->mNBackLinks && lnknode->mpBackLinks[j].pNode()->mAux==1; j++)
          {}
          
          if (j == lnknode->mNBackLinks) 
          {
            lastnode->mpNext = lnknode;
            lastnode  = lnknode;
            lnknode->mAux = 1;
            lnknode->mpNext = NULL;
          }
        }
      }
    }
    
    if (lastnode->mNLinks != 0) 
    {
      // There is a cycle in graph so we cannot sort nodes
      // topologicaly, so sort it at least somehow. :o|
      // Anyway this optimization algorithm is not optimal for graphs with cycles.
      for (node = pFirstNode; node != NULL; node = node->mpBackNext) 
      {
        node->mAux = 0;
      }
        
      pFirstNode->mAux = 1;
      for (lastnode = node = pFirstNode; node != NULL; node = node->mpNext) 
      {
        for (i=0; i < node->mNLinks; i++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* lnknode = node->mpLinks[i].pNode();
          if (lnknode->mAux == 0) 
          {
            lastnode->mpNext = lnknode;
            lastnode  = lnknode;
            lnknode->mAux = 1;
            lnknode->mpNext = NULL;
          }
        }
      }
      
      for (node=pFirstNode; node->mpNext->mNLinks != 0; node=node->mpNext)
      {}
  
      // Final node is not at the and of chain
      if (node->mpNext->mpNext) 
      { 
        lastnode->mpNext = node->mpNext;
        node->mpNext = node->mpNext->mpNext;
        lastnode = lastnode->mpNext;
        lastnode->mpNext = NULL;
      }
    }
  
    // !!! Unreachable nodes must be removed before sorting !!!
  
    for (node=pFirstNode; node != NULL; node=node->mpBackNext) 
    {
      while (node->mpBackNext && node->mpBackNext->mAux == 0) 
      {
        Node<NODE_REGULAR, LINK_REGULAR>* tnode = node->mpBackNext;
        node->mpBackNext = node->mpBackNext->mpBackNext;
        unreachable++;
        free(tnode->mpLinks);
        free(tnode->mpBackLinks);
        free(tnode);
      }
    }
  
    //  if (unreachable) Warning("Removing %d unreachable nodes", unreachable);
    if (unreachable) 
      Error("Networks contains unreachable nodes");
  
    pFirstNode->mpBackNext = NULL;
    
    for (node=pFirstNode; node->mpNext != NULL; node=node->mpNext) 
    {
      node->mpNext->mpBackNext = node;
    }
    
    for (i=1, node=pFirstNode; node != NULL; node = node->mpNext, i++) 
    {
      node->mAux=i;
    }
    
    for (;;) 
    {
      if (trace_flag & 2) 
      {
        for (i=0,node=pFirstNode; node; node=node->mpNext,i++)
        {}
        
        TraceLog("Forward pass.... (number of nodes: %d)", i);
      }
      
      LatticeLocalOptimization_ForwardPass(pFirstNode, strictTiming);
  
      if (trace_flag & 2) 
      {
        for (i=0,node=pFirstNode; node; node=node->mpNext,i++)
        {}
        
        TraceLog("Backward pass... (number of nodes: %d)", i);
      }
      
      if (!LatticeLocalOptimization_BackwardPass(pFirstNode, strictTiming)) 
        break;
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void ExpandMonophoneNetworkToTriphones(
    Node<NODE_REGULAR, LINK_REGULAR> *  pFirstNode,
    MyHSearchData *nonCDphones,
    MyHSearchData *CDphones)
  {
    Node<NODE_REGULAR, LINK_REGULAR> *  node;
    int     did_we_clone;
    int     i;
    int     j;
    int     k;
    // Find all Tee model, Word, and Null nodes (except the pFirstNode and last Null node)
    // and clone those not having single input and output
  
    do 
    {
      ENTRY     e    = {0}; //{0} is just to make compiler happy
      ENTRY *   ep;
      Node<NODE_REGULAR, LINK_REGULAR> *    prev = NULL;
      
      did_we_clone = 0;
      
      for (node = pFirstNode; node != NULL; prev = node, node = node->mpNext) 
      {
        if ((node->mNLinks == 0 || node->mNBackLinks == 0) ||
            (node->mNLinks == 1 && node->mNBackLinks == 1)) 
        {
          continue;
        }
        
        if (node->mType & NT_PHONE) 
        {
          e.key = node->mpName;
          my_hsearch_r(e, FIND, &ep, nonCDphones);
          if (ep == NULL || !reinterpret_cast<size_t>(ep->data))
            continue; // Node<NODE_REGULAR, LINK_REGULAR> is not a Tee model
        }
        
        did_we_clone = 1;
        assert(prev != NULL); //Otherwise pFirstNode node is not Null node
  
        // Remove links to current node form back-linked nodes and realloc
        // link arrays of back-linked nodes to hold node->mNLinks more links
        for (j=0; j < node->mNBackLinks; j++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR> *backnode = node->mpBackLinks[j].pNode();
          
          for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].pNode()!=node; k++)
          {}
          
          assert(k < backnode->mNLinks);
          
          // Otherwise link to 'node' is missing from which backlink exists
          backnode->mpLinks[k] = backnode->mpLinks[backnode->mNLinks-1];
  
          backnode->mpLinks = 
            (Link<NODE_REGULAR, LINK_REGULAR> *) realloc((backnode->mpLinks), 
                             (backnode->mNLinks-1+node->mNLinks)*sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
          
          if (backnode->mpLinks == NULL) 
            Error("Insufficient memory");
          
          backnode->mNLinks--;
        }
        
        // Remove backlinks to current node form linked nodes and realloc
        // backlink arrays of linked nodes to hold node->mNBackLinks more backlinks
        for (j=0; j < node->mNLinks; j++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* forwnode = node->mpLinks[j].pNode();
          
          for (k=0;k<forwnode->mNBackLinks&&forwnode->mpBackLinks[k].pNode()!=node;k++)
          {}
          
          assert(k < forwnode->mNBackLinks);
          // Otherwise link to 'node' is missing from which backlink exists
          forwnode->mpBackLinks[k] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
  
          forwnode->mpBackLinks = 
            (Link<NODE_REGULAR, LINK_REGULAR> *) realloc((forwnode->mpBackLinks),
                             (forwnode->mNBackLinks-1+node->mNBackLinks)*sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
                             
          if (forwnode->mpBackLinks == NULL) 
            Error("Insufficient memory");
            
          forwnode->mNBackLinks--;
        }
        // Alloc new node->mNLinks * node->mNBackLinks nodes and create new links
        // so that each backlinked node is conected with each linked node through
        // one new node.
        for (i=0; i < node->mNLinks; i++) 
        {
          for (j=0; j < node->mNBackLinks; j++) 
          {
            Node<NODE_REGULAR, LINK_REGULAR> *  tnode;
            Link<NODE_REGULAR, LINK_REGULAR>    forwlink = node->mpLinks[i];
            Link<NODE_REGULAR, LINK_REGULAR>    backlink = node->mpBackLinks[j];
  
            if ((tnode = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL)
              Error("Insufficient memory");
            
            tnode->mpName = NULL;
            *tnode = *node;
  
            if ((tnode->mpLinks     = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
              (tnode->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) 
            {
              Error("Insufficient memory");
            }
            
            tnode->mNLinks       = 1;
            tnode->mNBackLinks   = 1;
            tnode->mpLinks[0]     = forwlink;
            tnode->mpBackLinks[0] = backlink;
            
            forwlink.pNode()->mpBackLinks[forwlink.pNode()->mNBackLinks].SetNode(tnode);
            forwlink.pNode()->mpBackLinks[forwlink.pNode()->mNBackLinks].SetLmLike(forwlink.LmLike());
            forwlink.pNode()->mpBackLinks[forwlink.pNode()->mNBackLinks].SetAcousticLike(forwlink.AcousticLike());
            forwlink.pNode()->mNBackLinks++;

            backlink.pNode()->mpLinks    [backlink.pNode()->mNLinks    ].SetNode(tnode);
            backlink.pNode()->mpLinks    [backlink.pNode()->mNLinks    ].SetLmLike(backlink.LmLike());
            backlink.pNode()->mpLinks    [backlink.pNode()->mNLinks    ].SetAcousticLike(backlink.AcousticLike());
            backlink.pNode()->mNLinks++;
            
            prev->mpNext = tnode;
            prev = tnode;
          }
        }
        prev->mpNext = node->mpNext;
        free(node->mpLinks);
        free(node->mpBackLinks);
        free(node);
        node = prev;
      }
    } while (did_we_clone);
  
    // Assign to each node unique number, which will later allow to find groups of
    // expanded triphone nodes corresponding to original monophone nodes.
    int nbackmononodes;
    int nforwmononodes;
    int id = 0;
    
    for (node = pFirstNode; node != NULL; node = node->mpNext) 
      node->mAux = id++;
  
    // Expand monophone nodes to triphone nodes
    Node<NODE_REGULAR, LINK_REGULAR> *prev = NULL;
    for (node = pFirstNode; node != NULL; prev = node, node = node->mpNext) 
    {
      ENTRY e = {0}; //{0} is just to make compiler happy
      ENTRY * ep;
  
      if ((node->mType & NT_WORD) ||
          (node->mNLinks == 1 && node->mNBackLinks == 1)) 
      {
        continue;
      }
      
      assert(node->mType & NT_PHONE);
      e.key = node->mpName;
      my_hsearch_r(e, FIND, &ep, nonCDphones);
      if (ep != NULL && reinterpret_cast<size_t>(ep->data)) continue; // Node<NODE_REGULAR, LINK_REGULAR> is a Tee model
  
      assert(prev != NULL); //Otherwise first node is not Null node
  
      // Count groups of backlinked nodes corresponding to different monophones
      id = -1;
      nbackmononodes = 0;
      for (j=0; j < node->mNBackLinks; j++) 
      {
        if (node->mpBackLinks[j].pNode()->mAux != id) 
        {
          id = node->mpBackLinks[j].pNode()->mAux;
          nbackmononodes++;
        }
      }
      // Count groups of linked nodes corresponding to different monophones
      id = -1;
      nforwmononodes = 0;
      for (j=0; j < node->mNLinks; j++) 
      {
        if (node->mpLinks[j].pNode()->mAux != id) 
        {
          id = node->mpLinks[j].pNode()->mAux;
          nforwmononodes++;
        }
      }
  
      // Remove links to current node form backlinked nodes and realloc
      // link arrays of backlinked nodes to hold nforwmononodes more links
      for (j=0; j < node->mNBackLinks; j++) 
      {
        Node<NODE_REGULAR, LINK_REGULAR>* backnode = node->mpBackLinks[j].pNode();
        for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].pNode()!=node; k++);
        assert(k < backnode->mNLinks);
        // Otherwise link to 'node' is missing from which backlink exists
        memmove(backnode->mpLinks+k, backnode->mpLinks+k+1,
                (backnode->mNLinks-k-1) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
  
        backnode->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *)
          realloc(backnode->mpLinks,
                (backnode->mNLinks-1+nforwmononodes)*sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (backnode->mpLinks == NULL) Error("Insufficient memory");
        backnode->mNLinks--;
      }
  
      // Remove backlinks to current node form linked nodes and realloc
      // backlink arrays of linked nodes to hold nbackmononodes more backlinks
      for (j=0; j < node->mNLinks; j++) 
      {
        Node<NODE_REGULAR, LINK_REGULAR>* forwnode = node->mpLinks[j].pNode();
        for (k=0;k<forwnode->mNBackLinks&&forwnode->mpBackLinks[k].pNode()!=node;k++);
        assert(k < forwnode->mNBackLinks);
        // Otherwise link to 'node' is missing from which backlink exists
        memmove(forwnode->mpBackLinks+k, forwnode->mpBackLinks+k+1,
                (forwnode->mNBackLinks-k-1) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
  
        forwnode->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *)
          realloc(forwnode->mpBackLinks,
                (forwnode->mNBackLinks-1+nbackmononodes)*sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
        if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
        forwnode->mNBackLinks--;
      }
  
      // Alloc new nforwmononodes * nbackmononodes nodes and create new links
      // so that each backlinked node is conected through one new node with all
      // linked nodes belonging to one monophone group and vice versa each
      // linked node is conected through one new node with all backlinked nodes
      // belonging to one monophone group
      Link<NODE_REGULAR, LINK_REGULAR> *forwmono_start, *forwmono_end = node->mpLinks;
      for (i=0; i < nforwmononodes; i++) 
      {
        for (forwmono_start = forwmono_end;
            forwmono_end < node->mpLinks+node->mNLinks &&
            forwmono_start->pNode()->mAux == forwmono_end->pNode()->mAux;
            forwmono_end++)
        {}
  
        assert((i <  nforwmononodes-1 && forwmono_end <  node->mpLinks+node->mNLinks) ||
              (i == nforwmononodes-1 && forwmono_end == node->mpLinks+node->mNLinks));
  
        Link<NODE_REGULAR, LINK_REGULAR> *tlink, *backmono_start, *backmono_end = node->mpBackLinks;
        
        for (j=0; j < nbackmononodes; j++) 
        {
          for (backmono_start = backmono_end;
            backmono_end < node->mpBackLinks+node->mNBackLinks &&
            backmono_start->pNode()->mAux == backmono_end->pNode()->mAux;
            backmono_end++)
          {}
  
          assert((j <  nbackmononodes-1 && backmono_end <  node->mpBackLinks+node->mNBackLinks) ||
                 (j == nbackmononodes-1 && backmono_end == node->mpBackLinks+node->mNBackLinks));
  
          Node<NODE_REGULAR, LINK_REGULAR> * tnode;
          if ((tnode = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL)
            Error("Insufficient memory");
          
          *tnode = *node;
          tnode->mNLinks       = forwmono_end-forwmono_start;
          tnode->mNBackLinks   = backmono_end-backmono_start;
  
          if ((tnode->mpLinks =
              (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(tnode->mNLinks * sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
            (tnode->mpBackLinks =
              (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(tnode->mNBackLinks * sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) 
          {
            Error("Insufficient memory");
          }
          
          for (tlink = forwmono_start; tlink < forwmono_end; tlink++) 
          {
            tnode->mpLinks[tlink-forwmono_start] = *tlink;
            tlink->pNode()->mpBackLinks[tlink->pNode()->mNBackLinks].SetNode(tnode);
            tlink->pNode()->mpBackLinks[tlink->pNode()->mNBackLinks].SetLmLike(tlink->LmLike());
            tlink->pNode()->mpBackLinks[tlink->pNode()->mNBackLinks].SetAcousticLike(tlink->AcousticLike());
            tlink->pNode()->mNBackLinks++;
          }
          
          for (tlink = backmono_start; tlink < backmono_end; tlink++) 
          {
            tnode->mpBackLinks[tlink-backmono_start] = *tlink;
            tlink->pNode()->mpLinks[tlink->pNode()->mNLinks].SetNode(tnode);
            tlink->pNode()->mpLinks[tlink->pNode()->mNLinks].SetLmLike(tlink->LmLike());
            tlink->pNode()->mpLinks[tlink->pNode()->mNLinks].SetAcousticLike(tlink->AcousticLike());
            tlink->pNode()->mNLinks++;
          }
          
          prev->mpNext = tnode;
          prev = tnode;
        }
      }
      prev->mpNext = node->mpNext;
      free(node->mpLinks);
      free(node->mpBackLinks);
      free(node);
      node = prev;
    }
  
    // Give triphone names to phone nodes and create hash of these names
    for (node = pFirstNode; node != NULL; node = node->mpNext) 
    {
      ENTRY     e       = {0}; //{0} is just to make compiler happy
      ENTRY*    ep      = NULL;
      Node<NODE_REGULAR, LINK_REGULAR>*     lc      = NULL;
      Node<NODE_REGULAR, LINK_REGULAR>*     rc      = NULL;
      char*     lcname  = NULL;
      char*     rcname  = NULL;
      char*     triname = NULL;
      int       lcnlen  = 0;
      int       rcnlen  = 0;
  
      if (!(node->mType & NT_PHONE)) 
        continue;
  
      if (nonCDphones) 
      {
        e.key  = node->mpName;
        my_hsearch_r(e, FIND, &ep, nonCDphones);
      } 
      else 
      {
        ep = NULL;
      }
      
      if (ep != NULL) 
      {
        lc = rc = NULL;
      } 
      else 
      {
        for (lc = node;;) 
        {
          lc = lc->mNBackLinks ? lc->mpBackLinks[0].pNode() : NULL;
          
          if (lc == NULL)               break;
          if (!(lc->mType & NT_PHONE))  continue;
          if (nonCDphones == NULL)      break;
          
          e.key  = lc->mpName;
          my_hsearch_r(e, FIND, &ep, nonCDphones);
          if (ep == NULL || !reinterpret_cast<size_t>(ep->data)) break; // Node represents Tee model
        }
        
        for (rc = node;;) 
        {
          rc = rc->mNLinks ? rc->mpLinks[0].pNode() : NULL;
          
          if (rc == NULL)               break;
          if (!(rc->mType & NT_PHONE))  continue;
          if (nonCDphones == NULL)      break;
          
          e.key  = rc->mpName;
          my_hsearch_r(e, FIND, &ep, nonCDphones);
          if (ep == NULL || !reinterpret_cast<size_t>(ep->data)) break; // Node represents Tee model
        }
      }
      
      lcnlen = -1;
      if (lc != NULL) 
      {
        lcname = strrchr(lc->mpName, '-');
        
        if (lcname == NULL) 
          lcname = lc->mpName;
        else 
          lcname++;
          
        lcnlen = strcspn(lcname, "+");
      }
      
      rcnlen = -1;
      if (rc != NULL) 
      {
        rcname = strrchr(rc->mpName, '-');
        
        if (rcname == NULL) 
          rcname = rc->mpName;
        else 
          rcname++;
          
        rcnlen = strcspn(rcname, "+");
      }
      
      triname = (char *) malloc(lcnlen+1+strlen(node->mpName)+1+rcnlen+1);
      
      if (triname == NULL) 
        Error("Insufficient memory");
  
      triname[0] = '\0';
  
      if (lcnlen > 0) 
        strcat(strncat(triname, lcname, lcnlen), "-");
        
      strcat(triname, node->mpName);
      if (rcnlen > 0) 
        strncat(strcat(triname, "+"), rcname, rcnlen);
  
      e.key  = triname;
      my_hsearch_r(e, FIND, &ep, CDphones);
  
      if (ep == NULL) 
      {
        e.key  = triname;
        e.data = e.key;
  
        if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, CDphones))
          Error("Insufficient memory");
        
        node->mpName = triname;
      } 
      else 
      {
        free(triname);
        node->mpName = ep->key;
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  // Insert null node to self links. Self links are not supported
  // in network manipulation functions  
  void SelfLinksToNullNodes(Node<NODE_REGULAR, LINK_REGULAR> * pFirstNode)
  {
    int   i;
    int   j;
    Node<NODE_REGULAR, LINK_REGULAR>* node;
    Node<NODE_REGULAR, LINK_REGULAR>* tnode;
  
    for (node = pFirstNode; node != NULL; node = node->mpNext) 
    {
      for (i=0; i < node->mNLinks; i++) 
      {
        if (node->mpLinks[i].pNode() == node) 
        {
          if ((tnode           = (Node<NODE_REGULAR, LINK_REGULAR> *) calloc(1, sizeof(Node<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
              (tnode->mpLinks     = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL ||
              (tnode->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *) malloc(sizeof(Link<NODE_REGULAR, LINK_REGULAR>))) == NULL) 
          {
            Error("Insufficient memory");
          }
  
          tnode->mpName = NULL;
          node->mpLinks[i].SetNode(tnode);
          
          for (j=0; j<node->mNBackLinks && node->mpBackLinks[j].pNode()!=node; j++)
          {}
          
          assert(j<node->mNBackLinks);
          
          node->mpBackLinks[j].SetNode(tnode);
          node->mpBackLinks[j].SetLmLike(0.0);
          node->mpBackLinks[j].SetAcousticLike(0.0);
  
          tnode->mType       = NT_WORD;
          tnode->mpPronun     = NULL;
          tnode->mNLinks     = 1;
          tnode->mNBackLinks = 1;
          tnode->SetStart(UNDEF_TIME);
          tnode->SetStop (UNDEF_TIME);
  //        tnode->mpTokens     = NULL;
  //        tnode->mpExitToken  = NULL;
          tnode->mpLinks[0].SetNode(node);
          tnode->mpLinks[0].SetLmLike(0.0);
          tnode->mpLinks[0].SetAcousticLike(0.0);
          tnode->mpBackLinks[0].SetNode(node);
          tnode->mpBackLinks[0].SetLmLike(node->mpLinks[i].LmLike());
          tnode->mpBackLinks[0].SetAcousticLike(node->mpLinks[i].AcousticLike());
          tnode->mpNext = node->mpNext;
          node->mpNext = tnode;
        }
      }
    }
  }

  
  //***************************************************************************
  //***************************************************************************
  // Remove null nones having less than three predecessors or less than three successors
  int RemoveRedundantNullNodes(Node<NODE_REGULAR, LINK_REGULAR> *pFirstNode)
  {
    Node<NODE_REGULAR, LINK_REGULAR> *    node;
    Node<NODE_REGULAR, LINK_REGULAR> *    tnode;
    int       i;
    int       j;
    int       k;
    int       node_removed = 0;
  
    pFirstNode->mpBackNext = NULL;
    for (node = pFirstNode; node->mpNext != NULL; node = node->mpNext) {
      node->mpNext->mpBackNext = node;
    }
    for (node = pFirstNode; node != NULL; node = node->mpNext) {
      if (node->mType & NT_WORD && node->mpPronun == NULL &&
          node->mNLinks != 0 && node->mNBackLinks != 0  &&
        (node->mNLinks == 1 || node->mNBackLinks == 1 ||
        (node->mNLinks == 2 && node->mNBackLinks == 2))) {
  
      node_removed = 1;
  
      // Remove links to current node form backlinked nodes and realloc
      // link arrays of backlinked nodes to hold node->mNLinks more backlinks
        for (i = 0; i < node->mNBackLinks; i++) 
        {
          Node<NODE_REGULAR, LINK_REGULAR>* bakcnode = node->mpBackLinks[i].pNode();

          for (j=0; j<bakcnode->mNLinks && bakcnode->mpLinks[j].pNode()!=node; j++)
          {}

          assert(j < bakcnode->mNLinks); // Otherwise link to 'node' is missing
                                        // from which backlink exists
          bakcnode->mpLinks[j] = bakcnode->mpLinks[bakcnode->mNLinks-1];
  
          bakcnode->mpLinks = (Link<NODE_REGULAR, LINK_REGULAR> *)
            realloc(bakcnode->mpLinks,
                  (bakcnode->mNLinks - 1 + node->mNLinks) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
          if (bakcnode->mpLinks == NULL) Error("Insufficient memory");
          bakcnode->mNLinks--;// += word->npronuns-1;
        }
  
        // Remove backlinks to current node form linked nodes and realloc
        // backlink arrays of linked nodes to hold word->npronuns more backlinks
        for (i=0; i < node->mNLinks; i++) {
          Node<NODE_REGULAR, LINK_REGULAR> *forwnode = node->mpLinks[i].pNode();
          for (j=0;j<forwnode->mNBackLinks&&forwnode->mpBackLinks[j].pNode()!=node;j++);
          assert(j < forwnode->mNBackLinks);
          // Otherwise link to 'node' is missing from which backlink exists
          forwnode->mpBackLinks[j] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
  
          forwnode->mpBackLinks = (Link<NODE_REGULAR, LINK_REGULAR> *)
            realloc(forwnode->mpBackLinks,
                  (forwnode->mNBackLinks - 1 + node->mNBackLinks) * sizeof(Link<NODE_REGULAR, LINK_REGULAR>));
          if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
          forwnode->mNBackLinks--;
        }
        for (j = 0; j < node->mNBackLinks; j++) {
          Node<NODE_REGULAR, LINK_REGULAR> *backnode = node->mpBackLinks[j].pNode();
          int orig_nlinks = backnode->mNLinks;
  
          for (i=0; i < node->mNLinks; i++) {
            for(k = 0; k < orig_nlinks && backnode->mpLinks[k].pNode() != node->mpLinks[i].pNode(); k++);
            if(k < orig_nlinks) {
              // Link<NODE_REGULAR, LINK_REGULAR> which is to be created already exists. Its duplication must be avoided.
              backnode->mpLinks[k].SetLmLike(HIGHER_OF(backnode->mpLinks[k].LmLike(), 
                                             node->mpLinks[i].LmLike() + node->mpBackLinks[j].LmLike()));
            } else {
              backnode->mpLinks[backnode->mNLinks  ].SetNode(node->mpLinks[i].pNode());
              backnode->mpLinks[backnode->mNLinks++].SetLmLike(
                  node->mpLinks[i].LmLike() + node->mpBackLinks[j].LmLike());
            }
          }
        }
        for (j = 0; j < node->mNLinks; j++) {
          Node<NODE_REGULAR, LINK_REGULAR> *forwnode = node->mpLinks[j].pNode();
          int orig_nbacklinks = forwnode->mNBackLinks;
  
          for (i=0; i < node->mNBackLinks; i++) {
            for(k = 0; k < orig_nbacklinks && forwnode->mpBackLinks[k].pNode() != node->mpBackLinks[i].pNode(); k++);
            if (k < orig_nbacklinks) {
              // Link<NODE_REGULAR, LINK_REGULAR> which is to be created already exists. Its duplication must be avoided.
              forwnode->mpBackLinks[k].SetLmLike(HIGHER_OF(forwnode->mpBackLinks[k].LmLike(), 
                                                      node->mpBackLinks[i].LmLike() + node->mpLinks[j].LmLike()));
            } else {
              forwnode->mpBackLinks[forwnode->mNBackLinks  ].SetNode(node->mpBackLinks[i].pNode());
              forwnode->mpBackLinks[forwnode->mNBackLinks++].SetLmLike(
                  node->mpBackLinks[i].LmLike() + node->mpLinks[j].LmLike());
            }
          }
        }
        node->mpBackNext->mpNext = node->mpNext;
        node->mpNext->mpBackNext = node->mpBackNext;
        tnode = node;
        node = node->mpBackNext;
        free(tnode->mpLinks);
        free(tnode->mpBackLinks);
        free(tnode);
      }
    }
    return node_removed;
  }
  
  struct CorrPhnRec 
  {
    Node<NODE_REGULAR, LINK_REGULAR>      *   mpNode;
    long long     maxStopTimeTillNow;
    int           mId;
  };
  
  //***************************************************************************
  //***************************************************************************
  static int cmp_starts(const void *a, const void *b)
  {
    long long diff = ((struct CorrPhnRec *) a)->mpNode->Start()
                  - ((struct CorrPhnRec *) b)->mpNode->Start();
  
    if (diff != 0)
      return diff;
  
    return ((struct CorrPhnRec *) a)->mpNode->Stop()
        - ((struct CorrPhnRec *) b)->mpNode->Stop();
  }
  
  //***************************************************************************
  //***************************************************************************
  static int cmp_maxstop(const void *key, const void *elem)
  {
    struct CorrPhnRec *corr_phn = (struct CorrPhnRec *) elem;
  
    if (((Node<NODE_REGULAR, LINK_REGULAR> *) key)->Start() < corr_phn->maxStopTimeTillNow) {
      if (corr_phn->mId == 0 || // first fiead in the array
        ((Node<NODE_REGULAR, LINK_REGULAR> *) key)->Start() >= (corr_phn-1)->maxStopTimeTillNow)
          return  0;
      else return -1;
    } else return  1;
  }
  
  //***************************************************************************
  //***************************************************************************
  // Check whether strings a and b represents the same monophone
  int SamePhoneme(char *a, char *b)
  {
    char *  chptr;
    size_t  len;
  
    chptr = strrchr(a, '-');
    if (chptr) a = chptr+1;
    len = strcspn(a, "+");
  
    chptr = strrchr(b, '-');
    
    if (chptr) 
      b = chptr+1;
      
    if (strcspn(b, "+") != len) 
      return 0;
  
    return strncmp(a, b, len) == 0;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void ComputeAproximatePhoneAccuracy(Node<NODE_REGULAR, LINK_REGULAR> *pFirstNode, int type)
  {
    Node<NODE_REGULAR, LINK_REGULAR> *node;
    struct CorrPhnRec *corr_phn, *overlaped;
    int ncorr_phns = 0, i = 0;
    long long maxStopTime;
  
    for (node = pFirstNode; node != NULL; node = node->mpNext) {
      if (node->mType & NT_PHONE && node->mType & NT_TRUE) ncorr_phns++;
    }
    if (ncorr_phns == 0) Error("No correct phoneme node in network");
  
    corr_phn = (struct CorrPhnRec *)malloc(sizeof(struct CorrPhnRec)*ncorr_phns);
    if (corr_phn == NULL) Error("Insufficient memory");
  
    for (node = pFirstNode; node != NULL; node = node->mpNext) {
      if (node->mType & NT_PHONE && node->mType & NT_TRUE) {
        corr_phn[i++].mpNode = node;
      }
    }
  
    qsort(corr_phn, ncorr_phns, sizeof(struct CorrPhnRec), cmp_starts);
  
    maxStopTime = corr_phn[0].mpNode->Stop();
    for (i = 0; i < ncorr_phns; i++) {
      corr_phn[i].mId = i;
      maxStopTime = HIGHER_OF(maxStopTime, corr_phn[i].mpNode->Stop());
      corr_phn[i].maxStopTimeTillNow = maxStopTime;
    }
    for (node = pFirstNode; node != NULL; node = node->mpNext) {
      if (!(node->mType & NT_PHONE)) continue;
  
      if (node->Stop()  <= node->Start() ||
        !strcmp(node->mpName, "sil") ||
        !strcmp(node->mpName, "sp")) {
        //!!! List of ignored phonemes should be provided by some switch !!!
        node->mPhoneAccuracy = 0.0;
      } else {
        node->mPhoneAccuracy = -1.0;
        overlaped = (CorrPhnRec*)bsearch(node, corr_phn, ncorr_phns,
                            sizeof(struct CorrPhnRec), 
          cmp_maxstop);
  
        if (overlaped) {
          for (; overlaped < corr_phn + ncorr_phns &&
                overlaped->mpNode->Start() < node->Stop(); overlaped++) {
            if (overlaped->mpNode->Stop()  <= overlaped->mpNode->Start() ||
              overlaped->mpNode->Stop()  <= node->Start()) continue;
  
            node->mPhoneAccuracy =
              HIGHER_OF(node->mPhoneAccuracy,
                        (SamePhoneme(overlaped->mpNode->mpName, node->mpName) + 1.0) *
                        (LOWER_OF(overlaped->mpNode->Stop(), node->Stop()) -
                         HIGHER_OF(overlaped->mpNode->Start(), node->Start())) /
                        (overlaped->mpNode->Stop() - overlaped->mpNode->Start()) - 1.0);
  
            if (node->mPhoneAccuracy >= 1.0) break;
          }
        }
      }
    }
    free(corr_phn);
  }  

#ifndef NDEBUG
  // Debug function showing network using AT&T dot utility
  void dnet(Node<NODE_REGULAR, LINK_REGULAR> *net, int nAuxNodePtrs, ...)
  {
    static int dnetcnt=1;
    va_list ap;
    Node<NODE_REGULAR, LINK_REGULAR> *node;
    int i = 1;
  
    FILE *fp = popen("cat | (tf=`mktemp /tmp/netps.XXXXXX`;"
                    "dot -Tps > $tf; gv -scale -4 $tf; rm $tf)",
                    "w");
  //  FILE *fp = stdout;
  
    if (fp == NULL) return;
  
    for (node = net; node != NULL; node = node->mpNext) {
      node->mEmittingStateId = i++;
    }
    fprintf(fp, "digraph \"dnet%d\" {\nrankdir=LR\n", dnetcnt++);
  
  
    for (node = net; node != NULL; node = node->mpNext) {
      fprintf(fp, "n%d [shape=%s,label=\"%d:%s", node->mEmittingStateId,
              node->mType & NT_WORD ? "box" : "ellipse", node->mEmittingStateId,
              node->mType & NT_WORD ? (node->mpPronun ?
                                      node->mpPronun->mpWord->mpName : "-"):
              node->mType & NT_PHONE? node->mpName :
              node->mType & NT_MODEL? node->mpHmm->mpMacro->mpName : "???");
  
      if (node->mType & NT_WORD && node->mpPronun != NULL) {
        if (node->mpPronun != node->mpPronun->mpWord->pronuns[0]) {
          fprintf(fp, ":%d", node->mpPronun->variant_no);
        }
        fprintf(fp, "\\n");
  
        if (node->mpPronun->outSymbol != node->mpPronun->mpWord->mpName) {
          fprintf(fp, "[%s]", node->mpPronun->outSymbol ?
                              node->mpPronun->outSymbol : "");
        }
        if (node->mpPronun->prob != 0.0) {
          fprintf(fp, " "FLOAT_FMT, node->mpPronun->prob);
        }
      }
      fprintf(fp, "\"];\n");
  //    if (node->mpNext != NULL) {
  //     fprintf(fp,"n%d -> n%d [color=black,weight=1]\n",
  //             node->mEmittingStateId, node->mpNext->mEmittingStateId);
  //    }
  //    if (node->mpBackNext != NULL) {
  //     fprintf(fp,"n%d -> n%d [color=gray,weight=1]\n",
  //             node->mEmittingStateId, node->mpBackNext->mEmittingStateId);
  //    }
    }
    for (node = net; node != NULL; node = node->mpNext) {
      for (i = 0; i < node->mNLinks; i++) {
        fprintf(fp,"n%d -> n%d [color=blue,weight=1",
                node->mEmittingStateId,node->mpLinks[i].pNode()->mEmittingStateId);
        if (node->mpLinks[i].LmLike() != 0.0) {
          fprintf(fp,",label=\""FLOAT_FMT"\"", node->mpLinks[i].LmLike());
        }
        fprintf(fp,"];\n");
      }
  //    for (i = 0; i < node->mNBackLinks; i++) {
  //      fprintf(fp,"n%d -> n%d [color=red,weight=1",
  //              node->mEmittingStateId,node->mpBackLinks[i].pNode()->mEmittingStateId);
  //      if (node->mpBackLinks[i].mLmLike != 0.0) {
  //        fprintf(fp,",label=\""FLOAT_FMT"\"", node->mpBackLinks[i].mLmLike);
  //      }
  //      fprintf(fp,"];\n");
  //    }
    }
    va_start(ap, nAuxNodePtrs);
    typedef  Node<NODE_REGULAR, LINK_REGULAR>  my_node;
    for (i = 0; i < nAuxNodePtrs; i++) {
      my_node* ptr = va_arg(ap, my_node* );
      fprintf(fp, "AuxPtr%d [shape=plaintext];\nAuxPtr%d -> n%d\n",
              i, i, ptr->mEmittingStateId);
    }
    va_end(ap);
  
    fprintf(fp, "}\n");
    pclose(fp);
  }
#endif

  
  //****************************************************************************
  //****************************************************************************
  void NetworkExpansionsAndOptimizations(
    Node<NODE_REGULAR, LINK_REGULAR>*                   node,
    ExpansionOptions        expOptions,
    STKNetworkOutputFormat  out_net_fmt,
    MyHSearchData *         wordHash,
    MyHSearchData *         nonCDphHash,
    MyHSearchData *         triphHash)
  {

    if (expOptions.mNoWordExpansion  && !expOptions.mCDPhoneExpansion &&
        expOptions.mNoOptimization   && !out_net_fmt.mNoLMLikes &&
        !out_net_fmt.mNoTimes        && !out_net_fmt.mNoWordNodes &&
        !out_net_fmt.mNoModelNodes   && !out_net_fmt.mNoPronunVars) 
    {
      return;
    }
    
    SelfLinksToNullNodes(node);

    if (!expOptions.mNoWordExpansion) {
      if (!expOptions.mNoOptimization) {
        LatticeLocalOptimization(node, expOptions.mStrictTiming, expOptions.mTraceFlag);
      }
      assert(wordHash != NULL);
      ExpandWordNetworkByDictionary(node, wordHash, !expOptions.mRemoveWordsNodes,
                                                    !expOptions.mRespectPronunVar);
    }
    if (expOptions.mCDPhoneExpansion) {
      if (!expOptions.mNoOptimization) {
        LatticeLocalOptimization(node, expOptions.mStrictTiming, expOptions.mTraceFlag);
      }
      assert(triphHash != NULL && nonCDphHash != NULL);
      ExpandMonophoneNetworkToTriphones(node, nonCDphHash, triphHash);
    }
    DiscardUnwantedInfoInNetwork(node, out_net_fmt);
  
    if (!expOptions.mNoOptimization) {
      LatticeLocalOptimization(node, expOptions.mStrictTiming, expOptions.mTraceFlag);
    }

    RemoveRedundantNullNodes(node);
  } // void NetworkExpansionsAndOptimizations( )

} // namespace STK
