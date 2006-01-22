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
#define VERSION "0.2 "__TIME__" "__DATE__
#include "Models.h"
#include "labels.h"
#include "dict.h"
#include "Net.h"
//#include "Viterbi.h"
#include "common.h"

#ifndef WIN32
#include <unistd.h>
#else
#include "getopt.h"
#endif
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <malloc.h>
#include <assert.h>
#include <ctype.h>
#include <stdarg.h>

#define SIGNIFICANT_PROB_DIFFERENCE (0.01)

namespace STK
{
  void FreeNetwork(Node *node) 
  {
    Node *tnode;
    while (node) {
      tnode = node->mpNext;
      free(node->mpLinks);
      free(node->mpBackLinks);
      free(node->mpName);
      free(node);
      node = tnode;
    }
  }
  
  Node *MakeNetworkFromLabels(Label *labels, enum NodeType node_type)
  {
    Label *lp;
    Node *first, *last = NULL, *node;
  
    if ((first           = (Node *) calloc(1, sizeof(Node))) == NULL ||
      (last            = (Node *) calloc(1, sizeof(Node))) == NULL ||
      (first->mpLinks    = (Link *) malloc(sizeof(Link))) == NULL ||
      (last->mpBackLinks = (Link *) malloc(sizeof(Link))) == NULL) {
      Error("Insufficient memory");
    }
    first->mType         = last->mType         = NT;
    first->mpPronun        = last->mpPronun        = NULL;
    first->mNLinks        = last->mNBackLinks    = 1;
    first->mNBackLinks    = last->mNLinks        = 0;
    first->mpBackLinks     = last->mpLinks         = NULL;
    first->mStart         = last->mStart         = UNDEF_TIME;
    first->mStop          = last->mStop          = UNDEF_TIME;
  //  first->tokens        = last->tokens        = NULL;
  //  first->exitToken     = last->exitToken     = NULL;
  
    node = first;
    for (lp = labels; lp != NULL; lp = lp->mpNext) {
      Node *tnode;
  
      if ((tnode            = (Node *) calloc(1, sizeof(Node))) == NULL ||
        (tnode->mpLinks     = (Link *) malloc(sizeof(Link))) == NULL ||
        (tnode->mpBackLinks = (Link *) malloc(sizeof(Link))) == NULL) {
        Error("Insufficient memory");
      }
      node->mpLinks[0].mpNode    = tnode;
      node->mpLinks[0].mLike    = 0.0;
      switch (node_type) {
        case NT:  tnode->mpPronun = ((Word *) lp->mpData)->pronuns[0]; break;
        case NT_Model: tnode->mpHmm    =   (Hmm *) lp->mpData;              break;
        case NT_Phone: tnode->mpName   =  (char *) lp->mpData;              break;
        default:       Error("Fatal: Invalid node type");
      }
      tnode->mType       = node_type;
      tnode->mNLinks     = 1;
      tnode->mNBackLinks = 1;
      tnode->mStart      = lp->mStart;
      tnode->mStop       = lp->mStop;
      tnode->mpBackLinks[0].mpNode = node;
      tnode->mpBackLinks[0].mLike = 0.0;
      node->mpNext = tnode;
      node = tnode;
    }
    node->mpNext = last;
    node->mpLinks[0].mpNode    = last;
    node->mpLinks[0].mLike    = 0.0;
    last->mpBackLinks[0].mpNode = node;
    last->mpBackLinks[0].mLike = 0.0;
    last->mpNext = NULL;
    return first;
  }
  
  
  void ExpandWordNetworkByDictionary(
    Node *first,
    struct my_hsearch_data *dict,
    int keep_word_nodes,
    int multiple_pronun)
  {
    Node *node, *prev = NULL;
    int i, j;
  
    Pronun singlePronun, *singlePronunPtr;
    Word singlePronunWrd;
    singlePronunWrd.npronuns = 1;
    singlePronunWrd.pronuns  = &singlePronunPtr;
    singlePronunPtr = &singlePronun;
  
    assert(first != NULL || first->mType & NT || first->mpPronun == NULL);
  
    for (node = first; node != NULL; prev = node, node = node->mpNext) {
  
      if (!(node->mType & NT)) continue;
  
      if (node->mpPronun == NULL) continue;
      Word *word = node->mpPronun->word;
  
      //Do not expand non-dictionary words, which where added by ReadSTKNetwork
      if (word->npronunsInDict == 0) continue;
  
      if (!multiple_pronun) {
        singlePronunWrd.mpName = node->mpPronun->word->mpName;
        word = &singlePronunWrd;
        *word->pronuns = node->mpPronun;
      }
  
      // Remove links to current node form backlinked nodes and realloc
      // link arrays of backlinked nodes to hold word->npronuns more backlinks
      for (i = 0; i < node->mNBackLinks; i++) {
        Node *bakcnode = node->mpBackLinks[i].mpNode;
        for (j=0; j<bakcnode->mNLinks && bakcnode->mpLinks[j].mpNode!=node; j++);
        assert(j < bakcnode->mNLinks); // Otherwise link to 'node' is missing
                                      // from which backlink exists
        bakcnode->mpLinks[j] = bakcnode->mpLinks[bakcnode->mNLinks-1];
  
        bakcnode->mpLinks = (Link *)
          realloc(bakcnode->mpLinks,
                (bakcnode->mNLinks - 1 + word->npronuns) * sizeof(Link));
        if (bakcnode->mpLinks == NULL) Error("Insufficient memory");
        bakcnode->mNLinks--;// += word->npronuns-1;
      }
  
      // Remove backlinks to current node form linked nodes and realloc
      // backlink arrays of linked nodes to hold word->npronuns more backlinks
      for (i=0; i < node->mNLinks; i++) {
        Node *forwnode = node->mpLinks[i].mpNode;
        for (j=0;j<forwnode->mNBackLinks&&forwnode->mpBackLinks[j].mpNode!=node;j++);
        assert(j < forwnode->mNBackLinks);
        // Otherwise link to 'node' is missing from which backlink exists
        forwnode->mpBackLinks[j] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
  
        forwnode->mpBackLinks = (Link *)
          realloc(forwnode->mpBackLinks,
                (forwnode->mNBackLinks - 1 + word->npronuns) * sizeof(Link));
        if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
        forwnode->mNBackLinks--;
      }
      for (i = 0; i < word->npronuns; i++) {
        Pronun *pronun = word->pronuns[i];
        Node *pronun_first = NULL, *pronun_prev = NULL, *tnode;
  
        for (j = 0; j < pronun->nmodels; j++) {
          tnode = (Node *) calloc(1, sizeof(Node));
          if (tnode == NULL) Error("Insufficient memory");
  
          tnode->mType       = NT_Phone | (node->mType & NT_True);
          tnode->mpName       = pronun->model[j].mpName;
          tnode->mStart      = node->mStart;
          tnode->mStop       = node->mStop;
          tnode->phoneAccuracy = 1.0;
  
          if (j == 0) {
            pronun_first = tnode;
          } else {
            if ((pronun_prev->mpLinks = (Link *) malloc(sizeof(Link))) == NULL ||
              (tnode->mpBackLinks   = (Link *) malloc(sizeof(Link))) == NULL) {
              Error("Insufficient memory");
            }
            tnode->mNBackLinks          = 1;
            tnode->mpBackLinks[0].mpNode   = pronun_prev;
            tnode->mpBackLinks[0].mLike   = 0.0;
            pronun_prev->mNLinks        = 1;
            pronun_prev->mpLinks[0].mpNode = tnode;
            pronun_prev->mpLinks[0].mLike = 0.0;
            pronun_prev->mpNext          = tnode;
          }
          pronun_prev = tnode;
        }
        if (keep_word_nodes || j == 0) {
          tnode = (Node *) calloc(1, sizeof(Node));
          if (tnode == NULL) Error("Insufficient memory");
  
          tnode->mType       = NT | (node->mType & NT_True);
          tnode->mpPronun     = keep_word_nodes ? word->pronuns[i] : NULL;
          tnode->mStart      = node->mStart;
          tnode->mStop       = node->mStop;
  
          if (j == 0) {
            pronun_first = tnode;
          } else {
            if ((pronun_prev->mpLinks = (Link *) malloc(sizeof(Link))) == NULL ||
              (tnode->mpBackLinks   = (Link *) malloc(sizeof(Link))) == NULL) {
              Error("Insufficient memory");
            }
            tnode->mNBackLinks          = 1;
            tnode->mpBackLinks[0].mpNode   = pronun_prev;
            tnode->mpBackLinks[0].mLike   = 0.0;
            pronun_prev->mNLinks        = 1;
            pronun_prev->mpLinks[0].mpNode = tnode;
            pronun_prev->mpLinks[0].mLike = 0.0;
            pronun_prev->mpNext          = tnode;
          }
          pronun_prev = tnode;
        }
        if ((pronun_prev->mpLinks =
              (Link *) malloc(sizeof(Link) * node->mNLinks))==NULL ||
          (pronun_first->mpBackLinks =
              (Link *) malloc(sizeof(Link) * node->mNBackLinks)) == NULL) {
          Error("Insufficient memory");
        }
        pronun_prev->mNLinks      = node->mNLinks;
        pronun_first->mNBackLinks = node->mNBackLinks;
  
        for (j = 0; j < node->mNBackLinks; j++) {
          Node *backnode = node->mpBackLinks[j].mpNode;
          backnode->mpLinks[backnode->mNLinks  ].mpNode = pronun_first;
          backnode->mpLinks[backnode->mNLinks++].mLike = node->mpBackLinks[j].mLike;
          pronun_first->mpBackLinks[j] = node->mpBackLinks[j];
        }
        for (j=0; j < node->mNLinks; j++) {
          Node *forwnode = node->mpLinks[j].mpNode;
          forwnode->mpBackLinks[forwnode->mNBackLinks  ].mpNode = pronun_prev;
          forwnode->mpBackLinks[forwnode->mNBackLinks++].mLike = node->mpLinks[j].mLike;
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
  
  
  Node *DiscardUnwantedInfoInNetwork(Node *first, STKNetworkOutputFormat format)
  // The function discard the information in network records that is not to be
  // saved to the output. This should allow for more effective network
  // optimization, which will be run after calling this function and before
  // saving network to file.
  {
    Node *node;
    int  i;
  
    for (node = first; node != NULL; node = node->mpNext)  {
      if (format.no_LM_likes) {
        for (i=0; i < node->mNLinks;     i++) node->mpLinks    [i].mLike = 0.0;
        for (i=0; i < node->mNBackLinks; i++) node->mpBackLinks[i].mLike = 0.0;
      }
      if (format.no_times) {
        node->mStop = node->mStart = UNDEF_TIME;
      }
      if (format.no_word_nodes && node->mType & NT) {
        node->mpPronun = NULL;
      }
      if (format.no_model_nodes && (node->mType&NT_Model || node->mType&NT_Phone)) {
        node->mType = NT;
        node->mpPronun = NULL;
      }
      if (format.no_pronun_vars && node->mType & NT && node->mpPronun != NULL) {
        node->mpPronun = node->mpPronun->word->pronuns[0];
      }
    }
    return first;
  }
  
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
  
  void WriteSTKNetwork(
    FILE                   *lfp,
    Node                   *first,
    STKNetworkOutputFormat format,
    long                   sampPeriod,
    const char             *net_file,
    const char             *out_MNF)
  {
    int n, l=0;
    Node *node;
  
    for (n = 0, node = first; node != NULL; node = node->mpNext, n++)  {
      node->mAux = n;
      l += node->mNLinks;
    }
    fprintf(lfp,"N=%d L=%d\n", n, l);
    for (node = first; node != NULL; node = node->mpNext) {
      int j;
  
      if (format.all_field_names) fputs("I=", lfp);
      if (format.base62_labels) fprintBase62(lfp, node->mAux);
      else                     fprintf(lfp,"%d", node->mAux);
  
      if (!format.no_times && node->mStop != UNDEF_TIME) {
        fputs(" t=", lfp);
  
        if (node->mStart != UNDEF_TIME && format.start_times) {
          fprintf(lfp,"%g,", node->mStart * 1.0e-7 * sampPeriod);
        }
        fprintf(  lfp,"%g",  node->mStop  * 1.0e-7 * sampPeriod);
      }
      if (!(node->mType & NT && node->mpPronun == NULL)
        || !format.no_defaults) {
        putc(' ', lfp);
        putc(node->mType & NT   ? 'W' :
            node->mType & NT_Subnet ? 'S' :
                                      'M', lfp); // NT_Model, NT_Phone
        putc('=', lfp);
        fprintHTKstr(lfp, node->mType & NT_Model   ? node->mpHmm->mpMacro->mpName   :
                          node->mType & NT    ? (!node->mpPronun ? "!NULL" :
                                                    node->mpPronun->word->mpName) :
                                                    node->mpName); // NT_PHONE (NT_Subnet)
      }
      if (!format.no_pronun_vars && node->mType & NT
      && node->mpPronun != NULL && node->mpPronun->word->npronuns > 1
      && (node->mpPronun->variant_no > 1 || !format.no_defaults)) {
        fprintf(lfp," v=%d", node->mpPronun->variant_no);
      }
      if (node->mType & NT_True || node->mType & NT_Sticky) {
        fputs(" f=", lfp);
        if (node->mType & NT_True)   putc('T', lfp);
        if (node->mType & NT_Sticky) putc('K', lfp);
      }
      if (node->mType & NT_Phone && node->phoneAccuracy != 1.0) {
        fprintf(lfp," p="FLOAT_FMT, node->phoneAccuracy);
      }
      if (!format.arc_defs_to_end) {
        if (format.all_field_names) fprintf(lfp," J=%d", node->mNLinks);
  
        for (j = 0; j < node->mNLinks; j ++) {
          putc(' ', lfp);
          if (format.all_field_names) fputs("E=", lfp);
          if (format.base62_labels) fprintBase62(lfp, node->mpLinks[j].mpNode->mAux);
          else                     fprintf(lfp,"%d", node->mpLinks[j].mpNode->mAux);
          if (node->mpLinks[j].mLike != 0.0 && !format.no_LM_likes) {
            fprintf(lfp," l="FLOAT_FMT, node->mpLinks[j].mLike);
          }
        }
      }
      fputs("\n", lfp);
      if (ferror(lfp)) {
        Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
      }
    }
  
    if (format.arc_defs_to_end) {
      l = 0;
      for (node = first; node != NULL; node = node->mpNext) {
        int j;
  
        for (j = 0; j < node->mNLinks; j ++) {
          if (format.all_field_names) {
            fprintf(lfp, format.arc_defs_with_J ? "J=%d S=" : "I=", l++);
          }
          if (format.base62_labels) fprintBase62(lfp, node->mAux);
          else                     fprintf(lfp,"%d", node->mAux);
          putc(' ', lfp); // space = ' ';
          if (format.all_field_names) fputs("E=", lfp);
          if (format.base62_labels) fprintBase62(lfp, node->mpLinks[j].mpNode->mAux);
          else                     fprintf(lfp,"%d", node->mpLinks[j].mpNode->mAux);
          if (node->mpLinks[j].mLike != 0.0 && !format.no_LM_likes) {
            fprintf(lfp," l="FLOAT_FMT, node->mpLinks[j].mLike);
          }
          fputs("\n", lfp);
          if (ferror(lfp)) {
            Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
          }
        }
      }
    }
  }
  
  
  static int lnkcmp(const void *a, const void *b)
  {
  //  return ((Link *) a)->mpNode - ((Link *) b)->mpNode;
  //  Did not work with gcc, probably bug in gcc pointer arithmetic
    return (char *)((Link *) a)->mpNode - (char *)((Link *) b)->mpNode;
  }
  
  
  static int LatticeLocalOptimization_ForwardPass(Node *first, int strictTiming)
  {
    int i, j, k, l, m, rep;
    Node *node, *tnode;
    int node_removed = 0;
    FLOAT tlike;
    for (node = first; node != NULL; node = node->mpNext) {
      for (i = 0; i < node->mNLinks; i++) {
      //for (tnode = inode; tnode != NULL; tnode = (tnode == inode ? jnode : NULL)) {
        tnode = node->mpLinks[i].mpNode;
        if (tnode->mNLinks == 0) continue;
  
        // Weight pushing
        tlike = tnode->mpBackLinks[0].mLike;
        for (l=1; l <  tnode->mNBackLinks; l++) {
          tlike = HIGHER_OF(tlike, tnode->mpBackLinks[l].mLike);
        }
        for (l=0; l < tnode->mNBackLinks; l++) {
          Node *backnode = tnode->mpBackLinks[l].mpNode;
          tnode->mpBackLinks[l].mLike -= tlike;
          for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].mpNode!=tnode; k++);
          assert(k < backnode->mNLinks);
          backnode->mpLinks[k].mLike -= tlike;
      #ifndef NDEBUG
          for (k++; k<backnode->mNLinks && backnode->mpLinks[k].mpNode!=tnode; k++);
      #endif
          assert(k == backnode->mNLinks);
        }
        for (l=0; l < tnode->mNLinks; l++) {
          Node *forwnode = tnode->mpLinks[l].mpNode;
          tnode->mpLinks[l].mLike += tlike;
          for (k=0; k<forwnode->mNBackLinks && forwnode->mpBackLinks[k].mpNode!=tnode;k++);
          assert(k < forwnode->mNBackLinks);
          forwnode->mpBackLinks[k].mLike += tlike;
      #ifndef NDEBUG
          for (k++; k<forwnode->mNBackLinks && forwnode->mpBackLinks[k].mpNode!=tnode;k++);
      #endif
          assert(k == forwnode->mNBackLinks);
        }
      }
  //dnet(first, 1, node);
  
// For current node 'node', check for each possible pair of its successors 
// ('inode' and 'jnode') whether the pair may be merged to single node.
      for (i = 0; i < node->mNLinks-1; i++) {
        for (j = i+1; j < node->mNLinks; j++) {
          Node *inode = node->mpLinks[i].mpNode;
          Node *jnode = node->mpLinks[j].mpNode;

// Final node may be never merged.
        if (inode->nlinks == 0 || jnode->nlinks == 0) continue;
  
// Two nodes ('inode' and 'jnode') may be mergeg if they are of the same type, name, ...
// with the same predecessors and with the same weights on the links from predecesors.
          if ((inode->mType & ~NT_True) != (jnode->mType & ~NT_True)
          || ( inode->mType & NT_Phone && inode->mpName   != jnode->mpName)
          || ( inode->mType & NT  && inode->mpPronun != jnode->mpPronun)
  //          &&  (inode->mpPronun == NULL ||
  //             jnode->mpPronun == NULL ||
  //             inode->mpPronun->word       != jnode->mpPronun->word ||
  //             inode->mpPronun->outSymbol  != jnode->mpPronun->outSymbol ||
  //             inode->mpPronun->variant_no != jnode->mpPronun->variant_no ||
  //             inode->mpPronun->prob       != jnode->mpPronun->prob)
          || (inode->mNBackLinks != jnode->mNBackLinks)) {
          continue;
          }
          if (strictTiming && (inode->mStart != jnode->mStart
                          ||  inode->mStop  != jnode->mStop)) {
            continue;
          }

//Weights on the links from predecesors does not have to be exactely the same, but the must not
//differ more than by SIGNIFICANT_PROB_DIFFERENCE
          for (l=0; l < inode->mNBackLinks; l++) {
            if (inode->mpBackLinks[l].mpNode != jnode->mpBackLinks[l].mpNode) break;
            FLOAT ldiff =  inode->mpBackLinks[l].mLike - jnode->mpBackLinks[l].mLike;
            if (ldiff < -SIGNIFICANT_PROB_DIFFERENCE ||
              ldiff >  SIGNIFICANT_PROB_DIFFERENCE) break;
          }
          if (l < inode->mNBackLinks) continue;
  
  /*        if (memcmp(inode->mpBackLinks, jnode->mpBackLinks,
                    inode->mNBackLinks * sizeof(Link))) {
            continue;
          }*/
            // inode and jnode are the same nodes with the same predeccessors
            // Remove jnode and add its links to inode
  
          assert(inode->mNLinks && jnode->mNLinks);
  
            // Remove links to jnode form predeccessors
          for (l=0; l < jnode->mNBackLinks; l++) {
            Node *backnode = jnode->mpBackLinks[l].mpNode;
            for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].mpNode!=jnode; k++);
            assert(k < backnode->mNLinks);
            // Otherwise link to 'node' is missing from which backlink exists
            memmove(backnode->mpLinks+k, backnode->mpLinks+k+1,
                    (backnode->mNLinks-k-1) * sizeof(Link));
            backnode->mNLinks--;
          }
  
          // Merge jnode's links with inode links
  
          //Count jnode's links not present among inode's links
          rep = l = k = 0;
          while (k < jnode->mNLinks) {
            Link *ill = inode->mpLinks+l;
            Link *jlk = jnode->mpLinks+k;
            if (l == inode->mNLinks || ill->mpNode > jlk->mpNode){
              // k-th link of jnode will be included among inode's links.
              // Redirect corresponding baclink to inode
              for (m = 0; m < jlk->mpNode->mNBackLinks
                          && jlk->mpNode->mpBackLinks[m].mpNode != jnode; m++);
              assert(m < jlk->mpNode->mNBackLinks);
              jlk->mpNode->mpBackLinks[m].mpNode = inode;
              qsort(jlk->mpNode->mpBackLinks, jlk->mpNode->mNBackLinks,
                    sizeof(Link), lnkcmp);
              k++;
            } else  if (ill->mpNode == jlk->mpNode) {
              // l-th link of inode and k-th link of jnode points to
              // the same node. Link from jnode is redundant.
              // Remove backlinks to jnode form jnode's succesors
              for (m = 0; m < jlk->mpNode->mNBackLinks
                        && jlk->mpNode->mpBackLinks[m].mpNode != jnode; m++);
              assert(m < jlk->mpNode->mNBackLinks);
              memmove(jlk->mpNode->mpBackLinks+m, jlk->mpNode->mpBackLinks+m+1,
                      (jlk->mpNode->mNBackLinks-m-1) * sizeof(Link));
              jlk->mpNode->mNBackLinks--;
  
              ill->mLike = HIGHER_OF(ill->mLike, jlk->mLike);
              jlk->mpNode = NULL; // Mark link to be removed
              rep++; k++, l++;
            } else {
              l++;
            }
          }
          l = inode->mNLinks;
          inode->mNLinks += jnode->mNLinks-rep;
          inode->mpLinks = (Link *) realloc(inode->mpLinks,
                                          inode->mNLinks * sizeof(Link));
          if (inode->mpLinks == NULL) Error("Insufficient memory");
  
          for (k = 0; k < jnode->mNLinks; k++) {
            if (jnode->mpLinks[k].mpNode != NULL) {
              inode->mpLinks[l++] = jnode->mpLinks[k];
            }
          }
          qsort(inode->mpLinks, inode->mNLinks, sizeof(Link), lnkcmp);
  
          inode->mStart = inode->mStart == UNDEF_TIME || jnode->mStart == UNDEF_TIME
                        ? UNDEF_TIME : LOWER_OF(inode->mStart, jnode->mStart);
  
          inode->mStop  = inode->mStop == UNDEF_TIME || jnode->mStop == UNDEF_TIME
                        ? UNDEF_TIME : HIGHER_OF(inode->mStop, jnode->mStop);
  
          if (inode->mAux > jnode->mAux) {
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
          } else {
            jnode->mpNext->mpBackNext = jnode->mpBackNext;
            jnode->mpBackNext->mpNext = jnode->mpNext;
          }
          inode->mType |= jnode->mType & NT_True;
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
  
  Node *ReverseNetwork(Node *first)
  {
    Node *node, *last = NULL;
    for (node = first; node != NULL; node = node->mpBackNext) {
      Link *links  = node->mpLinks;
      int  nlinks  = node->mNLinks;
      node->mpLinks  = node->mpBackLinks;
      node->mNLinks = node->mNBackLinks;
      node->mpBackLinks  = links;
      node->mNBackLinks = nlinks;
      Node *next     = node->mpNext;
      node->mpNext     = node->mpBackNext;
      node->mpBackNext = next;
      node->mAux      = -node->mAux;
      last = node;
    }
    return last;
  }
  
  static int LatticeLocalOptimization_BackwardPass(Node *first, int strictTiming)
  {
    int node_removed;
    Node *last = ReverseNetwork(first);
    node_removed = LatticeLocalOptimization_ForwardPass(last, strictTiming);
  //  if (!node_removed) dnet(last, 0);
    ReverseNetwork(last);
  //  if (!node_removed) dnet(first, 0);
    return node_removed;
  }
  
  
  int nbacklinkscmp(const void *a, const void *b) {
    return (*(Node **) a)->mNBackLinks - (*(Node **) b)->mNBackLinks;
  }
  
  void LatticeLocalOptimization(Node *first, int strictTiming, int trace_flag)
  {
  //  Node *last, *node, *unsorted;
  //  int topolOrd = 1;
  
    Node *node, *lastnode;
    int i, j, unreachable = 0;
    // For each node, sort links by pointer value to allow
    // for easy comparison whether two nodes have the same set of links
  
    for (node = first; node != NULL; node = node->mpNext)  {
      node->mAux = 0;
      node->mpBackNext = node->mpNext;
      qsort(node->mpLinks, node->mNLinks, sizeof(Link), lnkcmp);
      qsort(node->mpBackLinks, node->mNBackLinks, sizeof(Link), lnkcmp);
    }
  
  // Sort nodes in topological order
  // printf("Sorting nodes...\n");
  
    first->mAux = 1;
    for (lastnode = node = first; node != NULL; node = node->mpNext) {
      for (i=0; i < node->mNLinks; i++) {
        Node *lnknode = node->mpLinks[i].mpNode;
        if (lnknode->mAux == 0) {
          for (j=0; j<lnknode->mNBackLinks && lnknode->mpBackLinks[j].mpNode->mAux==1; j++);
          if (j == lnknode->mNBackLinks) {
            lastnode->mpNext = lnknode;
            lastnode  = lnknode;
            lnknode->mAux = 1;
            lnknode->mpNext = NULL;
          }
        }
      }
    }
    if (lastnode->mNLinks != 0) {
      // There is a cycle in graph so we cannot sort nodes
      // topologicaly, so sort it at least somehow. :o|
      // Anyway this optimization algorithm is not optimal for graphs with cycles.
      for (node = first; node != NULL; node = node->mpBackNext) node->mAux = 0;
      first->mAux = 1;
      for (lastnode = node = first; node != NULL; node = node->mpNext) {
        for (i=0; i < node->mNLinks; i++) {
          Node *lnknode = node->mpLinks[i].mpNode;
          if (lnknode->mAux == 0) {
            lastnode->mpNext = lnknode;
            lastnode  = lnknode;
            lnknode->mAux = 1;
            lnknode->mpNext = NULL;
          }
        }
      }
      for (node=first; node->mpNext->mNLinks != 0; node=node->mpNext);
  
      if (node->mpNext->mpNext) { // Final node is not at the and of chain
        lastnode->mpNext = node->mpNext;
        node->mpNext = node->mpNext->mpNext;
        lastnode = lastnode->mpNext;
        lastnode->mpNext = NULL;
      }
    }
  
    // !!! Unreachable nodes must be removed before sorting !!!
  
    for (node=first; node != NULL; node=node->mpBackNext) {
      while (node->mpBackNext && node->mpBackNext->mAux == 0) {
        Node *tnode = node->mpBackNext;
        node->mpBackNext = node->mpBackNext->mpBackNext;
        unreachable++;
        free(tnode->mpLinks);
        free(tnode->mpBackLinks);
        free(tnode);
      }
    }
  
  //  if (unreachable) Warning("Removing %d unreachable nodes", unreachable);
    if (unreachable) Error("Networks contains unreachable nodes");
  
    first->mpBackNext = NULL;
    for (node=first; node->mpNext != NULL; node=node->mpNext) {
      node->mpNext->mpBackNext = node;
    }
    for (i=1, node=first; node != NULL; node = node->mpNext, i++) {
      node->mAux=i;
    }
    for (;;) {
      if (trace_flag & 2) {
        for (i=0,node=first; node; node=node->mpNext,i++);
        TraceLog("Forward pass.... (number of nodes: %d)", i);
      }
      LatticeLocalOptimization_ForwardPass(first, strictTiming);
  
      if (trace_flag & 2) {
        for (i=0,node=first; node; node=node->mpNext,i++);
        TraceLog("Backward pass... (number of nodes: %d)", i);
      }
      if (!LatticeLocalOptimization_BackwardPass(first, strictTiming)) break;
    }
  }
  
  void ExpandMonophoneNetworkToTriphones(
    Node *first,
    struct my_hsearch_data *nonCDphones,
    struct my_hsearch_data *CDphones)
  {
    Node *node;
    int did_we_clone, i, j, k;
    // Find all Tee model, Word, and Null nodes (except the first and last Null node)
    // and clone those not having single input and output
  
    do {
      ENTRY e, *ep;
      Node *prev = NULL;
      did_we_clone = 0;
      for (node = first; node != NULL; prev = node, node = node->mpNext) {
        if (node->mNLinks == 0 || node->mNBackLinks == 0 ||
          (node->mNLinks == 1 && node->mNBackLinks == 1)) {
          continue;
        }
        if (node->mType & NT_Phone) {
          e.key = node->mpName;
          my_hsearch_r(e, FIND, &ep, nonCDphones);
          if (ep == NULL || !(int) ep->data) continue; // Node is not a Tee model
        }
        did_we_clone = 1;
        assert(prev != NULL); //Otherwise first node is not Null node
  
        // Remove links to current node form back-linked nodes and realloc
        // link arrays of back-linked nodes to hold node->mNLinks more links
        for (j=0; j < node->mNBackLinks; j++) {
          Node *backnode = node->mpBackLinks[j].mpNode;
          for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].mpNode!=node; k++);
          assert(k < backnode->mNLinks);
          // Otherwise link to 'node' is missing from which backlink exists
          backnode->mpLinks[k] = backnode->mpLinks[backnode->mNLinks-1];
  
          backnode->mpLinks = (Link *)
            realloc(backnode->mpLinks,
                  (backnode->mNLinks-1+node->mNLinks)*sizeof(Link));
          if (backnode->mpLinks == NULL) Error("Insufficient memory");
          backnode->mNLinks--;
        }
        // Remove backlinks to current node form linked nodes and realloc
        // backlink arrays of linked nodes to hold node->mNBackLinks more backlinks
        for (j=0; j < node->mNLinks; j++) {
          Node *forwnode = node->mpLinks[j].mpNode;
          for (k=0;k<forwnode->mNBackLinks&&forwnode->mpBackLinks[k].mpNode!=node;k++);
          assert(k < forwnode->mNBackLinks);
          // Otherwise link to 'node' is missing from which backlink exists
          forwnode->mpBackLinks[k] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
  
          forwnode->mpBackLinks = (Link *)
            realloc(forwnode->mpBackLinks,
                  (forwnode->mNBackLinks-1+node->mNBackLinks)*sizeof(Link));
          if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
          forwnode->mNBackLinks--;
        }
        // Alloc new node->mNLinks * node->mNBackLinks nodes and create new links
        // so that each backlinked node is conected with each linked node through
        // one new node.
        for (i=0; i < node->mNLinks; i++) 
        {
          for (j=0; j < node->mNBackLinks; j++) 
          {
            Node *tnode;
            Link forwlink = node->mpLinks[i];
            Link backlink = node->mpBackLinks[j];
  
            if ((tnode           = (Node *) calloc(1, sizeof(Node))) == NULL)
              Error("Insufficient memory");
            
            *tnode = *node;
  
            if ((tnode->mpLinks     = (Link *) malloc(sizeof(Link))) == NULL ||
              (tnode->mpBackLinks = (Link *) malloc(sizeof(Link))) == NULL) 
            {
              Error("Insufficient memory");
            }
            
            tnode->mNLinks       = 1;
            tnode->mNBackLinks   = 1;
            tnode->mpLinks[0]     = forwlink;
            tnode->mpBackLinks[0] = backlink;
            
            forwlink.mpNode->mpBackLinks[forwlink.mpNode->mNBackLinks  ].mpNode = tnode;
            forwlink.mpNode->mpBackLinks[forwlink.mpNode->mNBackLinks++].mLike  = forwlink.mLike;
            backlink.mpNode->mpLinks    [backlink.mpNode->mNLinks      ].mpNode = tnode;
            backlink.mpNode->mpLinks    [backlink.mpNode->mNLinks++    ].mLike  = backlink.mLike;
            
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
    int nbackmononodes, nforwmononodes, id = 0;
    for (node = first; node != NULL; node = node->mpNext) node->mAux = id++;
  
    // Expand monophone nodes to triphone nodes
    Node *prev = NULL;
    for (node = first; node != NULL; prev = node, node = node->mpNext) {
      ENTRY e, *ep;
  
      if (node->mType & NT ||
        (node->mNLinks == 1 && node->mNBackLinks == 1)) {
        continue;
      }
      assert(node->mType & NT_Phone);
      e.key = node->mpName;
      my_hsearch_r(e, FIND, &ep, nonCDphones);
      if (ep != NULL && (int) ep->data) continue; // Node is a Tee model
  
      assert(prev != NULL); //Otherwise first node is not Null node
  
      // Count groups of backlinked nodes corresponding to different monophones
      id = -1;
      nbackmononodes = 0;
      for (j=0; j < node->mNBackLinks; j++) {
        if (node->mpBackLinks[j].mpNode->mAux != id) {
          id = node->mpBackLinks[j].mpNode->mAux;
          nbackmononodes++;
        }
      }
      // Count groups of linked nodes corresponding to different monophones
      id = -1;
      nforwmononodes = 0;
      for (j=0; j < node->mNLinks; j++) {
        if (node->mpLinks[j].mpNode->mAux != id) {
          id = node->mpLinks[j].mpNode->mAux;
          nforwmononodes++;
        }
      }
  
      // Remove links to current node form backlinked nodes and realloc
      // link arrays of backlinked nodes to hold nforwmononodes more links
      for (j=0; j < node->mNBackLinks; j++) {
        Node *backnode = node->mpBackLinks[j].mpNode;
        for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].mpNode!=node; k++);
        assert(k < backnode->mNLinks);
        // Otherwise link to 'node' is missing from which backlink exists
        memmove(backnode->mpLinks+k, backnode->mpLinks+k+1,
                (backnode->mNLinks-k-1) * sizeof(Link));
  
        backnode->mpLinks = (Link *)
          realloc(backnode->mpLinks,
                (backnode->mNLinks-1+nforwmononodes)*sizeof(Link));
        if (backnode->mpLinks == NULL) Error("Insufficient memory");
        backnode->mNLinks--;
      }
  
      // Remove backlinks to current node form linked nodes and realloc
      // backlink arrays of linked nodes to hold nbackmononodes more backlinks
      for (j=0; j < node->mNLinks; j++) {
        Node *forwnode = node->mpLinks[j].mpNode;
        for (k=0;k<forwnode->mNBackLinks&&forwnode->mpBackLinks[k].mpNode!=node;k++);
        assert(k < forwnode->mNBackLinks);
        // Otherwise link to 'node' is missing from which backlink exists
        memmove(forwnode->mpBackLinks+k, forwnode->mpBackLinks+k+1,
                (forwnode->mNBackLinks-k-1) * sizeof(Link));
  
        forwnode->mpBackLinks = (Link *)
          realloc(forwnode->mpBackLinks,
                (forwnode->mNBackLinks-1+nbackmononodes)*sizeof(Link));
        if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
        forwnode->mNBackLinks--;
      }
  
      // Alloc new nforwmononodes * nbackmononodes nodes and create new links
      // so that each backlinked node is conected through one new node with all
      // linked nodes belonging to one monophone group and vice versa each
      // linked node is conected through one new node with all backlinked nodes
      // belonging to one monophone group
      Link *forwmono_start, *forwmono_end = node->mpLinks;
      for (i=0; i < nforwmononodes; i++) {
        for (forwmono_start = forwmono_end;
            forwmono_end < node->mpLinks+node->mNLinks &&
            forwmono_start->mpNode->mAux == forwmono_end->mpNode->mAux;
            forwmono_end++);
  
        assert((i <  nforwmononodes-1 && forwmono_end <  node->mpLinks+node->mNLinks) ||
              (i == nforwmononodes-1 && forwmono_end == node->mpLinks+node->mNLinks));
  
        Link *tlink, *backmono_start, *backmono_end = node->mpBackLinks;
        for (j=0; j < nbackmononodes; j++) {
          for (backmono_start = backmono_end;
            backmono_end < node->mpBackLinks+node->mNBackLinks &&
            backmono_start->mpNode->mAux == backmono_end->mpNode->mAux;
            backmono_end++);
  
            assert((j <  nbackmononodes-1 && backmono_end <  node->mpBackLinks+node->mNBackLinks) ||
                  (j == nbackmononodes-1 && backmono_end == node->mpBackLinks+node->mNBackLinks));
  
          Node *tnode;
          if ((tnode           = (Node *) calloc(1, sizeof(Node))) == NULL) {
            Error("Insufficient memory");
          }
          *tnode = *node;
          tnode->mNLinks       = forwmono_end-forwmono_start;
          tnode->mNBackLinks   = backmono_end-backmono_start;
  
          if ((tnode->mpLinks =
              (Link *) malloc(tnode->mNLinks * sizeof(Link))) == NULL ||
            (tnode->mpBackLinks =
              (Link *) malloc(tnode->mNBackLinks * sizeof(Link))) == NULL) {
            Error("Insufficient memory");
          }
          for (tlink = forwmono_start; tlink < forwmono_end; tlink++) {
            tnode->mpLinks[tlink-forwmono_start] = *tlink;
            tlink->mpNode->mpBackLinks[tlink->mpNode->mNBackLinks  ].mpNode = tnode;
            tlink->mpNode->mpBackLinks[tlink->mpNode->mNBackLinks++].mLike = tlink->mLike;
          }
          for (tlink = backmono_start; tlink < backmono_end; tlink++) {
            tnode->mpBackLinks[tlink-backmono_start] = *tlink;
            tlink->mpNode->mpLinks[tlink->mpNode->mNLinks  ].mpNode = tnode;
            tlink->mpNode->mpLinks[tlink->mpNode->mNLinks++].mLike = tlink->mLike;
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
    for (node = first; node != NULL; node = node->mpNext) {
      ENTRY e, *ep;
      Node *lc, *rc;
      char *lcname, *rcname, *triname;
      int   lcnlen,  rcnlen;
  
      if (!(node->mType & NT_Phone)) continue;
  
      if (nonCDphones) {
        e.key  = node->mpName;
        my_hsearch_r(e, FIND, &ep, nonCDphones);
      } else {
        ep = NULL;
      }
      if (ep != NULL) {
        lc = rc = NULL;
      } else {
        for (lc = node;;) {
          lc = lc->mNBackLinks ? lc->mpBackLinks[0].mpNode : NULL;
          if (lc == NULL)           break;
          if (!(lc->mType & NT_Phone)) continue;
          if (nonCDphones == NULL)  break;
          e.key  = lc->mpName;
          my_hsearch_r(e, FIND, &ep, nonCDphones);
          if (ep == NULL || !(int) ep->data) break; // Node represents Tee model
        }
        for (rc = node;;) {
          rc = rc->mNLinks ? rc->mpLinks[0].mpNode : NULL;
          if (rc == NULL)           break;
          if (!(rc->mType & NT_Phone)) continue;
          if (nonCDphones == NULL)  break;
          e.key  = rc->mpName;
          my_hsearch_r(e, FIND, &ep, nonCDphones);
          if (ep == NULL || !(int) ep->data) break; // Node represents Tee model
        }
      }
      lcnlen = -1;
      if (lc != NULL) {
        lcname = strrchr(lc->mpName, '-');
        if (lcname == NULL) lcname = lc->mpName;
        else lcname++;
        lcnlen = strcspn(lcname, "+");
      }
      rcnlen = -1;
      if (rc != NULL) {
        rcname = strrchr(rc->mpName, '-');
        if (rcname == NULL) rcname = rc->mpName;
        else rcname++;
        rcnlen = strcspn(rcname, "+");
      }
      triname = (char *) malloc(lcnlen+1+strlen(node->mpName)+1+rcnlen+1);
      if (triname == NULL) Error("Insufficient memory");
  
      triname[0] = '\0';
  
      if (lcnlen > 0) strcat(strncat(triname, lcname, lcnlen), "-");
      strcat(triname, node->mpName);
      if (rcnlen > 0) strncat(strcat(triname, "+"), rcname, rcnlen);
  
      e.key  = triname;
      my_hsearch_r(e, FIND, &ep, CDphones);
  
      if (ep == NULL) {
        e.key  = triname;
        e.data = e.key;
  
        if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, CDphones)) {
          Error("Insufficient memory");
        }
        node->mpName = triname;
      } else {
        free(triname);
        node->mpName = ep->key;
      }
    }
  }
  
  static int getInteger(char *str, char **endPtr,
                        const char *file_name, int line_no)
  {
    long l = strtoul(str, endPtr, 10);
  
    if (str == *endPtr || (**endPtr && !isspace(**endPtr))) {
      Error("Invalid integral value (%s:%d)", file_name, line_no);
    }
    while (isspace(**endPtr)) ++*endPtr;
  
    return l;
  }
  
  static float getFloat(char *str, char **endPtr,
                        const char *file_name, int line_no)
  {
    double d = strtod(str, endPtr);
  
    if (str == *endPtr || (**endPtr && !isspace(**endPtr))) {
      Error("Invalid float value (%s:%d)", file_name, line_no);
    }
    while (isspace(**endPtr)) ++*endPtr;
  
    return d;
  }
  
  
  static int getNodeNumber(int nnodes, char *str, char **endPtr,
                          const char *file_name,int line_no)
  {
    long node_id = getInteger(str, endPtr, file_name, line_no);
  
    if (node_id < 0 || node_id >= nnodes) {
      Error("Node number out of range (%s:%d)", file_name, line_no);
    }
    return node_id;
  }
  
  // Insert null node to self links. Self links are not supported
  // in network manipulation functions
  
  void SelfLinksToNullNodes(Node *first)
  {
    int i, j;
    Node *node, *tnode;
  
    for (node = first; node != NULL; node = node->mpNext) {
      for (i=0; i < node->mNLinks; i++) {
        if (node->mpLinks[i].mpNode == node) {
          if ((tnode           = (Node *) calloc(1, sizeof(Node))) == NULL ||
            (tnode->mpLinks     = (Link *) malloc(sizeof(Link))) == NULL ||
            (tnode->mpBackLinks = (Link *) malloc(sizeof(Link))) == NULL) {
            Error("Insufficient memory");
          }
  
          node->mpLinks[i].mpNode = tnode;
          for (j=0; j<node->mNBackLinks && node->mpBackLinks[j].mpNode!=node; j++);
          assert(j<node->mNBackLinks);
          node->mpBackLinks[j].mpNode = tnode;
          node->mpBackLinks[j].mLike = 0.0;
  
          tnode->mType       = NT;
          tnode->mpPronun     = NULL;
          tnode->mNLinks     = 1;
          tnode->mNBackLinks = 1;
          tnode->mStart      = UNDEF_TIME;
          tnode->mStop       = UNDEF_TIME;
  //        tnode->tokens     = NULL;
  //        tnode->exitToken  = NULL;
          tnode->mpLinks[0].mpNode     = node;
          tnode->mpLinks[0].mLike     = 0.0;
          tnode->mpBackLinks[0].mpNode = node;
          tnode->mpBackLinks[0].mLike = node->mpLinks[i].mLike;
          tnode->mpNext = node->mpNext;
          node->mpNext = tnode;
        }
      }
    }
  }
  
  // Remove null nones having less than three predecessors or less than three successors
  int RemoveRedundantNullNodes(Node *first)
  {
    Node *    node;
    Node *    tnode;
    int       i;
    int       j;
    int       k;
    int       node_removed = 0;
  
    first->mpBackNext = NULL;
    for (node = first; node->mpNext != NULL; node = node->mpNext) {
      node->mpNext->mpBackNext = node;
    }
    for (node = first; node != NULL; node = node->mpNext) {
      if (node->mType & NT && node->mpPronun == NULL &&
          node->mNLinks != 0 && node->mNBackLinks != 0  &&
        (node->mNLinks == 1 || node->mNBackLinks == 1 ||
        (node->mNLinks == 2 && node->mNBackLinks == 2))) {
  
      node_removed = 1;
  
      // Remove links to current node form backlinked nodes and realloc
      // link arrays of backlinked nodes to hold node->mNLinks more backlinks
        for (i = 0; i < node->mNBackLinks; i++) {
          Node *bakcnode = node->mpBackLinks[i].mpNode;
          for (j=0; j<bakcnode->mNLinks && bakcnode->mpLinks[j].mpNode!=node; j++);
          assert(j < bakcnode->mNLinks); // Otherwise link to 'node' is missing
                                        // from which backlink exists
          bakcnode->mpLinks[j] = bakcnode->mpLinks[bakcnode->mNLinks-1];
  
          bakcnode->mpLinks = (Link *)
            realloc(bakcnode->mpLinks,
                  (bakcnode->mNLinks - 1 + node->mNLinks) * sizeof(Link));
          if (bakcnode->mpLinks == NULL) Error("Insufficient memory");
          bakcnode->mNLinks--;// += word->npronuns-1;
        }
  
        // Remove backlinks to current node form linked nodes and realloc
        // backlink arrays of linked nodes to hold word->npronuns more backlinks
        for (i=0; i < node->mNLinks; i++) {
          Node *forwnode = node->mpLinks[i].mpNode;
          for (j=0;j<forwnode->mNBackLinks&&forwnode->mpBackLinks[j].mpNode!=node;j++);
          assert(j < forwnode->mNBackLinks);
          // Otherwise link to 'node' is missing from which backlink exists
          forwnode->mpBackLinks[j] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
  
          forwnode->mpBackLinks = (Link *)
            realloc(forwnode->mpBackLinks,
                  (forwnode->mNBackLinks - 1 + node->mNBackLinks) * sizeof(Link));
          if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
          forwnode->mNBackLinks--;
        }
        for (j = 0; j < node->mNBackLinks; j++) {
          Node *backnode = node->mpBackLinks[j].mpNode;
          int orig_nlinks = backnode->mNLinks;
  
          for (i=0; i < node->mNLinks; i++) {
            for(k = 0; k < orig_nlinks && backnode->mpLinks[k].mpNode != node->mpLinks[i].mpNode; k++);
            if(k < orig_nlinks) {
              // Link which is to be created already exists. Its duplication must be avoided.
              backnode->mpLinks[k].mLike = HIGHER_OF(backnode->mpLinks[k].mLike, 
                                                  node->mpLinks[i].mLike + node->mpBackLinks[j].mLike);
            } else {
              backnode->mpLinks[backnode->mNLinks  ].mpNode = node->mpLinks[i].mpNode;
              backnode->mpLinks[backnode->mNLinks++].mLike
                = node->mpLinks[i].mLike + node->mpBackLinks[j].mLike;
            }
          }
        }
        for (j = 0; j < node->mNLinks; j++) {
          Node *forwnode = node->mpLinks[j].mpNode;
          int orig_nbacklinks = forwnode->mNBackLinks;
  
          for (i=0; i < node->mNBackLinks; i++) {
            for(k = 0; k < orig_nbacklinks && forwnode->mpBackLinks[k].mpNode != node->mpBackLinks[i].mpNode; k++);
            if (k < orig_nbacklinks) {
              // Link which is to be created already exists. Its duplication must be avoided.
              forwnode->mpBackLinks[k].mLike = HIGHER_OF(forwnode->mpBackLinks[k].mLike, 
                                                      node->mpBackLinks[i].mLike + node->mpLinks[j].mLike);
            } else {
              forwnode->mpBackLinks[forwnode->mNBackLinks  ].mpNode = node->mpBackLinks[i].mpNode;
              forwnode->mpBackLinks[forwnode->mNBackLinks++].mLike
                = node->mpBackLinks[i].mLike + node->mpLinks[j].mLike;
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
    Node      *   mpNode;
    long long     maxStopTimeTillNow;
    int           mId;
  };
  
  static int cmp_starts(const void *a, const void *b)
  {
    long long diff = ((struct CorrPhnRec *) a)->mpNode->mStart
                  - ((struct CorrPhnRec *) b)->mpNode->mStart;
  
    if (diff != 0)
      return diff;
  
    return ((struct CorrPhnRec *) a)->mpNode->mStop
        - ((struct CorrPhnRec *) b)->mpNode->mStop;
  }
  
  static int cmp_maxstop(const void *key, const void *elem)
  {
    struct CorrPhnRec *corr_phn = (struct CorrPhnRec *) elem;
  
    if (((Node *) key)->mStart < corr_phn->maxStopTimeTillNow) {
      if (corr_phn->mId == 0 || // first fiead in the array
        ((Node *) key)->mStart >= (corr_phn-1)->maxStopTimeTillNow)
          return  0;
      else return -1;
    } else return  1;
  }
  
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
  
  
  void ComputeAproximatePhoneAccuracy(Node *first, int type)
  {
    Node *node;
    struct CorrPhnRec *corr_phn, *overlaped;
    int ncorr_phns = 0, i = 0;
    long long maxStopTime;
  
    for (node = first; node != NULL; node = node->mpNext) {
      if (node->mType & NT_Phone && node->mType & NT_True) ncorr_phns++;
    }
    if (ncorr_phns == 0) Error("No correct phoneme node in network");
  
    corr_phn = (struct CorrPhnRec *)malloc(sizeof(struct CorrPhnRec)*ncorr_phns);
    if (corr_phn == NULL) Error("Insufficient memory");
  
    for (node = first; node != NULL; node = node->mpNext) {
      if (node->mType & NT_Phone && node->mType & NT_True) {
        corr_phn[i++].mpNode = node;
      }
    }
  
    qsort(corr_phn, ncorr_phns, sizeof(struct CorrPhnRec), cmp_starts);
  
    maxStopTime = corr_phn[0].mpNode->mStop;
    for (i = 0; i < ncorr_phns; i++) {
      corr_phn[i].mId = i;
      maxStopTime = HIGHER_OF(maxStopTime, corr_phn[i].mpNode->mStop);
      corr_phn[i].maxStopTimeTillNow = maxStopTime;
    }
    for (node = first; node != NULL; node = node->mpNext) {
      if (!(node->mType & NT_Phone)) continue;
  
      if (node->mStop  <= node->mStart ||
        !strcmp(node->mpName, "sil") ||
        !strcmp(node->mpName, "sp")) {
        //!!! List of ignored phonemes should be provided by some switch !!!
        node->phoneAccuracy = 0.0;
      } else {
        node->phoneAccuracy = -1.0;
        overlaped = (CorrPhnRec*)bsearch(node, corr_phn, ncorr_phns,
                            sizeof(struct CorrPhnRec), 
          cmp_maxstop);
  
        if (overlaped) {
          for (; overlaped < corr_phn + ncorr_phns &&
                overlaped->mpNode->mStart < node->mStop; overlaped++) {
            if (overlaped->mpNode->mStop  <= overlaped->mpNode->mStart ||
              overlaped->mpNode->mStop  <= node->mStart) continue;
  
            node->phoneAccuracy =
              HIGHER_OF(node->phoneAccuracy,
                        (SamePhoneme(overlaped->mpNode->mpName, node->mpName) + 1.0) *
                        (LOWER_OF(overlaped->mpNode->mStop, node->mStop) -
                          HIGHER_OF(overlaped->mpNode->mStart, node->mStart)) /
                        (overlaped->mpNode->mStop - overlaped->mpNode->mStart) - 1.0);
  
            if (node->phoneAccuracy >= 1.0) break;
          }
        }
      }
    }
    free(corr_phn);
  }
  
  #define INIT_NODE_HASH_SIZE 1000
  Node *find_or_create_node(struct my_hsearch_data *node_hash, char *node_id, Node **last)
  // Auxiliary function used by ReadHTKLattice_new. (Optionally initialize
  // uninitialized node_hash) Search for node record at key node_id. If found,
  // pointer to this record is returned, otherwise new node record is allocated
  // with type set to NT_UNDEF and entered to has at key node_id and pointer
  // to this new record is returned.
  {
    Node *node;
    ENTRY e, *ep;
  
    if (node_hash->mTabSize == 0 && !my_hcreate_r(INIT_NODE_HASH_SIZE, node_hash)) {
      Error("Insufficient memory");
    }
    e.key = node_id;
    my_hsearch_r(e, FIND, &ep, node_hash);
  
    if (ep != NULL) return (Node *) ep->data;
  
  
    node = (Node *) calloc(1, sizeof(Node));
    if (node == NULL) Error("Insufficient memory");
  
    node->mNLinks     = 0;
    node->mpLinks      = NULL;
    node->mNBackLinks = 0;
    node->mpBackLinks  = NULL;
    node->mType       = NT;
    node->mStart      = UNDEF_TIME;
    node->mStop       = UNDEF_TIME;
    node->mpPronun     = NULL;
    node->mpBackNext   = *last;
    node->phoneAccuracy = 1.0;
    *last  = node;
    e.key  = strdup(node_id);
    e.data = node;
  
    if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, node_hash)) {
      Error("Insufficient memory");
    }
    return node;
  }
  
  int RemoveCommentLines(FILE *fp)
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
  
  Node *ReadSTKNetworkInOldFormat(
    FILE *                    lfp,
    struct my_hsearch_data *  word_hash,
    struct my_hsearch_data *  phone_hash,
    LabelFormat               labelFormat,
    long                      sampPeriod,
    const char *              file_name,
    const char *              in_MLF)
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
    Node *    node;
    Node **   nodes;
  
    RemoveCommentLines(lfp);
  
    if (fscanf(lfp," %1023[^0-9]", wordOrModelName) == 1) {
      for (i=0; i<strlen(wordOrModelName); i++) {
        wordOrModelName[i] = toupper(wordOrModelName[i]);
      }
      while (--i>=0 && (wordOrModelName[i] == '='||isspace(wordOrModelName[i]))) {
        wordOrModelName[i] = '\0';
      }
    }
    if ((strcmp(wordOrModelName, "NUMNODES:") &&
        strcmp(wordOrModelName, "NUMBEROFNODES")) ||        //Obsolete NumerOfNodes
      fscanf(lfp," %d NumberOfArcs=%d", &numOfNodes, &i)<1){//Obsolete NumerOfArcs
      Error("Syntax error in file %s\nKeyword NumNodes: is missing", file_name);
    }
    if ((nodes = (Node **) calloc(numOfNodes, sizeof(Node *))) == NULL) {
      Error("Insufficient memory");
    }
    for (i=0; i < numOfNodes; i++) {
      if ((nodes[i] = (Node *) calloc(1, sizeof(Node))) == NULL) {
        Error("Insufficient memory");
      }
      nodes[i]->mType = NT_Undef;
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
              Error("Node %d is not defined in file %s", j, file_name);
            }
          }
          
        default:
          Error("Invalid syntax in definition of node %d in file %s",
                nodeId, file_name);
      }
      
      if (nodeId >= numOfNodes) {
        Error("Invalid definition of node %d in file %s.\n"
              "Node Id is bigger than number of nodes", nodeId, file_name);
      }
      node = nodes[nodeId];
  
      if (node->mType != NT_Undef)
        Error("Redefinition of node %d in file %s", nodeId, file_name);
      
      if (toupper(nodeType) != nodeType) {
        nodeType = toupper(nodeType);
        node->mType = 0;
      } else {
        node->mType = NT_True;
      }
      if (nodeType != 'M' && nodeType != 'W' && nodeType != 'N' &&
        nodeType != 'S' && nodeType != 'K' && nodeType != 'F') {
        Error("Invalid definition of node %d in file %s.\n"
              "Supported values for node type are: M - model, W - word, N - null, S - subnet, K - keyword, F - filler",
              nodeId, file_name);
      }
      node->mStart = node->mStop = UNDEF_TIME;
      node->phoneAccuracy = 1.0;
  //    node->mAux = *totalNumOfNodes;
  //    ++*totalNumOfNodes;
  
      if (nodeType == 'S') {
        FILE *snfp;
        Node *subnetFirst;
  
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
        node->mType |= NT_Phone;
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
        node->mType |= NT;
  
        if (nodeType == 'K' || nodeType == 'F') {
        node->mType |= NT_Sticky;
        nodeType = nodeType == 'K' ? 'W' : 'N';
        }
        if (nodeType == 'W') {
          ENTRY e, *ep;
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
          node->mStart = (start - center_shift - labelFormat.left_extent)  / sampPeriod;
          node->mStop  = (stop  + center_shift + labelFormat.right_extent) / sampPeriod;
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
        if ((node->mpLinks = (Link *) malloc(numOfLinks * sizeof(Link))) == NULL) {
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
  
      for (j=0; j < numOfLinks; j++) {
        if (fscanf(lfp, "%d ", &linkId) != 1) {
          Error("Invalid definition of node %d in file %s.\n"
                "Link Id is expected in link list", nodeId, file_name);
        }
        if (linkId >= numOfNodes) {
          Error("Invalid definition of node %d in file %s.\n"
                "Link Id is bigger than number of nodes", nodeId, file_name);
        }
        if (fscanf(lfp, "{%lf} ", &linkLike) != 1) {
          linkLike = 0.0;
        }
        node->mpLinks[j].mpNode = nodes[linkId];
        node->mpLinks[j].mLike = linkLike;
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
    if (!(nodes[0]           ->mType & NT) || nodes[0]           ->mpPronun != NULL ||
      !(nodes[numOfNodes-1]->mType & NT) || nodes[numOfNodes-1]->mpPronun != NULL) {
      Error("Start node and final node must be Null nodes (%s)", file_name);
    }
    for (i = 0; i < numOfNodes-1; i++) {
      nodes[i]->mpNext = nodes[i+1];
    }
  
    // create back links
    for (i = 0; i < numOfNodes; i++) {
      if (!nodes[i]->mpBackLinks) // Could be already alocated for subnetwork
        nodes[i]->mpBackLinks = (Link *) malloc(nodes[i]->mNBackLinks * sizeof(Link));
      if (nodes[i]->mpBackLinks == NULL) Error("Insufficient memory");
      nodes[i]->mNBackLinks = 0;
    }
    for (i = 0; i < numOfNodes; i++) {
      for (j=0; j < nodes[i]->mNLinks; j++) {
        Node *forwNode = nodes[i]->mpLinks[j].mpNode;
        forwNode->mpBackLinks[forwNode->mNBackLinks].mpNode = nodes[i];
        forwNode->mpBackLinks[forwNode->mNBackLinks].mLike = nodes[i]->mpLinks[j].mLike;
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
  
  Node *ReadSTKNetwork(
    FILE *lfp,
    struct my_hsearch_data *word_hash,
    struct my_hsearch_data *phone_hash,
    int notInDict,
    LabelFormat labelFormat,
    long sampPeriod,
    const char *file_name,
    const char *in_MLF)
  {
    Node *node, *enode = NULL,
        *first = NULL, *last = NULL,
        *fnode = NULL, *lnode;
    char *line;
    int line_no   =  0;
    char *chptr, *valptr, *phn_marks = NULL;
    Word *word     = NULL;
    int  i, pron_var  = 1;
    enum {LINE_START, AFTER_J, HEADER_DEF, ARC_DEF, NODE_DEF} state;
    struct my_hsearch_data node_hash = {0};
    struct readline_data   rld       = {0};
  
    for (;;) {
      do {
        if ((chptr = line = readline(lfp, &rld)) == NULL) break;
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
      while (*chptr) {
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
            if (getHTKstr(valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            node = find_or_create_node(&node_hash, valptr, &last);
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
            int nnodes  = getInteger(valptr, &chptr, file_name, line_no);
  
            if (node_hash.mTabSize == 0 && !my_hcreate_r(nnodes, &node_hash)) {
              Error("Insufficient memory");
            }
          } else { // Skip unknown header term
            if (getHTKstr(valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
          }
          continue;
        }
        if (state == NODE_DEF) {
          if ((!strcmp(chptr, "time") || !strcmp(chptr, "t")) && !(labelFormat.TIMES_OFF)) {
            char *colonptr=valptr;
            while (*colonptr && !isspace(*colonptr) && *colonptr != ',') colonptr++;
  
            if (*colonptr == ',') {
              if (colonptr != valptr) {
                *colonptr = ' ';
                node->mStart = 100 * (long long) (0.5 + 1e5 *
                              getFloat(valptr, &chptr, file_name, line_no));
              }
              valptr = colonptr+1;
            }
            node->mStop = 100 * (long long) (0.5 + 1e5 *
                        getFloat(valptr, &chptr, file_name, line_no));
          } else if (!strcmp(chptr, "var") || !strcmp(chptr, "v")) {
            pron_var = getInteger(valptr, &chptr, file_name, line_no);
            if (pron_var < 1) {
              Error("Invalid pronunciation variant (%s:%d)", file_name, line_no);
            }
          } else if (!strcmp(chptr, "p")) {
            node->phoneAccuracy = getFloat(valptr, &chptr, file_name, line_no);
          } else if (!strcmp(chptr, "flag") || !strcmp(chptr, "f")) {
            if (getHTKstr(valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            for (; *valptr; valptr++) {
              switch (toupper(*valptr)) {
                case 'K':
                case 'F':  node->mType |= NT_Sticky; break;
                case 'T':  node->mType |= NT_True;   break;
                default:
                  Error("Invalid flag '%c' (%s:%d)", *valptr, file_name, line_no);
              }
            }
          } else if (!strcmp(chptr, "L")) {
            Error("Sub-lattice nodes are not yet supported (%s:%d)",
                  *valptr, file_name, line_no);
          } else if (!strcmp(chptr, "WORD") || !strcmp(chptr, "W")) {
            ENTRY e, *ep;
            if (getHTKstr(e.key = valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            if (!strcmp(e.key, "!NULL")) {
              word = NULL;
            } else {
              my_hsearch_r(e, FIND, &ep, word_hash);
  
              if (ep != NULL) {
                word = (Word *) ep->data;
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
            node->mType &= ~(NT_Model | NT_Phone);
            node->mType |= NT;
          } else if (!strcmp(chptr, "MODEL") || !strcmp(chptr, "M")) {
            ENTRY e, *ep;
  
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
            node->mType &= ~NT;
            node->mType |= NT_Phone;
          } else if (*chptr=='\0' || !strcmp(chptr,"END") || !strcmp(chptr,"E")) {
            state = ARC_DEF;
          } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
            Error("%s (%s:%d)", chptr, file_name, line_no);
          }
          if (state == ARC_DEF || *chptr == '\0') {
            // Node definition is over. For NT, select right pronun according to
            // word and pron_var; and continue with parsing the arc definition below
            if (node->mType & NT && word != NULL) {
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
  
                  word->pronuns[i]->word       = word;
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
          }
        }
        if (state == ARC_DEF) {
          if (!*chptr || !strcmp(chptr, "END") || !strcmp(chptr, "E")) {
            if (getHTKstr(valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            enode = find_or_create_node(&node_hash, valptr, &last);
  
            int nl = ++node->mNLinks;
            node->mpLinks = (Link *) realloc(node->mpLinks, nl * sizeof(Link));
            if (node->mpLinks == NULL) Error("Insufficient memory");
            node->mpLinks[nl-1].mpNode = enode;
            node->mpLinks[nl-1].mLike = 0.0;
  
            nl = ++enode->mNBackLinks;
            enode->mpBackLinks = (Link *) realloc(enode->mpBackLinks, nl * sizeof(Link));
            if (enode->mpBackLinks == NULL) Error("Insufficient memory");
            enode->mpBackLinks[nl-1].mpNode = node;
            enode->mpBackLinks[nl-1].mLike = 0.0;
  
          } else if (!strcmp(chptr, "language") || !strcmp(chptr, "l")) {
            FLOAT mLike = getFloat(valptr, &chptr, file_name, line_no);
            //Set LM score to link pointing to enode. This link can possibly start
            //from a phone node already inserted (div=) between 'node' ans 'enode'
            Node *last = enode->mpBackLinks[enode->mNBackLinks-1].mpNode;
            last->mpLinks[last->mNLinks-1].mLike = mLike;
            enode->mpBackLinks[enode->mNBackLinks-1].mLike = mLike;
          } else if (!strcmp(chptr, "div") || !strcmp(chptr, "d")) {
            ENTRY e, *ep;
            char  name[1024];
            float time;
            int   n;
            Node  *last = node;
            FLOAT mLike  = node->mpLinks[node->mNLinks-1].mLike;
  
            if (node->mpLinks[node->mNLinks-1].mpNode != enode) {
              Error("Redefinition of  (%s:%d)", chptr, file_name, line_no);
            }
            if (getHTKstr(phn_marks=valptr, &chptr)) {
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
            time = -FLT_MAX;
            while (sscanf(phn_marks, ":%[^,:]%n,%f%n", name, &n, &time, &n) > 0) {
              Node *tnode;
              phn_marks+=n;
  
              if ((tnode            = (Node *) calloc(1, sizeof(Node))) == NULL ||
                (tnode->mpLinks     = (Link *) malloc(sizeof(Link))) == NULL ||
                (tnode->mpBackLinks = (Link *) malloc(sizeof(Link))) == NULL) {
                Error("Insufficient memory");
              }
  
              //Use special type to mark nodes inserted by d=..., they will need
              //special treatment. Later, they will become ordinary NT_Phone nodes
              tnode->mType      = NT_Phone | NT_Model;
              tnode->mNLinks    = tnode->mNBackLinks = 1;
              tnode->mpBackNext  = enode->mpBackNext;
              enode->mpBackNext  = tnode;
              tnode->phoneAccuracy = 1.0;
              //Store phone durations now. Will be replaced by absolute times below.
              tnode->mStart  = time != -FLT_MAX ? 100 * (long long) (0.5 + 1e5 * time) : UNDEF_TIME;
              tnode->mStop   = UNDEF_TIME;
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
              last->mpLinks[last->mNLinks-1].mpNode = tnode;
              last->mpLinks[last->mNLinks-1].mLike = 0.0;
              tnode->mpBackLinks[0].mpNode = last;
              tnode->mpBackLinks[0].mLike = 0.0;
              last = tnode;
            }
            if (strcmp(phn_marks,":")) {
              Error("Invalid specification of phone marks (d=) (%s:%d)",
                    file_name, line_no);
            }
            last->mpLinks[last->mNLinks-1].mpNode = enode;
            last->mpLinks[last->mNLinks-1].mLike = mLike;
            enode->mpBackLinks[enode->mNBackLinks-1].mpNode = last;
            enode->mpBackLinks[enode->mNBackLinks-1].mLike = mLike;
  
          } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
            Error("%s (%s:%d)", chptr, file_name, line_no);
          }
        }
      }
    }
    my_hdestroy_r(&node_hash, 1);
    lnode = last;
    first = last = NULL;
    if (lnode) lnode->mpNext = NULL;
    for (node = lnode; node != NULL; fnode = node, node = node->mpBackNext) {
      if (node->mpBackNext) node->mpBackNext->mpNext = node;
  
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
      //If only stop time is specified, set start time to lowest predecessor stop time
      if (node->mStart == UNDEF_TIME && node->mStop != UNDEF_TIME) {
        int i;
        for (i = 0; i < node->mNBackLinks; i++) {
          Node *backnode = node->mpBackLinks[i].mpNode;
          //When seraring predecessors, skip nodes inserted by d=...
          while (backnode->mType == (NT_Phone | NT_Model)) {
            assert(backnode->mNBackLinks == 1);
            backnode = backnode->mpBackLinks[0].mpNode;
          }
          if (backnode->mStop != UNDEF_TIME) {
            node->mStart = node->mStart == UNDEF_TIME
                          ?          backnode->mStop
                          : LOWER_OF(backnode->mStop, node->mStart);
          }
        }
        if (node->mStart == UNDEF_TIME) node->mStart = 0;
      }
      //For model nodes defined by d=... (NT_Phone | NT_Model), node->mStart contains
      //only phone durations. Absolute times must be computed derived starting from
      //the end time of the node to which arc with d=... definition points.
      if (node->mType == (NT_Phone | NT_Model)) {
        assert(node->mNLinks == 1);
        node->mStop = node->mpLinks[0].mpNode->mType == (NT_Phone | NT_Model)
                    && node->mStart != UNDEF_TIME
                    ? node->mpLinks[0].mpNode->mStart : node->mpLinks[0].mpNode->mStop;
        node->mStart = node->mStart != UNDEF_TIME && node->mStop != UNDEF_TIME
                      ? node->mStop - node->mStart : node->mpLinks[0].mpNode->mStart;
      }
    }
    if (!first || !last) {
      Error("Network contain no start node or no final node (%s)", file_name);
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
      if (node->mType == (NT_Phone | NT_Model)) {
        node->mType = NT_Phone;
      }
      if (node->mStart != UNDEF_TIME) {
        node->mStart = (node->mStart - labelFormat.left_extent) / sampPeriod;
      }
      if (node->mStop  != UNDEF_TIME) {
        node->mStop  = (node->mStop + labelFormat.right_extent) / sampPeriod;
      }
    }
    if (first->mpPronun != NULL) {
      node = (Node *) calloc(1, sizeof(Node));
      if (node == NULL) Error("Insufficient memory");
      node->mpNext       = first;
      node->mpBackNext   = NULL;
      first->mpBackNext  = node;
      node->mType       = NT;
      node->mpPronun     = NULL;
      node->mStart      = UNDEF_TIME;
      node->mStop       = UNDEF_TIME;
      node->mNBackLinks = 0;
      node->mpBackLinks  = NULL;
      node->mNLinks     = 1;
      node->mpLinks      = (Link*) malloc(sizeof(Link));
      if (node->mpLinks == NULL) Error("Insufficient memory");
      node->mpLinks[0].mLike = 0.0;
      node->mpLinks[0].mpNode = first;
      first->mNBackLinks = 1;
      first->mpBackLinks  = (Link*) malloc(sizeof(Link));
      if (first->mpBackLinks == NULL) Error("Insufficient memory");
      first->mpBackLinks[0].mLike = 0.0;
      first->mpBackLinks[0].mpNode = node;
      first = node;
    }
    if (last->mpPronun != NULL) {
      node = (Node *) calloc(1, sizeof(Node));
      if (node == NULL) Error("Insufficient memory");
      last->mpNext      = node;
      node->mpNext      = NULL;
      node->mpBackNext  = last;
      node->mType      = NT;
      node->mpPronun    = NULL;
      node->mStart     = UNDEF_TIME;
      node->mStop      = UNDEF_TIME;
      node->mNLinks    = 0;
      node->mpLinks     = NULL;
      last->mNLinks = 1;
      last->mpLinks  = (Link*) malloc(sizeof(Link));
      if (last->mpLinks == NULL) Error("Insufficient memory");
      last->mpLinks[0].mLike = 0.0;
      last->mpLinks[0].mpNode = node;
      node->mNBackLinks = 1;
      node->mpBackLinks  = (Link*) malloc(sizeof(Link));
      if (node->mpBackLinks == NULL) Error("Insufficient memory");
      node->mpBackLinks[0].mLike = 0.0;
      node->mpBackLinks[0].mpNode = last;
    }
    return first;
  }
  
  
  // Debug function showing network using AT&T dot utility
  void dnet(Node *net, int nAuxNodePtrs, ...)
  {
    static int dnetcnt=1;
    va_list ap;
    Node *node;
    int i = 1;
  
    FILE *fp = popen("cat | (tf=`mktemp /tmp/netps.XXXXXX`;"
                    "dot -Tps > $tf; gv -scale -4 $tf; rm $tf)",
                    "w");
  //  FILE *fp = stdout;
  
    if (fp == NULL) return;
  
    for (node = net; node != NULL; node = node->mpNext) {
      node->estate_id = i++;
    }
    fprintf(fp, "digraph \"dnet%d\" {\nrankdir=LR\n", dnetcnt++);
  
  
    for (node = net; node != NULL; node = node->mpNext) {
      fprintf(fp, "n%d [shape=%s,label=\"%d:%s", node->estate_id,
              node->mType & NT ? "box" : "ellipse", node->estate_id,
              node->mType & NT ? (node->mpPronun ?
                                      node->mpPronun->word->mpName : "-"):
              node->mType & NT_Phone? node->mpName :
              node->mType & NT_Model? node->mpHmm->mpMacro->mpName : "???");
  
      if (node->mType & NT && node->mpPronun != NULL) {
        if (node->mpPronun != node->mpPronun->word->pronuns[0]) {
          fprintf(fp, ":%d", node->mpPronun->variant_no);
        }
        fprintf(fp, "\\n");
  
        if (node->mpPronun->outSymbol != node->mpPronun->word->mpName) {
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
  //             node->estate_id, node->mpNext->estate_id);
  //    }
  //    if (node->mpBackNext != NULL) {
  //     fprintf(fp,"n%d -> n%d [color=gray,weight=1]\n",
  //             node->estate_id, node->mpBackNext->estate_id);
  //    }
    }
    for (node = net; node != NULL; node = node->mpNext) {
      for (i = 0; i < node->mNLinks; i++) {
        fprintf(fp,"n%d -> n%d [color=blue,weight=1",
                node->estate_id,node->mpLinks[i].mpNode->estate_id);
        if (node->mpLinks[i].mLike != 0.0) {
          fprintf(fp,",label=\""FLOAT_FMT"\"", node->mpLinks[i].mLike);
        }
        fprintf(fp,"];\n");
      }
  //    for (i = 0; i < node->mNBackLinks; i++) {
  //      fprintf(fp,"n%d -> n%d [color=red,weight=1",
  //              node->estate_id,node->mpBackLinks[i].mpNode->estate_id);
  //      if (node->mpBackLinks[i].mLike != 0.0) {
  //        fprintf(fp,",label=\""FLOAT_FMT"\"", node->mpBackLinks[i].mLike);
  //      }
  //      fprintf(fp,"];\n");
  //    }
    }
    va_start(ap, nAuxNodePtrs);
    for (i = 0; i < nAuxNodePtrs; i++) {
      Node *ptr = va_arg(ap, Node *);
      fprintf(fp, "AuxPtr%d [shape=plaintext];\nAuxPtr%d -> n%d\n",
              i, i, ptr->estate_id);
    }
    va_end(ap);
  
    fprintf(fp, "}\n");
    pclose(fp);
  }
  
  
  Node *ReadHTKLattice(
    FILE *lfp,
    struct my_hsearch_data *word_hash,
    struct my_hsearch_data *phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char *file_name)
  {
    Node **nodes = NULL, *node,
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
          Error("'=' expected (%s:%d)");
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
            nodes  = (Node **) calloc(nnodes, sizeof(Node *));
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
            nodes[node_id] = (Node *) calloc(1, sizeof(Node));
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
            ENTRY e, *ep;
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
        nodes[node_id]->mType       = NT;
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
        nodes[node_id]->mStart     = UNDEF_TIME;
        nodes[node_id]->mStop      = node_time + labelFormat.right_extent;
      } else if (state == ARC_DEF) {
        if (arc_start == -1 || arc_end == -1) {
          Error("Start node or end node not defined (%s:%d)", file_name, line_no);
        }
  
        int linkId = nodes[arc_start]->mNLinks++;
        nodes[arc_start]->mpLinks =
          (Link *) realloc(nodes[arc_start]->mpLinks, (linkId+1) * sizeof(Link));
  
        if (nodes[arc_start]->mpLinks == NULL) Error("Insufficient memory");
  
        last = nodes[arc_start];
  
        if (phn_marks) {
          ENTRY e, *ep;
          char  name[1024];
          float time;
          int   n;
  
          while (sscanf(phn_marks, ":%[^,],%f%n", name, &time, &n) > 1) {
            phn_marks+=n;
  
            if ((node            = (Node *) calloc(1, sizeof(Node))) == NULL ||
              (node->mpLinks     = (Link *) malloc(sizeof(Link))) == NULL ||
              (node->mpBackLinks = (Link *) malloc(sizeof(Link))) == NULL) {
              Error("Insufficient memory");
            }
            node->mType = NT_Phone;
            node->mNLinks = node->mNBackLinks = 1;
            node->mpNext   = last->mpNext;
            last->mpNext = node;
            node->phoneAccuracy = 1.0;
  
            if (!(labelFormat.TIMES_OFF)) {
              node->mStart = last->mStop - labelFormat.left_extent - labelFormat.right_extent;
              node->mStop  = last->mStop + 100 * (long long) (0.5 + 1e5 * time)
                                      + labelFormat.right_extent;
            } else {
              node->mStart  = UNDEF_TIME;
              node->mStop   = UNDEF_TIME;
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
            last->mpLinks[last->mNLinks-1].mpNode = node;
            last->mpLinks[last->mNLinks-1].mLike = 0.0;
            node->mpBackLinks[0].mpNode = last;
            node->mpBackLinks[0].mLike = 0.0;
            last = node;
          }
          if (strcmp(phn_marks,":")) {
            Error("Invalid specification of phone marks (d=) (%s:%d)",
                  file_name, line_no);
          }
        }
        linkId = nodes[arc_end]->mNBackLinks++;
        nodes[arc_end]->mpBackLinks =
          (Link *) realloc(nodes[arc_end]->mpBackLinks, (linkId+1) * sizeof(Link));
  
        if (nodes[arc_end]->mpBackLinks == NULL) Error("Insufficient memory");
  
        last->mpLinks[last->mNLinks-1].mpNode = nodes[arc_end];
        last->mpLinks[last->mNLinks-1].mLike = arc_like;
        nodes[arc_end]->mpBackLinks[linkId].mpNode = last;
        nodes[arc_end]->mpBackLinks[linkId].mLike = arc_like;
  
        if (nodes[arc_start]->mStop != UNDEF_TIME) {
          if (nodes[arc_end]->mStart == UNDEF_TIME) {
            nodes[arc_end]->mStart = nodes[arc_start]->mStop - labelFormat.left_extent
                                                          - labelFormat.right_extent;
          } else {
            nodes[arc_end]->mStart = LOWER_OF(nodes[arc_end]->mStart,
                                            nodes[arc_start]->mStop - labelFormat.left_extent
                                                                    - labelFormat.right_extent);
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
      node = (Node *) calloc(1, sizeof(Node));
      if (node == NULL) Error("Insufficient memory");
      node->mpNext       = first;
      node->mpBackNext   = NULL;
      first->mpBackNext  = node;
      node->mType       = NT;
      node->mpPronun     = NULL;
      node->mStart      = UNDEF_TIME;
      node->mStop       = UNDEF_TIME;
      node->mNBackLinks = 0;
      node->mpBackLinks  = NULL;
      node->mNLinks     = 1;
      node->mpLinks      = (Link*) malloc(sizeof(Link));
      if (node->mpLinks == NULL) Error("Insufficient memory");
      node->mpLinks[0].mLike = 0.0;
      node->mpLinks[0].mpNode = first;
  
      first->mNBackLinks = 1;
      first->mpBackLinks  = (Link*) malloc(sizeof(Link));
      if (first->mpBackLinks == NULL) Error("Insufficient memory");
      first->mpBackLinks[0].mLike = 0.0;
      first->mpBackLinks[0].mpNode = node;
      first = node;
    }
    if (last->mpPronun != NULL) {
      node = (Node *) calloc(1, sizeof(Node));
      if (node == NULL) Error("Insufficient memory");
      last->mpNext      = node;
      node->mpNext      = NULL;
      node->mpBackNext  = last;
  
      node->mType      = NT;
      node->mpPronun    = NULL;
      node->mStart     = UNDEF_TIME;
      node->mStop      = UNDEF_TIME;
      node->mNLinks    = 0;
      node->mpLinks     = NULL;
  
      last->mNLinks = 1;
      last->mpLinks  = (Link*) malloc(sizeof(Link));
      if (last->mpLinks == NULL) Error("Insufficient memory");
      last->mpLinks[0].mLike = 0.0;
      last->mpLinks[0].mpNode = node;
  
      node->mNBackLinks = 1;
      node->mpBackLinks  = (Link*) malloc(sizeof(Link));
      if (node->mpBackLinks == NULL) Error("Insufficient memory");
      node->mpBackLinks[0].mLike = 0.0;
      node->mpBackLinks[0].mpNode = last;
    }
    return first;
  }
  
  void WriteSTKNetworkInOldFormat(
    FILE        *lfp,
    Node       *first,
    LabelFormat labelFormat,
    long        sampPeriod,
    const char  *net_file,
    const char  *out_MNF)
  {
    int i;
    Node *node;
  
    for (i = 0, node = first; node != NULL; node = node->mpNext, i++)  {
      node->mAux = i;
    }
    fprintf(lfp,"NUMNODES: %d\n", i);
    for (i = 0, node = first; node != NULL; node = node->mpNext, i++) {
      int j,
      type = node->mType & NT_Model       ? 'M'  :
            node->mType & NT_Phone       ? 'M'  :
            node->mType & NT_Subnet      ? 'S'  :
            node->mType & NT        ?
              (node->mpPronun == NULL     ?
                (node->mType & NT_Sticky ? 'F'  :
                                          'N') :
                (node->mType & NT_Sticky ? 'K'  :
                                          'W')):
                                          '?';
      if (!(node->mType & NT_True)) {
        type = tolower(type);
      }
      fprintf(lfp,"%d\t%c %s",
              i, type,
              node->mType & NT_Model   ? node->mpHmm->mpMacro->mpName :
              node->mType & NT_Phone   ? node->mpName :
              node->mType & NT_Subnet  ? node->mpName :
              node->mType & NT    ?
                (node->mpPronun == NULL ? "-" :
                                        node->mpPronun->word->mpName):
                                        "?");
      if (node->mType & NT && node->mpPronun) {
        if (node->mpPronun->word->mpName != node->mpPronun->outSymbol) {
          fprintf(lfp," [%s]", node->mpPronun->outSymbol);
        }
        if (node->mpPronun->prob != 0.0 || node->mpPronun->word->npronuns > 1) {
          fprintf(lfp," {%d "FLOAT_FMT"}",
                  node->mpPronun->variant_no,
                  node->mpPronun->prob);
        }
      }
      if (node->mType & NT_Phone && node->phoneAccuracy != 1.0) {
        fprintf(lfp," {"FLOAT_FMT"}", node->phoneAccuracy);
      }
      if (!(labelFormat.TIMES_OFF) &&
        node->mStart != UNDEF_TIME && node->mStop != UNDEF_TIME) {
        int ctm = labelFormat.CENTRE_TM;
        fprintf(lfp," (%lld %lld)",
                    (long long) sampPeriod * (2 * node->mStart + ctm) / 2 - labelFormat.left_extent,
                    (long long) sampPeriod * (2 * node->mStop - ctm)  / 2 + labelFormat.right_extent);
      }
      fprintf(lfp,"\t%d", node->mNLinks);
      for (j = 0; j < node->mNLinks; j ++) {
        fprintf(lfp," %d",node->mpLinks[j].mpNode->mAux);
        if (node->mpLinks[j].mLike != 0.0) {
          fprintf(lfp," {"FLOAT_FMT"}", node->mpLinks[j].mLike);
        }
      }
      fputs("\n", lfp);
      if (ferror(lfp)) {
        Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
      }
    }
  }
  
  void NetworkExpansionsAndOptimizations(
    Node *node,
    ExpansionOptions expOptions,
    STKNetworkOutputFormat out_net_fmt,
    struct my_hsearch_data *wordHash,
    struct my_hsearch_data *nonCDphHash,
    struct my_hsearch_data *triphHash)
  {
    if (expOptions.no_word_expansion && !expOptions.CD_phone_expansion
    && expOptions.no_optimization    && !out_net_fmt.no_LM_likes
    &&!out_net_fmt.no_times          && !out_net_fmt.no_word_nodes &&
      !out_net_fmt.no_model_nodes    && !out_net_fmt.no_pronun_vars) return;
    SelfLinksToNullNodes(node);
    if (!expOptions.no_word_expansion) {
      if (!expOptions.no_optimization) {
        LatticeLocalOptimization(node, expOptions.strict_timing, expOptions.trace_flag);
      }
      ExpandWordNetworkByDictionary(node, wordHash, !expOptions.remove_words_nodes,
                                                    !expOptions.respect_pronun_var);
    }
    if (expOptions.CD_phone_expansion) {
      if (!expOptions.no_optimization) {
        LatticeLocalOptimization(node, expOptions.strict_timing, expOptions.trace_flag);
      }
      ExpandMonophoneNetworkToTriphones(node, nonCDphHash, triphHash);
    }
    DiscardUnwantedInfoInNetwork(node, out_net_fmt);
  
    if (!expOptions.no_optimization) {
      LatticeLocalOptimization(node, expOptions.strict_timing, expOptions.trace_flag);
    }
    RemoveRedundantNullNodes(node);
  }

} // namespace STK
