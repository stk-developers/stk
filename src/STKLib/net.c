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
#include "net.h"
//#include "viterbi.h"
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

void FreeNetwork(Node *node) {
  Node *tnode;
  while (node) {
    tnode = node->next;
    free(node->links);
    free(node->backlinks);
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
     (first->links    = (Link *) malloc(sizeof(Link))) == NULL ||
     (last->backlinks = (Link *) malloc(sizeof(Link))) == NULL) {
    Error("Insufficient memory");
  }
  first->mType          = last->mType          = NT;
  first->pronun        = last->pronun        = NULL;
  first->nlinks        = last->nbacklinks    = 1;
  first->nbacklinks    = last->nlinks        = 0;
  first->backlinks     = last->links         = NULL;
  first->start         = last->start         = UNDEF_TIME;
  first->stop          = last->stop          = UNDEF_TIME;
//  first->tokens        = last->tokens        = NULL;
//  first->exitToken     = last->exitToken     = NULL;

  node = first;
  for (lp = labels; lp != NULL; lp = lp->next) {
    Node *tnode;

    if ((tnode            = (Node *) calloc(1, sizeof(Node))) == NULL ||
       (tnode->links     = (Link *) malloc(sizeof(Link))) == NULL ||
       (tnode->backlinks = (Link *) malloc(sizeof(Link))) == NULL) {
      Error("Insufficient memory");
    }
    node->links[0].node    = tnode;
    node->links[0].like    = 0.0;
    switch (node_type) {
      case NT:  tnode->pronun = ((Word *) lp->data)->pronuns[0]; break;
      case NT_Model: tnode->hmm    =   (Hmm *) lp->data;              break;
      case NT_Phone: tnode->mpName   =  (char *) lp->data;              break;
      default:       Error("Fatal: Invalid node type");
    }
    tnode->mType       = node_type;
    tnode->nlinks     = 1;
    tnode->nbacklinks = 1;
    tnode->start      = lp->start;
    tnode->stop       = lp->stop;
    tnode->backlinks[0].node = node;
    tnode->backlinks[0].like = 0.0;
    node->next = tnode;
    node = tnode;
  }
  node->next = last;
  node->links[0].node    = last;
  node->links[0].like    = 0.0;
  last->backlinks[0].node = node;
  last->backlinks[0].like = 0.0;
  last->next = NULL;
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

  assert(first != NULL || first->mType & NT || first->pronun == NULL);

  for (node = first; node != NULL; prev = node, node = node->next) {

    if (!(node->mType & NT)) continue;

    if (node->pronun == NULL) continue;
    Word *word = node->pronun->word;

    //Do not expand non-dictionary words, which where added by ReadSTKNetwork
    if (word->npronunsInDict == 0) continue;

    if (!multiple_pronun) {
      singlePronunWrd.mpName = node->pronun->word->mpName;
      word = &singlePronunWrd;
      *word->pronuns = node->pronun;
    }

    // Remove links to current node form backlinked nodes and realloc
    // link arrays of backlinked nodes to hold word->npronuns more backlinks
    for (i = 0; i < node->nbacklinks; i++) {
      Node *bakcnode = node->backlinks[i].node;
      for (j=0; j<bakcnode->nlinks && bakcnode->links[j].node!=node; j++);
      assert(j < bakcnode->nlinks); // Otherwise link to 'node' is missing
                                    // from which backlink exists
      bakcnode->links[j] = bakcnode->links[bakcnode->nlinks-1];

      bakcnode->links = (Link *)
        realloc(bakcnode->links,
               (bakcnode->nlinks - 1 + word->npronuns) * sizeof(Link));
      if (bakcnode->links == NULL) Error("Insufficient memory");
      bakcnode->nlinks--;// += word->npronuns-1;
    }

    // Remove backlinks to current node form linked nodes and realloc
    // backlink arrays of linked nodes to hold word->npronuns more backlinks
    for (i=0; i < node->nlinks; i++) {
      Node *forwnode = node->links[i].node;
      for (j=0;j<forwnode->nbacklinks&&forwnode->backlinks[j].node!=node;j++);
      assert(j < forwnode->nbacklinks);
      // Otherwise link to 'node' is missing from which backlink exists
      forwnode->backlinks[j] = forwnode->backlinks[forwnode->nbacklinks-1];

      forwnode->backlinks = (Link *)
        realloc(forwnode->backlinks,
              (forwnode->nbacklinks - 1 + word->npronuns) * sizeof(Link));
      if (forwnode->backlinks == NULL) Error("Insufficient memory");
      forwnode->nbacklinks--;
    }
    for (i = 0; i < word->npronuns; i++) {
      Pronun *pronun = word->pronuns[i];
      Node *pronun_first = NULL, *pronun_prev = NULL, *tnode;

      for (j = 0; j < pronun->nmodels; j++) {
        tnode = (Node *) calloc(1, sizeof(Node));
        if (tnode == NULL) Error("Insufficient memory");

        tnode->mType       = NT_Phone | (node->mType & NT_True);
        tnode->mpName       = pronun->model[j].mpName;
        tnode->start      = node->start;
        tnode->stop       = node->stop;
        tnode->phoneAccuracy = 1.0;

        if (j == 0) {
          pronun_first = tnode;
        } else {
          if ((pronun_prev->links = (Link *) malloc(sizeof(Link))) == NULL ||
             (tnode->backlinks   = (Link *) malloc(sizeof(Link))) == NULL) {
            Error("Insufficient memory");
          }
          tnode->nbacklinks          = 1;
          tnode->backlinks[0].node   = pronun_prev;
          tnode->backlinks[0].like   = 0.0;
          pronun_prev->nlinks        = 1;
          pronun_prev->links[0].node = tnode;
          pronun_prev->links[0].like = 0.0;
          pronun_prev->next          = tnode;
        }
        pronun_prev = tnode;
      }
      if (keep_word_nodes || j == 0) {
        tnode = (Node *) calloc(1, sizeof(Node));
        if (tnode == NULL) Error("Insufficient memory");

        tnode->mType       = NT | (node->mType & NT_True);
        tnode->pronun     = keep_word_nodes ? word->pronuns[i] : NULL;
        tnode->start      = node->start;
        tnode->stop       = node->stop;

        if (j == 0) {
          pronun_first = tnode;
        } else {
          if ((pronun_prev->links = (Link *) malloc(sizeof(Link))) == NULL ||
             (tnode->backlinks   = (Link *) malloc(sizeof(Link))) == NULL) {
            Error("Insufficient memory");
          }
          tnode->nbacklinks          = 1;
          tnode->backlinks[0].node   = pronun_prev;
          tnode->backlinks[0].like   = 0.0;
          pronun_prev->nlinks        = 1;
          pronun_prev->links[0].node = tnode;
          pronun_prev->links[0].like = 0.0;
          pronun_prev->next          = tnode;
        }
        pronun_prev = tnode;
      }
      if ((pronun_prev->links =
            (Link *) malloc(sizeof(Link) * node->nlinks))==NULL ||
         (pronun_first->backlinks =
            (Link *) malloc(sizeof(Link) * node->nbacklinks)) == NULL) {
        Error("Insufficient memory");
      }
      pronun_prev->nlinks      = node->nlinks;
      pronun_first->nbacklinks = node->nbacklinks;

      for (j = 0; j < node->nbacklinks; j++) {
        Node *backnode = node->backlinks[j].node;
        backnode->links[backnode->nlinks  ].node = pronun_first;
        backnode->links[backnode->nlinks++].like = node->backlinks[j].like;
        pronun_first->backlinks[j] = node->backlinks[j];
      }
      for (j=0; j < node->nlinks; j++) {
        Node *forwnode = node->links[j].node;
        forwnode->backlinks[forwnode->nbacklinks  ].node = pronun_prev;
        forwnode->backlinks[forwnode->nbacklinks++].like = node->links[j].like;
        pronun_prev->links[j] = node->links[j];
      }
      if (prev != NULL) prev->next = pronun_first;
      prev = pronun_prev;
    }
    prev->next = node->next;
    free(node->links);
    free(node->backlinks);
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

  for (node = first; node != NULL; node = node->next)  {
    if (format.no_LM_likes) {
      for (i=0; i < node->nlinks;     i++) node->links    [i].like = 0.0;
      for (i=0; i < node->nbacklinks; i++) node->backlinks[i].like = 0.0;
    }
    if (format.no_times) {
      node->stop = node->start = UNDEF_TIME;
    }
    if (format.no_word_nodes && node->mType & NT) {
      node->pronun = NULL;
    }
    if (format.no_model_nodes && (node->mType&NT_Model || node->mType&NT_Phone)) {
      node->mType = NT;
      node->pronun = NULL;
    }
    if (format.no_pronun_vars && node->mType & NT && node->pronun != NULL) {
      node->pronun = node->pronun->word->pronuns[0];
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

  for (n = 0, node = first; node != NULL; node = node->next, n++)  {
    node->aux = n;
    l += node->nlinks;
  }
  fprintf(lfp,"N=%d L=%d\n", n, l);
  for (node = first; node != NULL; node = node->next) {
    int j;

    if (format.all_field_names) fputs("I=", lfp);
    if (format.base62_labels) fprintBase62(lfp, node->aux);
    else                     fprintf(lfp,"%d", node->aux);

    if (!format.no_times && node->stop != UNDEF_TIME) {
      fputs(" t=", lfp);

      if (node->start != UNDEF_TIME && format.start_times) {
        fprintf(lfp,"%g,", node->start * 1.0e-7 * sampPeriod);
      }
      fprintf(  lfp,"%g",  node->stop  * 1.0e-7 * sampPeriod);
    }
    if (!(node->mType & NT && node->pronun == NULL)
       || !format.no_defaults) {
      putc(' ', lfp);
      putc(node->mType & NT   ? 'W' :
           node->mType & NT_Subnet ? 'S' :
                                    'M', lfp); // NT_Model, NT_Phone
      putc('=', lfp);
      fprintHTKstr(lfp, node->mType & NT_Model   ? node->hmm->mpMacro->mpName   :
                        node->mType & NT    ? (!node->pronun ? "!NULL" :
                                                   node->pronun->word->mpName) :
                                                  node->mpName); // NT_PHONE (NT_Subnet)
    }
    if (!format.no_pronun_vars && node->mType & NT
    && node->pronun != NULL && node->pronun->word->npronuns > 1
    && (node->pronun->variant_no > 1 || !format.no_defaults)) {
      fprintf(lfp," v=%d", node->pronun->variant_no);
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
      if (format.all_field_names) fprintf(lfp," J=%d", node->nlinks);

      for (j = 0; j < node->nlinks; j ++) {
        putc(' ', lfp);
        if (format.all_field_names) fputs("E=", lfp);
        if (format.base62_labels) fprintBase62(lfp, node->links[j].node->aux);
        else                     fprintf(lfp,"%d", node->links[j].node->aux);
        if (node->links[j].like != 0.0 && !format.no_LM_likes) {
          fprintf(lfp," l="FLOAT_FMT, node->links[j].like);
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
    for (node = first; node != NULL; node = node->next) {
      int j;

      for (j = 0; j < node->nlinks; j ++) {
        if (format.all_field_names) {
          fprintf(lfp, format.arc_defs_with_J ? "J=%d S=" : "I=", l++);
        }
        if (format.base62_labels) fprintBase62(lfp, node->aux);
        else                     fprintf(lfp,"%d", node->aux);
        putc(' ', lfp); // space = ' ';
        if (format.all_field_names) fputs("E=", lfp);
        if (format.base62_labels) fprintBase62(lfp, node->links[j].node->aux);
        else                     fprintf(lfp,"%d", node->links[j].node->aux);
        if (node->links[j].like != 0.0 && !format.no_LM_likes) {
          fprintf(lfp," l="FLOAT_FMT, node->links[j].like);
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
//  return ((Link *) a)->node - ((Link *) b)->node;
//  Did not work with gcc, probably bug in gcc pointer arithmetic
  return (char *)((Link *) a)->node - (char *)((Link *) b)->node;
}


static int LatticeLocalOptimization_ForwardPass(Node *first, int strictTiming)
{
  int i, j, k, l, m, rep;
  Node *node, *tnode;
  int node_removed = 0;
  FLOAT tlike;
  for (node = first; node != NULL; node = node->next) {
    for (i = 0; i < node->nlinks; i++) {
    //for (tnode = inode; tnode != NULL; tnode = (tnode == inode ? jnode : NULL)) {
      tnode = node->links[i].node;
      if (tnode->nlinks == 0) continue;

      // Weight pushing
      tlike = tnode->backlinks[0].like;
      for (l=1; l <  tnode->nbacklinks; l++) {
        tlike = HIGHER_OF(tlike, tnode->backlinks[l].like);
      }
      for (l=0; l < tnode->nbacklinks; l++) {
        Node *backnode = tnode->backlinks[l].node;
        tnode->backlinks[l].like -= tlike;
        for (k=0; k<backnode->nlinks && backnode->links[k].node!=tnode; k++);
        assert(k < backnode->nlinks);
        backnode->links[k].like -= tlike;
    #ifndef NDEBUG
        for (k++; k<backnode->nlinks && backnode->links[k].node!=tnode; k++);
    #endif
        assert(k == backnode->nlinks);
      }
      for (l=0; l < tnode->nlinks; l++) {
        Node *forwnode = tnode->links[l].node;
        tnode->links[l].like += tlike;
        for (k=0; k<forwnode->nbacklinks && forwnode->backlinks[k].node!=tnode;k++);
        assert(k < forwnode->nbacklinks);
        forwnode->backlinks[k].like += tlike;
    #ifndef NDEBUG
        for (k++; k<forwnode->nbacklinks && forwnode->backlinks[k].node!=tnode;k++);
    #endif
        assert(k == forwnode->nbacklinks);
      }
    }
//dnet(first, 1, node);

    for (i = 0; i < node->nlinks-1; i++) {
      for (j = i+1; j < node->nlinks; j++) {
        Node *inode = node->links[i].node;
        Node *jnode = node->links[j].node;

        if ((inode->mType & ~NT_True) != (jnode->mType & ~NT_True)
        || ( inode->mType & NT_Phone && inode->mpName   != jnode->mpName)
        || ( inode->mType & NT  && inode->pronun != jnode->pronun)
//          &&  (inode->pronun == NULL ||
//             jnode->pronun == NULL ||
//             inode->pronun->word       != jnode->pronun->word ||
//             inode->pronun->outSymbol  != jnode->pronun->outSymbol ||
//             inode->pronun->variant_no != jnode->pronun->variant_no ||
//             inode->pronun->prob       != jnode->pronun->prob)
        || (inode->nbacklinks != jnode->nbacklinks)) {
         continue;
        }
        if (strictTiming && (inode->start != jnode->start
                        ||  inode->stop  != jnode->stop)) {
          continue;
        }
        for (l=0; l < inode->nbacklinks; l++) {
          if (inode->backlinks[l].node != jnode->backlinks[l].node) break;
          FLOAT ldiff =  inode->backlinks[l].like - jnode->backlinks[l].like;
          if (ldiff < -SIGNIFICANT_PROB_DIFFERENCE ||
             ldiff >  SIGNIFICANT_PROB_DIFFERENCE) break;
        }
        if (l < inode->nbacklinks) continue;

/*        if (memcmp(inode->backlinks, jnode->backlinks,
                  inode->nbacklinks * sizeof(Link))) {
          continue;
        }*/
          // inode and jnode are the same nodes with the same predeccessors
          // Remove jnode and add its links to inode

        assert(inode->nlinks && jnode->nlinks);

          // Remove links to jnode form predeccessors
        for (l=0; l < jnode->nbacklinks; l++) {
          Node *backnode = jnode->backlinks[l].node;
          for (k=0; k<backnode->nlinks && backnode->links[k].node!=jnode; k++);
          assert(k < backnode->nlinks);
          // Otherwise link to 'node' is missing from which backlink exists
          memmove(backnode->links+k, backnode->links+k+1,
                  (backnode->nlinks-k-1) * sizeof(Link));
          backnode->nlinks--;
        }

        // Merge jnode's links with inode links

        //Count jnode's links not present among inode's links
        rep = l = k = 0;
        while (k < jnode->nlinks) {
          Link *ill = inode->links+l;
          Link *jlk = jnode->links+k;
          if (l == inode->nlinks || ill->node > jlk->node){
            // k-th link of jnode will be included among inode's links.
            // Redirect corresponding baclink to inode
            for (m = 0; m < jlk->node->nbacklinks
                        && jlk->node->backlinks[m].node != jnode; m++);
            assert(m < jlk->node->nbacklinks);
            jlk->node->backlinks[m].node = inode;
            qsort(jlk->node->backlinks, jlk->node->nbacklinks,
                  sizeof(Link), lnkcmp);
            k++;
          } else  if (ill->node == jlk->node) {
            // l-th link of inode and k-th link of jnode points to
            // the same node. Link from jnode is redundant.
            // Remove backlinks to jnode form jnode's succesors
            for (m = 0; m < jlk->node->nbacklinks
                      && jlk->node->backlinks[m].node != jnode; m++);
            assert(m < jlk->node->nbacklinks);
            memmove(jlk->node->backlinks+m, jlk->node->backlinks+m+1,
                    (jlk->node->nbacklinks-m-1) * sizeof(Link));
            jlk->node->nbacklinks--;

            ill->like = HIGHER_OF(ill->like, jlk->like);
            jlk->node = NULL; // Mark link to be removed
            rep++; k++, l++;
          } else {
            l++;
          }
        }
        l = inode->nlinks;
        inode->nlinks += jnode->nlinks-rep;
        inode->links = (Link *) realloc(inode->links,
                                        inode->nlinks * sizeof(Link));
        if (inode->links == NULL) Error("Insufficient memory");

        for (k = 0; k < jnode->nlinks; k++) {
          if (jnode->links[k].node != NULL) {
            inode->links[l++] = jnode->links[k];
          }
        }
        qsort(inode->links, inode->nlinks, sizeof(Link), lnkcmp);

        inode->start = inode->start == UNDEF_TIME || jnode->start == UNDEF_TIME
                       ? UNDEF_TIME : LOWER_OF(inode->start, jnode->start);

        inode->stop  = inode->stop == UNDEF_TIME || jnode->stop == UNDEF_TIME
                       ? UNDEF_TIME : HIGHER_OF(inode->stop, jnode->stop);

        if (inode->aux > jnode->aux) {
        // Make sure that topological order of new inode's links
        // (inherited from jnode) is higher than inode's order.
        // In the 'next' list, move inode to jnode's lower position
          inode->backnext->next = inode->next;
          inode->next->backnext = inode->backnext;
          inode->next           = jnode->next;
          inode->backnext       = jnode->backnext;
          inode->backnext->next = inode;
          inode->next->backnext = inode;
          inode->aux = jnode->aux;
        } else {
          jnode->next->backnext = jnode->backnext;
          jnode->backnext->next = jnode->next;
        }
        inode->mType |= jnode->mType & NT_True;
        free(jnode->links);
        free(jnode->backlinks);
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
  for (node = first; node != NULL; node = node->backnext) {
    Link *links  = node->links;
    int  nlinks  = node->nlinks;
    node->links  = node->backlinks;
    node->nlinks = node->nbacklinks;
    node->backlinks  = links;
    node->nbacklinks = nlinks;
    Node *next     = node->next;
    node->next     = node->backnext;
    node->backnext = next;
    node->aux      = -node->aux;
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
  return (*(Node **) a)->nbacklinks - (*(Node **) b)->nbacklinks;
}

void LatticeLocalOptimization(Node *first, int strictTiming, int trace_flag)
{
//  Node *last, *node, *unsorted;
//  int topolOrd = 1;

  Node *node, *lastnode;
  int i, j, unreachable = 0;
  // For each node, sort links by pointer value to allow
  // for easy comparison whether two nodes have the same set of links

  for (node = first; node != NULL; node = node->next)  {
    node->aux = 0;
    node->backnext = node->next;
    qsort(node->links, node->nlinks, sizeof(Link), lnkcmp);
    qsort(node->backlinks, node->nbacklinks, sizeof(Link), lnkcmp);
  }

// Sort nodes in topological order
// printf("Sorting nodes...\n");

  first->aux = 1;
  for (lastnode = node = first; node != NULL; node = node->next) {
    for (i=0; i < node->nlinks; i++) {
      Node *lnknode = node->links[i].node;
      if (lnknode->aux == 0) {
        for (j=0; j<lnknode->nbacklinks && lnknode->backlinks[j].node->aux==1; j++);
        if (j == lnknode->nbacklinks) {
          lastnode->next = lnknode;
          lastnode  = lnknode;
          lnknode->aux = 1;
          lnknode->next = NULL;
        }
      }
    }
  }
  if (lastnode->nlinks != 0) {
    // There is a cycle in graph so we cannot sort nodes
    // topologicaly, so sort it at least somehow. :o|
    // Anyway this optimization algorithm is not optimal for graphs with cycles.
    for (node = first; node != NULL; node = node->backnext) node->aux = 0;
    first->aux = 1;
    for (lastnode = node = first; node != NULL; node = node->next) {
      for (i=0; i < node->nlinks; i++) {
        Node *lnknode = node->links[i].node;
        if (lnknode->aux == 0) {
          lastnode->next = lnknode;
          lastnode  = lnknode;
          lnknode->aux = 1;
          lnknode->next = NULL;
        }
      }
    }
    for (node=first; node->next->nlinks != 0; node=node->next);

    if (node->next->next) { // Final node is not at the and of chain
      lastnode->next = node->next;
      node->next = node->next->next;
      lastnode = lastnode->next;
      lastnode->next = NULL;
    }
  }

  // !!! Unreachable nodes must be removed before sorting !!!

  for (node=first; node != NULL; node=node->backnext) {
    while (node->backnext && node->backnext->aux == 0) {
      Node *tnode = node->backnext;
      node->backnext = node->backnext->backnext;
      unreachable++;
      free(tnode->links);
      free(tnode->backlinks);
      free(tnode);
    }
  }

//  if (unreachable) Warning("Removing %d unreachable nodes", unreachable);
  if (unreachable) Error("Networks contains unreachable nodes");

  first->backnext = NULL;
  for (node=first; node->next != NULL; node=node->next) {
    node->next->backnext = node;
  }
  for (i=1, node=first; node != NULL; node = node->next, i++) {
    node->aux=i;
  }
  for (;;) {
    if (trace_flag & 2) {
      for (i=0,node=first; node; node=node->next,i++);
      TraceLog("Forward pass.... (number of nodes: %d)", i);
    }
    LatticeLocalOptimization_ForwardPass(first, strictTiming);

    if (trace_flag & 2) {
      for (i=0,node=first; node; node=node->next,i++);
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
    for (node = first; node != NULL; prev = node, node = node->next) {
      if (node->nlinks == 0 || node->nbacklinks == 0 ||
         (node->nlinks == 1 && node->nbacklinks == 1)) {
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
      // link arrays of back-linked nodes to hold node->nlinks more links
      for (j=0; j < node->nbacklinks; j++) {
        Node *backnode = node->backlinks[j].node;
        for (k=0; k<backnode->nlinks && backnode->links[k].node!=node; k++);
        assert(k < backnode->nlinks);
        // Otherwise link to 'node' is missing from which backlink exists
        backnode->links[k] = backnode->links[backnode->nlinks-1];

        backnode->links = (Link *)
          realloc(backnode->links,
                 (backnode->nlinks-1+node->nlinks)*sizeof(Link));
        if (backnode->links == NULL) Error("Insufficient memory");
        backnode->nlinks--;
      }
      // Remove backlinks to current node form linked nodes and realloc
      // backlink arrays of linked nodes to hold node->nbacklinks more backlinks
      for (j=0; j < node->nlinks; j++) {
        Node *forwnode = node->links[j].node;
        for (k=0;k<forwnode->nbacklinks&&forwnode->backlinks[k].node!=node;k++);
        assert(k < forwnode->nbacklinks);
        // Otherwise link to 'node' is missing from which backlink exists
        forwnode->backlinks[k] = forwnode->backlinks[forwnode->nbacklinks-1];

        forwnode->backlinks = (Link *)
          realloc(forwnode->backlinks,
                 (forwnode->nbacklinks-1+node->nbacklinks)*sizeof(Link));
        if (forwnode->backlinks == NULL) Error("Insufficient memory");
        forwnode->nbacklinks--;
      }
      // Alloc new node->nlinks * node->nbacklinks nodes and create new links
      // so that each backlinked node is conected with each linked node through
      // one new node.
      for (i=0; i < node->nlinks; i++) {
        for (j=0; j < node->nbacklinks; j++) {
          Node *tnode;
          Link forwlink = node->links[i];
          Link backlink = node->backlinks[j];

          if ((tnode           = (Node *) calloc(1, sizeof(Node))) == NULL) {
            Error("Insufficient memory");
          }
          *tnode = *node;

          if ((tnode->links     = (Link *) malloc(sizeof(Link))) == NULL ||
             (tnode->backlinks = (Link *) malloc(sizeof(Link))) == NULL) {
            Error("Insufficient memory");
          }
          tnode->nlinks       = 1;
          tnode->nbacklinks   = 1;
          tnode->links[0]     = forwlink;
          tnode->backlinks[0] = backlink;
          forwlink.node->backlinks[forwlink.node->nbacklinks  ].node = tnode;
          forwlink.node->backlinks[forwlink.node->nbacklinks++].like = forwlink.like;
          backlink.node->    links[backlink.node->nlinks  ].node     = tnode;
          backlink.node->    links[backlink.node->nlinks++].like     = backlink.like;
          prev->next = tnode;
          prev = tnode;
        }
      }
      prev->next = node->next;
      free(node->links);
      free(node->backlinks);
      free(node);
      node = prev;
    }
  } while (did_we_clone);

  // Assign to each node unique number, which will later allow to find groups of
  // expanded triphone nodes corresponding to original monophone nodes.
  int nbackmononodes, nforwmononodes, id = 0;
  for (node = first; node != NULL; node = node->next) node->aux = id++;

  // Expand monophone nodes to triphone nodes
  Node *prev = NULL;
  for (node = first; node != NULL; prev = node, node = node->next) {
    ENTRY e, *ep;

    if (node->mType & NT ||
      (node->nlinks == 1 && node->nbacklinks == 1)) {
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
    for (j=0; j < node->nbacklinks; j++) {
      if (node->backlinks[j].node->aux != id) {
        id = node->backlinks[j].node->aux;
        nbackmononodes++;
      }
    }
    // Count groups of linked nodes corresponding to different monophones
    id = -1;
    nforwmononodes = 0;
    for (j=0; j < node->nlinks; j++) {
      if (node->links[j].node->aux != id) {
        id = node->links[j].node->aux;
        nforwmononodes++;
      }
    }

    // Remove links to current node form backlinked nodes and realloc
    // link arrays of backlinked nodes to hold nforwmononodes more links
    for (j=0; j < node->nbacklinks; j++) {
      Node *backnode = node->backlinks[j].node;
      for (k=0; k<backnode->nlinks && backnode->links[k].node!=node; k++);
      assert(k < backnode->nlinks);
      // Otherwise link to 'node' is missing from which backlink exists
      memmove(backnode->links+k, backnode->links+k+1,
              (backnode->nlinks-k-1) * sizeof(Link));

      backnode->links = (Link *)
        realloc(backnode->links,
              (backnode->nlinks-1+nforwmononodes)*sizeof(Link));
      if (backnode->links == NULL) Error("Insufficient memory");
      backnode->nlinks--;
    }

    // Remove backlinks to current node form linked nodes and realloc
    // backlink arrays of linked nodes to hold nbackmononodes more backlinks
    for (j=0; j < node->nlinks; j++) {
      Node *forwnode = node->links[j].node;
      for (k=0;k<forwnode->nbacklinks&&forwnode->backlinks[k].node!=node;k++);
      assert(k < forwnode->nbacklinks);
      // Otherwise link to 'node' is missing from which backlink exists
      memmove(forwnode->backlinks+k, forwnode->backlinks+k+1,
              (forwnode->nbacklinks-k-1) * sizeof(Link));

      forwnode->backlinks = (Link *)
        realloc(forwnode->backlinks,
              (forwnode->nbacklinks-1+nbackmononodes)*sizeof(Link));
      if (forwnode->backlinks == NULL) Error("Insufficient memory");
      forwnode->nbacklinks--;
    }

    // Alloc new nforwmononodes * nbackmononodes nodes and create new links
    // so that each backlinked node is conected through one new node with all
    // linked nodes belonging to one monophone group and vice versa each
    // linked node is conected through one new node with all backlinked nodes
    // belonging to one monophone group
    Link *forwmono_start, *forwmono_end = node->links;
    for (i=0; i < nforwmononodes; i++) {
      for (forwmono_start = forwmono_end;
          forwmono_end < node->links+node->nlinks &&
          forwmono_start->node->aux == forwmono_end->node->aux;
          forwmono_end++);

      assert((i <  nforwmononodes-1 && forwmono_end <  node->links+node->nlinks) ||
            (i == nforwmononodes-1 && forwmono_end == node->links+node->nlinks));

      Link *tlink, *backmono_start, *backmono_end = node->backlinks;
      for (j=0; j < nbackmononodes; j++) {
        for (backmono_start = backmono_end;
          backmono_end < node->backlinks+node->nbacklinks &&
          backmono_start->node->aux == backmono_end->node->aux;
          backmono_end++);

          assert((j <  nbackmononodes-1 && backmono_end <  node->backlinks+node->nbacklinks) ||
                 (j == nbackmononodes-1 && backmono_end == node->backlinks+node->nbacklinks));

        Node *tnode;
        if ((tnode           = (Node *) calloc(1, sizeof(Node))) == NULL) {
          Error("Insufficient memory");
        }
        *tnode = *node;
        tnode->nlinks       = forwmono_end-forwmono_start;
        tnode->nbacklinks   = backmono_end-backmono_start;

        if ((tnode->links =
             (Link *) malloc(tnode->nlinks * sizeof(Link))) == NULL ||
          (tnode->backlinks =
             (Link *) malloc(tnode->nbacklinks * sizeof(Link))) == NULL) {
          Error("Insufficient memory");
        }
        for (tlink = forwmono_start; tlink < forwmono_end; tlink++) {
          tnode->links[tlink-forwmono_start] = *tlink;
          tlink->node->backlinks[tlink->node->nbacklinks  ].node = tnode;
          tlink->node->backlinks[tlink->node->nbacklinks++].like = tlink->like;
        }
        for (tlink = backmono_start; tlink < backmono_end; tlink++) {
          tnode->backlinks[tlink-backmono_start] = *tlink;
          tlink->node->links[tlink->node->nlinks  ].node = tnode;
          tlink->node->links[tlink->node->nlinks++].like = tlink->like;
        }
        prev->next = tnode;
        prev = tnode;
      }
    }
    prev->next = node->next;
    free(node->links);
    free(node->backlinks);
    free(node);
    node = prev;
  }

  // Give triphone names to phone nodes and create hash of these names
  for (node = first; node != NULL; node = node->next) {
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
        lc = lc->nbacklinks ? lc->backlinks[0].node : NULL;
        if (lc == NULL)           break;
        if (!(lc->mType & NT_Phone)) continue;
        if (nonCDphones == NULL)  break;
        e.key  = lc->mpName;
        my_hsearch_r(e, FIND, &ep, nonCDphones);
        if (ep == NULL || !(int) ep->data) break; // Node represents Tee model
      }
      for (rc = node;;) {
        rc = rc->nlinks ? rc->links[0].node : NULL;
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

  for (node = first; node != NULL; node = node->next) {
    for (i=0; i < node->nlinks; i++) {
      if (node->links[i].node == node) {
        if ((tnode           = (Node *) calloc(1, sizeof(Node))) == NULL ||
          (tnode->links     = (Link *) malloc(sizeof(Link))) == NULL ||
          (tnode->backlinks = (Link *) malloc(sizeof(Link))) == NULL) {
          Error("Insufficient memory");
        }

        node->links[i].node = tnode;
        for (j=0; j<node->nbacklinks && node->backlinks[j].node!=node; j++);
        assert(j<node->nbacklinks);
        node->backlinks[j].node = tnode;
        node->backlinks[j].like = 0.0;

        tnode->mType       = NT;
        tnode->pronun     = NULL;
        tnode->nlinks     = 1;
        tnode->nbacklinks = 1;
        tnode->start      = UNDEF_TIME;
        tnode->stop       = UNDEF_TIME;
//        tnode->tokens     = NULL;
//        tnode->exitToken  = NULL;
        tnode->links[0].node     = node;
        tnode->links[0].like     = 0.0;
        tnode->backlinks[0].node = node;
        tnode->backlinks[0].like = node->links[i].like;
        tnode->next = node->next;
        node->next = tnode;
      }
    }
  }
}

// Remove null nones having less than three predecessors or less than three successors
int RemoveRedundantNullNodes(Node *first)
{
  Node *node, *tnode;
  int i, j, node_removed = 0;

  first->backnext = NULL;
  for (node = first; node->next != NULL; node = node->next) {
    node->next->backnext = node;
  }
  for (node = first; node != NULL; node = node->next) {
    if (node->mType & NT && node->pronun == NULL &&
        node->nlinks != 0 && node->nbacklinks != 0  &&
       (node->nlinks == 1 || node->nbacklinks == 1 ||
       (node->nlinks == 2 && node->nbacklinks == 2))) {

     node_removed = 1;

    // Remove links to current node form backlinked nodes and realloc
    // link arrays of backlinked nodes to hold node->nlinks more backlinks
      for (i = 0; i < node->nbacklinks; i++) {
        Node *bakcnode = node->backlinks[i].node;
        for (j=0; j<bakcnode->nlinks && bakcnode->links[j].node!=node; j++);
        assert(j < bakcnode->nlinks); // Otherwise link to 'node' is missing
                                      // from which backlink exists
        bakcnode->links[j] = bakcnode->links[bakcnode->nlinks-1];

        bakcnode->links = (Link *)
          realloc(bakcnode->links,
                (bakcnode->nlinks - 1 + node->nlinks) * sizeof(Link));
        if (bakcnode->links == NULL) Error("Insufficient memory");
        bakcnode->nlinks--;// += word->npronuns-1;
      }

      // Remove backlinks to current node form linked nodes and realloc
      // backlink arrays of linked nodes to hold word->npronuns more backlinks
      for (i=0; i < node->nlinks; i++) {
        Node *forwnode = node->links[i].node;
        for (j=0;j<forwnode->nbacklinks&&forwnode->backlinks[j].node!=node;j++);
        assert(j < forwnode->nbacklinks);
        // Otherwise link to 'node' is missing from which backlink exists
        forwnode->backlinks[j] = forwnode->backlinks[forwnode->nbacklinks-1];

        forwnode->backlinks = (Link *)
          realloc(forwnode->backlinks,
                (forwnode->nbacklinks - 1 + node->nbacklinks) * sizeof(Link));
        if (forwnode->backlinks == NULL) Error("Insufficient memory");
        forwnode->nbacklinks--;
      }
      for (j = 0; j < node->nbacklinks; j++) {
        Node *backnode = node->backlinks[j].node;

        for (i=0; i < node->nlinks; i++) {
          backnode->links[backnode->nlinks  ].node = node->links[i].node;
          backnode->links[backnode->nlinks++].like
            = node->links[i].like + node->backlinks[j].like;
        }
      }
      for (j = 0; j < node->nlinks; j++) {
        Node *forwnode = node->links[j].node;

        for (i=0; i < node->nbacklinks; i++) {
          forwnode->backlinks[forwnode->nbacklinks  ].node = node->backlinks[i].node;
          forwnode->backlinks[forwnode->nbacklinks++].like
            = node->backlinks[i].like + node->links[j].like;
        }
      }
      node->backnext->next = node->next;
      node->next->backnext = node->backnext;
      tnode = node;
      node = node->backnext;
      free(tnode->links);
      free(tnode->backlinks);
      free(tnode);
    }
  }
  return node_removed;
}

struct CorrPhnRec {
  Node      *node;
  long long maxStopTimeTillNow;
  int       id;
};

static int cmp_starts(const void *a, const void *b)
{
  long long diff = ((struct CorrPhnRec *) a)->node->start
                 - ((struct CorrPhnRec *) b)->node->start;

  if (diff != 0)
    return diff;

  return ((struct CorrPhnRec *) a)->node->stop
       - ((struct CorrPhnRec *) b)->node->stop;
}

static int cmp_maxstop(const void *key, const void *elem)
{
  struct CorrPhnRec *corr_phn = (struct CorrPhnRec *) elem;

  if (((Node *) key)->start < corr_phn->maxStopTimeTillNow) {
    if (corr_phn->id == 0 || // first fiead in the array
       ((Node *) key)->start >= (corr_phn-1)->maxStopTimeTillNow)
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

  for (node = first; node != NULL; node = node->next) {
    if (node->mType & NT_Phone && node->mType & NT_True) ncorr_phns++;
  }
  if (ncorr_phns == 0) Error("No correct phoneme node in network");

  corr_phn = (struct CorrPhnRec *)malloc(sizeof(struct CorrPhnRec)*ncorr_phns);
  if (corr_phn == NULL) Error("Insufficient memory");

  for (node = first; node != NULL; node = node->next) {
    if (node->mType & NT_Phone && node->mType & NT_True) {
      corr_phn[i++].node = node;
    }
  }

  qsort(corr_phn, ncorr_phns, sizeof(struct CorrPhnRec), cmp_starts);

  maxStopTime = corr_phn[0].node->stop;
  for (i = 0; i < ncorr_phns; i++) {
    corr_phn[i].id = i;
    maxStopTime = HIGHER_OF(maxStopTime, corr_phn[i].node->stop);
    corr_phn[i].maxStopTimeTillNow = maxStopTime;
  }
  for (node = first; node != NULL; node = node->next) {
    if (!(node->mType & NT_Phone)) continue;

    if (node->stop  <= node->start ||
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
              overlaped->node->start < node->stop; overlaped++) {
          if (overlaped->node->stop  <= overlaped->node->start ||
             overlaped->node->stop  <= node->start) continue;

          node->phoneAccuracy =
            HIGHER_OF(node->phoneAccuracy,
                      (SamePhoneme(overlaped->node->mpName, node->mpName) + 1.0) *
                      (LOWER_OF(overlaped->node->stop, node->stop) -
                        HIGHER_OF(overlaped->node->start, node->start)) /
                      (overlaped->node->stop - overlaped->node->start) - 1.0);

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

  if (node_hash->tabsize == 0 && !my_hcreate_r(INIT_NODE_HASH_SIZE, node_hash)) {
    Error("Insufficient memory");
  }
  e.key = node_id;
  my_hsearch_r(e, FIND, &ep, node_hash);

  if (ep != NULL) return (Node *) ep->data;


  node = (Node *) calloc(1, sizeof(Node));
  if (node == NULL) Error("Insufficient memory");

  node->nlinks     = 0;
  node->links      = NULL;
  node->nbacklinks = 0;
  node->backlinks  = NULL;
  node->mType       = NT;
  node->start      = UNDEF_TIME;
  node->stop       = UNDEF_TIME;
  node->pronun     = NULL;
  node->backnext   = *last;
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
    node->start = node->stop = UNDEF_TIME;
    node->phoneAccuracy = 1.0;
//    node->aux = *totalNumOfNodes;
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
      subnetFirst->nbacklinks = node->nbacklinks;
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
        node->pronun = word ? word->pronuns[pronunVar] : NULL;
      } else {
        node->pronun = NULL;
      }
    }
    if (nodeType != 'S') {
      if (fscanf(lfp, " (%lld %lld)", &start, &stop)==2 && !(labelFormat.TIMES_OFF)) {
        long center_shift = labelFormat.CENTRE_TM ? sampPeriod / 2 : 0;
        node->start = (start - center_shift - labelFormat.left_extent)  / sampPeriod;
        node->stop  = (stop  + center_shift + labelFormat.right_extent) / sampPeriod;
      }
    }
    if (fscanf(lfp, "%d ", &numOfLinks) != 1) {
      Error("Invalid definition of node %d in file %s.\n"
            "Number of links is expected", nodeId, file_name);
    }
    if (nodeType == 'S') { // Add links to the final node of the subnetwork
      while (node->next != NULL) node = node->next;
    }
    if (numOfLinks) {
      if ((node->links = (Link *) malloc(numOfLinks * sizeof(Link))) == NULL) {
        Error("Insufficient memory");
      }
    } else {
      if (nodeType == 'M') {
        Error("Invalid definition of node %d in file %s.\n"
              "Model node must have at least one link", nodeId, file_name);
      }
      node->links = NULL;
    }
    node->nlinks = numOfLinks;

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
      node->links[j].node = nodes[linkId];
      node->links[j].like = linkLike;
      nodes[linkId]->nbacklinks++;
    }
  }
  for (i = 1; i < numOfNodes-1; i++) {
    if (nodes[i]->nlinks == 0) {
      if (nodes[numOfNodes-1]->nlinks == 0) {
        Error("Network contains multiple nodes with no successors (%s)",
              file_name);
      }
      node = nodes[numOfNodes-1];
      nodes[numOfNodes-1] = nodes[i];
      nodes[i] = node;
    }
    if (nodes[i]->nbacklinks == 0) {
      if (nodes[0]->nbacklinks == 0) {
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
  if (nodes[0]->nbacklinks != 0 || nodes[numOfNodes-1]->nlinks != 0) {
    Error("Network contain no start node or no final node (%s)", file_name);
  }
  if (!(nodes[0]           ->mType & NT) || nodes[0]           ->pronun != NULL ||
     !(nodes[numOfNodes-1]->mType & NT) || nodes[numOfNodes-1]->pronun != NULL) {
    Error("Start node and final node must be Null nodes (%s)", file_name);
  }
  for (i = 0; i < numOfNodes-1; i++) {
    nodes[i]->next = nodes[i+1];
  }

  // create back links
  for (i = 0; i < numOfNodes; i++) {
    if (!nodes[i]->backlinks) // Could be already alocated for subnetwork
      nodes[i]->backlinks = (Link *) malloc(nodes[i]->nbacklinks * sizeof(Link));
    if (nodes[i]->backlinks == NULL) Error("Insufficient memory");
    nodes[i]->nbacklinks = 0;
  }
  for (i = 0; i < numOfNodes; i++) {
    for (j=0; j < nodes[i]->nlinks; j++) {
      Node *forwNode = nodes[i]->links[j].node;
      forwNode->backlinks[forwNode->nbacklinks].node = nodes[i];
      forwNode->backlinks[forwNode->nbacklinks].like = nodes[i]->links[j].like;
      forwNode->nbacklinks++;
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

          if (node_hash.tabsize == 0 && !my_hcreate_r(nnodes, &node_hash)) {
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
              node->start = 100 * (long long) (0.5 + 1e5 *
                            getFloat(valptr, &chptr, file_name, line_no));
            }
            valptr = colonptr+1;
          }
          node->stop = 100 * (long long) (0.5 + 1e5 *
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
            node->pronun = word->pronuns[pron_var-1];
          }
        }
      }
      if (state == ARC_DEF) {
        if (!*chptr || !strcmp(chptr, "END") || !strcmp(chptr, "E")) {
          if (getHTKstr(valptr, &chptr)) {
            Error("%s (%s:%d)", chptr, file_name, line_no);
          }
          enode = find_or_create_node(&node_hash, valptr, &last);

          int nl = ++node->nlinks;
          node->links = (Link *) realloc(node->links, nl * sizeof(Link));
          if (node->links == NULL) Error("Insufficient memory");
          node->links[nl-1].node = enode;
          node->links[nl-1].like = 0.0;

          nl = ++enode->nbacklinks;
          enode->backlinks = (Link *) realloc(enode->backlinks, nl * sizeof(Link));
          if (enode->backlinks == NULL) Error("Insufficient memory");
          enode->backlinks[nl-1].node = node;
          enode->backlinks[nl-1].like = 0.0;

        } else if (!strcmp(chptr, "language") || !strcmp(chptr, "l")) {
          FLOAT like = getFloat(valptr, &chptr, file_name, line_no);
          //Set LM score to link pointing to enode. This link can possibly start
          //from a phone node already inserted (div=) between 'node' ans 'enode'
          Node *last = enode->backlinks[enode->nbacklinks-1].node;
          last->links[last->nlinks-1].like = like;
          enode->backlinks[enode->nbacklinks-1].like = like;
        } else if (!strcmp(chptr, "div") || !strcmp(chptr, "d")) {
          ENTRY e, *ep;
          char  name[1024];
          float time;
          int   n;
          Node  *last = node;
          FLOAT like  = node->links[node->nlinks-1].like;

          if (node->links[node->nlinks-1].node != enode) {
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
              (tnode->links     = (Link *) malloc(sizeof(Link))) == NULL ||
              (tnode->backlinks = (Link *) malloc(sizeof(Link))) == NULL) {
              Error("Insufficient memory");
            }

            //Use special type to mark nodes inserted by d=..., they will need
            //special treatment. Later, they will become ordinary NT_Phone nodes
            tnode->mType      = NT_Phone | NT_Model;
            tnode->nlinks    = tnode->nbacklinks = 1;
            tnode->backnext  = enode->backnext;
            enode->backnext  = tnode;
            tnode->phoneAccuracy = 1.0;
            //Store phone durations now. Will be replaced by absolute times below.
            tnode->start  = time != -FLT_MAX ? 100 * (long long) (0.5 + 1e5 * time) : UNDEF_TIME;
            tnode->stop   = UNDEF_TIME;
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
            last->links[last->nlinks-1].node = tnode;
            last->links[last->nlinks-1].like = 0.0;
            tnode->backlinks[0].node = last;
            tnode->backlinks[0].like = 0.0;
            last = tnode;
          }
          if (strcmp(phn_marks,":")) {
            Error("Invalid specification of phone marks (d=) (%s:%d)",
                  file_name, line_no);
          }
          last->links[last->nlinks-1].node = enode;
          last->links[last->nlinks-1].like = like;
          enode->backlinks[enode->nbacklinks-1].node = last;
          enode->backlinks[enode->nbacklinks-1].like = like;

        } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
          Error("%s (%s:%d)", chptr, file_name, line_no);
        }
      }
    }
  }
  my_hdestroy_r(&node_hash, 1);
  lnode = last;
  first = last = NULL;
  if (lnode) lnode->next = NULL;
  for (node = lnode; node != NULL; fnode = node, node = node->backnext) {
    if (node->backnext) node->backnext->next = node;

    if (node->nlinks == 0) {
      if (last)
        Error("Network has multiple nodes with no successors (%s)", file_name);
      last = node;
    }
    if (node->nbacklinks == 0) {
      if (first)
        Error("Network has multiple nodes with no predecessor (%s)", file_name);
      first = node;
    }
    //If only stop time is specified, set start time to lowest predecessor stop time
    if (node->start == UNDEF_TIME && node->stop != UNDEF_TIME) {
      int i;
      for (i = 0; i < node->nbacklinks; i++) {
        Node *backnode = node->backlinks[i].node;
        //When seraring predecessors, skip nodes inserted by d=...
        while (backnode->mType == (NT_Phone | NT_Model)) {
          assert(backnode->nbacklinks == 1);
          backnode = backnode->backlinks[0].node;
        }
        if (backnode->stop != UNDEF_TIME) {
          node->start = node->start == UNDEF_TIME
                        ?          backnode->stop
                        : LOWER_OF(backnode->stop, node->start);
        }
      }
      if (node->start == UNDEF_TIME) node->start = 0;
    }
    //For model nodes defined by d=... (NT_Phone | NT_Model), node->start contains
    //only phone durations. Absolute times must be computed derived starting from
    //the end time of the node to which arc with d=... definition points.
    if (node->mType == (NT_Phone | NT_Model)) {
      assert(node->nlinks == 1);
      node->stop = node->links[0].node->mType == (NT_Phone | NT_Model)
                   && node->start != UNDEF_TIME
                   ? node->links[0].node->start : node->links[0].node->stop;
      node->start = node->start != UNDEF_TIME && node->stop != UNDEF_TIME
                    ? node->stop - node->start : node->links[0].node->start;
    }
  }
  if (!first || !last) {
    Error("Network contain no start node or no final node (%s)", file_name);
  }
  if (first != fnode) {
    if (first->next)     first->next->backnext = first->backnext;
    if (first->backnext) first->backnext->next = first->next;
    if (first == lnode)  lnode = first->backnext;

    first->backnext = NULL;
    fnode->backnext = first;
    first->next = fnode;
  }
  if (last != lnode) {
    if (last->next)     last->next->backnext = last->backnext;
    if (last->backnext) last->backnext->next = last->next;
    last->next = NULL;
    lnode->next = last;
    last->backnext = lnode;
  }
  for (node = first; node != NULL; node = node->next) {
    if (node->mType == (NT_Phone | NT_Model)) {
      node->mType = NT_Phone;
    }
    if (node->start != UNDEF_TIME) {
      node->start = (node->start - labelFormat.left_extent) / sampPeriod;
    }
    if (node->stop  != UNDEF_TIME) {
      node->stop  = (node->stop + labelFormat.right_extent) / sampPeriod;
    }
  }
  if (first->pronun != NULL) {
    node = (Node *) calloc(1, sizeof(Node));
    if (node == NULL) Error("Insufficient memory");
    node->next       = first;
    node->backnext   = NULL;
    first->backnext  = node;
    node->mType       = NT;
    node->pronun     = NULL;
    node->start      = UNDEF_TIME;
    node->stop       = UNDEF_TIME;
    node->nbacklinks = 0;
    node->backlinks  = NULL;
    node->nlinks     = 1;
    node->links      = (Link*) malloc(sizeof(Link));
    if (node->links == NULL) Error("Insufficient memory");
    node->links[0].like = 0.0;
    node->links[0].node = first;
    first->nbacklinks = 1;
    first->backlinks  = (Link*) malloc(sizeof(Link));
    if (first->backlinks == NULL) Error("Insufficient memory");
    first->backlinks[0].like = 0.0;
    first->backlinks[0].node = node;
    first = node;
  }
  if (last->pronun != NULL) {
    node = (Node *) calloc(1, sizeof(Node));
    if (node == NULL) Error("Insufficient memory");
    last->next      = node;
    node->next      = NULL;
    node->backnext  = last;
    node->mType      = NT;
    node->pronun    = NULL;
    node->start     = UNDEF_TIME;
    node->stop      = UNDEF_TIME;
    node->nlinks    = 0;
    node->links     = NULL;
    last->nlinks = 1;
    last->links  = (Link*) malloc(sizeof(Link));
    if (last->links == NULL) Error("Insufficient memory");
    last->links[0].like = 0.0;
    last->links[0].node = node;
    node->nbacklinks = 1;
    node->backlinks  = (Link*) malloc(sizeof(Link));
    if (node->backlinks == NULL) Error("Insufficient memory");
    node->backlinks[0].like = 0.0;
    node->backlinks[0].node = last;
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

  for (node = net; node != NULL; node = node->next) {
    node->estate_id = i++;
  }
  fprintf(fp, "digraph \"dnet%d\" {\nrankdir=LR\n", dnetcnt++);


  for (node = net; node != NULL; node = node->next) {
    fprintf(fp, "n%d [shape=%s,label=\"%d:%s", node->estate_id,
            node->mType & NT ? "box" : "ellipse", node->estate_id,
            node->mType & NT ? (node->pronun ?
                                    node->pronun->word->mpName : "-"):
            node->mType & NT_Phone? node->mpName :
            node->mType & NT_Model? node->hmm->mpMacro->mpName : "???");

    if (node->mType & NT && node->pronun != NULL) {
      if (node->pronun != node->pronun->word->pronuns[0]) {
        fprintf(fp, ":%d", node->pronun->variant_no);
      }
      fprintf(fp, "\\n");

      if (node->pronun->outSymbol != node->pronun->word->mpName) {
        fprintf(fp, "[%s]", node->pronun->outSymbol ?
                            node->pronun->outSymbol : "");
      }
      if (node->pronun->prob != 0.0) {
        fprintf(fp, " "FLOAT_FMT, node->pronun->prob);
      }
    }
    fprintf(fp, "\"];\n");
//    if (node->next != NULL) {
//     fprintf(fp,"n%d -> n%d [color=black,weight=1]\n",
//             node->estate_id, node->next->estate_id);
//    }
//    if (node->backnext != NULL) {
//     fprintf(fp,"n%d -> n%d [color=gray,weight=1]\n",
//             node->estate_id, node->backnext->estate_id);
//    }
  }
  for (node = net; node != NULL; node = node->next) {
    for (i = 0; i < node->nlinks; i++) {
      fprintf(fp,"n%d -> n%d [color=blue,weight=1",
              node->estate_id,node->links[i].node->estate_id);
      if (node->links[i].like != 0.0) {
        fprintf(fp,",label=\""FLOAT_FMT"\"", node->links[i].like);
      }
      fprintf(fp,"];\n");
    }
//    for (i = 0; i < node->nbacklinks; i++) {
//      fprintf(fp,"n%d -> n%d [color=red,weight=1",
//              node->estate_id,node->backlinks[i].node->estate_id);
//      if (node->backlinks[i].like != 0.0) {
//        fprintf(fp,",label=\""FLOAT_FMT"\"", node->backlinks[i].like);
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
            last->next = nodes[node_id];
          }
          last = nodes[node_id];
          nodes[node_id]->next = NULL;

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
      nodes[node_id]->nlinks     = 0;
      nodes[node_id]->nbacklinks = 0;
      nodes[node_id]->links      = NULL;
      nodes[node_id]->backlinks  = NULL;
      if (node_word && node_word->npronuns <= node_var) {
        Error("Word %s does not have pronunciation varian %d (%s:%d)",
              node_word->mpName, node_var+1, file_name, line_no);
      }
      nodes[node_id]->pronun     = node_word ? node_word->pronuns[node_var]
                                             : NULL;
      nodes[node_id]->start     = UNDEF_TIME;
      nodes[node_id]->stop      = node_time + labelFormat.right_extent;
    } else if (state == ARC_DEF) {
      if (arc_start == -1 || arc_end == -1) {
        Error("Start node or end node not defined (%s:%d)", file_name, line_no);
      }

      int linkId = nodes[arc_start]->nlinks++;
      nodes[arc_start]->links =
        (Link *) realloc(nodes[arc_start]->links, (linkId+1) * sizeof(Link));

      if (nodes[arc_start]->links == NULL) Error("Insufficient memory");

      last = nodes[arc_start];

      if (phn_marks) {
        ENTRY e, *ep;
        char  name[1024];
        float time;
        int   n;

        while (sscanf(phn_marks, ":%[^,],%f%n", name, &time, &n) > 1) {
          phn_marks+=n;

          if ((node            = (Node *) calloc(1, sizeof(Node))) == NULL ||
             (node->links     = (Link *) malloc(sizeof(Link))) == NULL ||
             (node->backlinks = (Link *) malloc(sizeof(Link))) == NULL) {
            Error("Insufficient memory");
          }
          node->mType = NT_Phone;
          node->nlinks = node->nbacklinks = 1;
          node->next   = last->next;
          last->next = node;
          node->phoneAccuracy = 1.0;

          if (!(labelFormat.TIMES_OFF)) {
            node->start = last->stop - labelFormat.left_extent - labelFormat.right_extent;
            node->stop  = last->stop + 100 * (long long) (0.5 + 1e5 * time)
                                     + labelFormat.right_extent;
          } else {
            node->start  = UNDEF_TIME;
            node->stop   = UNDEF_TIME;
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
          last->links[last->nlinks-1].node = node;
          last->links[last->nlinks-1].like = 0.0;
          node->backlinks[0].node = last;
          node->backlinks[0].like = 0.0;
          last = node;
        }
        if (strcmp(phn_marks,":")) {
          Error("Invalid specification of phone marks (d=) (%s:%d)",
                file_name, line_no);
        }
      }
      linkId = nodes[arc_end]->nbacklinks++;
      nodes[arc_end]->backlinks =
        (Link *) realloc(nodes[arc_end]->backlinks, (linkId+1) * sizeof(Link));

      if (nodes[arc_end]->backlinks == NULL) Error("Insufficient memory");

      last->links[last->nlinks-1].node = nodes[arc_end];
      last->links[last->nlinks-1].like = arc_like;
      nodes[arc_end]->backlinks[linkId].node = last;
      nodes[arc_end]->backlinks[linkId].like = arc_like;

      if (nodes[arc_start]->stop != UNDEF_TIME) {
        if (nodes[arc_end]->start == UNDEF_TIME) {
          nodes[arc_end]->start = nodes[arc_start]->stop - labelFormat.left_extent
                                                         - labelFormat.right_extent;
        } else {
          nodes[arc_end]->start = LOWER_OF(nodes[arc_end]->start,
                                           nodes[arc_start]->stop - labelFormat.left_extent
                                                                  - labelFormat.right_extent);
        }
      }
    }
  }
  free(nodes);
  fnode = first;
  first = last = NULL;
  if (fnode) fnode->backnext = NULL;
  for (node = fnode; node != NULL; lnode = node, node = node->next) {
    if (node->next) node->next->backnext = node;

    if (node->nlinks == 0) {
      if (last)
        Error("Network has multiple nodes with no successors (%s)", file_name);
      last = node;
    }
    if (node->nbacklinks == 0) {
      if (first)
        Error("Network has multiple nodes with no predecessor (%s)", file_name);
      first = node;
    }
  }
  if (!first || !last) {
    Error("Network contain no start node or no final node (%s)", file_name);
  }
  if (first != fnode) {
    if (first->next) first->next->backnext = first->backnext;
    first->backnext->next = first->next;
    first->backnext = NULL;
    fnode->backnext = first;
    first->next = fnode;
  }
  if (last != lnode) {
    if (last->backnext) last->backnext->next = last->next;
    last->next->backnext = last->backnext;
    last->next = NULL;
    lnode->next = last;
    last->backnext = lnode;
  }
  if (first->pronun != NULL) {
    node = (Node *) calloc(1, sizeof(Node));
    if (node == NULL) Error("Insufficient memory");
    node->next       = first;
    node->backnext   = NULL;
    first->backnext  = node;
    node->mType       = NT;
    node->pronun     = NULL;
    node->start      = UNDEF_TIME;
    node->stop       = UNDEF_TIME;
    node->nbacklinks = 0;
    node->backlinks  = NULL;
    node->nlinks     = 1;
    node->links      = (Link*) malloc(sizeof(Link));
    if (node->links == NULL) Error("Insufficient memory");
    node->links[0].like = 0.0;
    node->links[0].node = first;

    first->nbacklinks = 1;
    first->backlinks  = (Link*) malloc(sizeof(Link));
    if (first->backlinks == NULL) Error("Insufficient memory");
    first->backlinks[0].like = 0.0;
    first->backlinks[0].node = node;
    first = node;
  }
  if (last->pronun != NULL) {
    node = (Node *) calloc(1, sizeof(Node));
    if (node == NULL) Error("Insufficient memory");
    last->next      = node;
    node->next      = NULL;
    node->backnext  = last;

    node->mType      = NT;
    node->pronun    = NULL;
    node->start     = UNDEF_TIME;
    node->stop      = UNDEF_TIME;
    node->nlinks    = 0;
    node->links     = NULL;

    last->nlinks = 1;
    last->links  = (Link*) malloc(sizeof(Link));
    if (last->links == NULL) Error("Insufficient memory");
    last->links[0].like = 0.0;
    last->links[0].node = node;

    node->nbacklinks = 1;
    node->backlinks  = (Link*) malloc(sizeof(Link));
    if (node->backlinks == NULL) Error("Insufficient memory");
    node->backlinks[0].like = 0.0;
    node->backlinks[0].node = last;
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

  for (i = 0, node = first; node != NULL; node = node->next, i++)  {
    node->aux = i;
  }
  fprintf(lfp,"NUMNODES: %d\n", i);
  for (i = 0, node = first; node != NULL; node = node->next, i++) {
    int j,
    type = node->mType & NT_Model       ? 'M'  :
           node->mType & NT_Phone       ? 'M'  :
           node->mType & NT_Subnet      ? 'S'  :
           node->mType & NT        ?
             (node->pronun == NULL     ?
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
            node->mType & NT_Model   ? node->hmm->mpMacro->mpName :
            node->mType & NT_Phone   ? node->mpName :
            node->mType & NT_Subnet  ? node->mpName :
            node->mType & NT    ?
              (node->pronun == NULL ? "-" :
                                      node->pronun->word->mpName):
                                      "?");
    if (node->mType & NT && node->pronun) {
       if (node->pronun->word->mpName != node->pronun->outSymbol) {
         fprintf(lfp," [%s]", node->pronun->outSymbol);
       }
       if (node->pronun->prob != 0.0 || node->pronun->word->npronuns > 1) {
         fprintf(lfp," {%d "FLOAT_FMT"}",
                 node->pronun->variant_no,
                 node->pronun->prob);
       }
    }
    if (node->mType & NT_Phone && node->phoneAccuracy != 1.0) {
      fprintf(lfp," {"FLOAT_FMT"}", node->phoneAccuracy);
    }
    if (!(labelFormat.TIMES_OFF) &&
       node->start != UNDEF_TIME && node->stop != UNDEF_TIME) {
      int ctm = labelFormat.CENTRE_TM;
      fprintf(lfp," (%lld %lld)",
                  (long long) sampPeriod * (2 * node->start + ctm) / 2 - labelFormat.left_extent,
                  (long long) sampPeriod * (2 * node->stop - ctm)  / 2 + labelFormat.right_extent);
    }
    fprintf(lfp,"\t%d", node->nlinks);
    for (j = 0; j < node->nlinks; j ++) {
      fprintf(lfp," %d",node->links[j].node->aux);
      if (node->links[j].like != 0.0) {
        fprintf(lfp," {"FLOAT_FMT"}", node->links[j].like);
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
  if (expOptions.no_word_expansion &&  expOptions.CD_phone_expansion
  && expOptions.no_optimization   && !out_net_fmt.no_LM_likes
  &&!out_net_fmt.no_times         && !out_net_fmt.no_word_nodes &&
    !out_net_fmt.no_model_nodes   && !out_net_fmt.no_pronun_vars) return;
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
