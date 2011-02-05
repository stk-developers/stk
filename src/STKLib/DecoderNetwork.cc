
#include "DecoderNetwork.h"
#include <cstdarg>

namespace STK
{
  


  //***************************************************************************
  //***************************************************************************
//  int nbacklinkscmp(const void *a, const void *b) 
//  {
//    return (*(DecoderNetwork::NodeType **) a)->rNBackLinks() 
//         - (*(DecoderNetwork::NodeType **) b)->rNBackLinks();
//  }
  
  
  
  //***************************************************************************

  

  //***************************************************************************
  //***************************************************************************
  int 
  DecoderNetwork::
  lnkcmp(const void *a, const void *b)
  {
  //  return ((LinkType *) a)->pNode() - ((LinkType *) b)->pNode();
  //  Did not work with gcc, probably bug in gcc pointer arithmetic
    return (char *)((LinkType *) a)->pNode() - (char *)((LinkType *) b)->pNode();
  }


  //***************************************************************************
  //***************************************************************************
    void 
    DecoderNetwork::
    BuildFromLabels(const Label* pLabels, NodeKind nodeKind)
    {
      const Label*   p_lp;
      NodeType*      p_first;
      NodeType*      p_last = &(*end());
      NodeType*      p_node;
    
      // allocate the memory
      if ((p_first             = (NodeType *) calloc(1, sizeof(NodeType))) == NULL ||
          (p_last              = (NodeType *) calloc(1, sizeof(NodeType))) == NULL ||
          (p_first->rpLinks()    = (LinkType *) malloc(sizeof(LinkType))) == NULL    ||
          (p_last->rpBackLinks() = (LinkType *) malloc(sizeof(LinkType))) == NULL) 
      {
        Error("Insufficient memory");
      }

      p_first->mC.mpName       = p_last->mC.mpName      = NULL;     
      p_first->mC.mType        = p_last->mC.mType       = NT_WORD;
      p_first->mC.mpPronun     = p_last->mC.mpPronun    = NULL;
      p_first->rNLinks()      = p_last->rNBackLinks() = 1;
      p_first->rNBackLinks()  = p_last->rNLinks()     = 0;
      p_first->rpBackLinks()  = p_last->rpLinks()     = NULL;
      p_first->mC.SetStart(UNDEF_TIME);
      p_first->mC.SetStop(UNDEF_TIME);
      p_last->mC.SetStop(UNDEF_TIME);
      p_last->mC.SetStart(UNDEF_TIME);
      //p_first->mC.mpTokens        = p_last->mC.mpTokens        = NULL;
      //p_first->mC.mpExitToken     = p_last->mC.mpExitToken     = NULL;
    
      p_node = p_first;
      
      for (p_lp = pLabels; p_lp != NULL; p_lp = p_lp->mpNext) 
      {
        NodeType*  p_tmp_node;
    
        if ((p_tmp_node              = (NodeType *) calloc(1, sizeof(NodeType))) == NULL ||
            (p_tmp_node->rpLinks()     = (LinkType *) malloc(sizeof(LinkType))) == NULL    ||
            (p_tmp_node->rpBackLinks() = (LinkType *) malloc(sizeof(LinkType))) == NULL) 
        {
          Error("Insufficient memory");
        }
        
        p_tmp_node->mC.mpName = NULL;
        
        p_node->rpLinks()[0].SetNode(p_tmp_node);
        p_node->rpLinks()[0].SetLmLike(0.0);
        p_node->rpLinks()[0].SetAcousticLike(0.0);
        
        switch (nodeKind) 
        {
          case NT_WORD:  p_tmp_node->mC.mpPronun = ((Word *) p_lp->mpData)->pronuns[0]; break;
          case NT_MODEL: p_tmp_node->mC.mpHmm    =   (Hmm *) p_lp->mpData;              break;
          case NT_PHONE: p_tmp_node->mC.mpName   =  (char *) p_lp->mpData;              break;
          default:       Error("Fatal: Invalid node type");
        }
        
        p_tmp_node->mC.mType       = nodeKind;
        p_tmp_node->rNLinks()     = 1;
        p_tmp_node->rNBackLinks() = 1;
        p_tmp_node->mC.SetStart(p_lp->mStart);
        p_tmp_node->mC.SetStop(p_lp->mStop);
        p_tmp_node->rpBackLinks()[0].SetNode(p_node);
        p_tmp_node->rpBackLinks()[0].SetLmLike(0.0);
        p_tmp_node->rpBackLinks()[0].SetAcousticLike(0.0);
        p_node->mpNext = p_tmp_node;
        p_node = p_tmp_node;
      }
      
      p_node->mpNext = p_last;
      p_node->rpLinks()[0].SetNode(p_last);
      p_node->rpLinks()[0].SetLmLike(0.0);
      p_node->rpLinks()[0].SetAcousticLike(0.0);
      p_last->rpBackLinks()[0].SetNode(p_node);
      p_last->rpBackLinks()[0].SetLmLike(0.0);
      p_last->rpBackLinks()[0].SetAcousticLike(0.0);
      p_last->mpNext = &(*end());

      this->splice(p_first, p_last);
    }
  // Network::MakeNetworkFromLabels(Label * pLabels, NodeType nodeType)
  //***************************************************************************
  



  //****************************************************************************
  //****************************************************************************
    void 
    DecoderNetwork::
    ExpandByDictionary(
      MyHSearchData* pDict,
      bool keep_word_nodes,
      bool multiple_pronun)
    {
      NodeType* node; 
      NodeType* prev = NULL;
      int   i;
      int   j;
    
      Pronun    singlePronun;
      Pronun*   singlePronunPtr;
      Word      singlePronunWrd;

      singlePronunWrd.npronuns = 1;
      singlePronunWrd.pronuns  = &singlePronunPtr;
      singlePronunPtr = &singlePronun;
    
      assert(!IsEmpty() || pFirst()->mC.mType & NT_WORD || pFirst()->mC.mpPronun == NULL);
    
      for (node = pFirst(); node != NULL; prev = node, node = node->mpNext) 
      {
        if (!(node->mC.mType & NT_WORD)) continue;
    
        if (node->mC.mpPronun == NULL) continue;
        Word *word = node->mC.mpPronun->mpWord;
    
        //Do not expand non-dictionary words, which where added by ReadSTKNetwork
        if (word->npronunsInDict == 0) continue;
    
        if (!multiple_pronun) 
        {
          singlePronunWrd.mpName = node->mC.mpPronun->mpWord->mpName;
          word = &singlePronunWrd;
          *word->pronuns = node->mC.mpPronun;
        }
    
        // Remove links to current node form backlinked nodes and realloc
        // link arrays of backlinked nodes to hold word->npronuns more backlinks
        for (i = 0; i < node->rNBackLinks(); i++) 
        {
          NodeType *bakcnode = node->rpBackLinks()[i].pNode();

          for (j=0; j<bakcnode->NLinks() && bakcnode->rpLinks()[j].pNode()!=node; j++)
          {}

          assert(j < bakcnode->NLinks()); // Otherwise link to 'node' is missing
                                        // from which backlink exists
          bakcnode->rpLinks()[j] = bakcnode->rpLinks()[bakcnode->NLinks()-1];
    
          bakcnode->rpLinks() = (LinkType *)
            realloc(bakcnode->rpLinks(),
                  (bakcnode->NLinks() - 1 + word->npronuns) * sizeof(LinkType));
          if (bakcnode->rpLinks() == NULL) Error("Insufficient memory");
          bakcnode->rNLinks()--;// += word->npronuns-1;
        }
    
        // Remove backlinks to current node form linked nodes and realloc
        // backlink arrays of linked nodes to hold word->npronuns more backlinks
        for (i=0; i < node->NLinks(); i++) 
        {
          NodeType* forwnode = node->rpLinks()[i].pNode();

          for (j=0;j<forwnode->rNBackLinks()&&forwnode->rpBackLinks()[j].pNode()!=node;j++)
          {}

          assert(j < forwnode->rNBackLinks());
          // Otherwise link to 'node' is missing from which backlink exists
          forwnode->rpBackLinks()[j] = forwnode->rpBackLinks()[forwnode->rNBackLinks()-1];
    
          forwnode->rpBackLinks() = (LinkType *)
            realloc(forwnode->rpBackLinks(),
                  (forwnode->rNBackLinks() - 1 + word->npronuns) * sizeof(LinkType));

          if (forwnode->rpBackLinks() == NULL) 
            Error("Insufficient memory");

          forwnode->rNBackLinks()--;
        }

        for (i = 0; i < word->npronuns; i++) 
        {
          Pronun *pronun = word->pronuns[i];
          NodeType *pronun_first = NULL, *pronun_prev = NULL, *tnode;
    
          for (j = 0; j < pronun->nmodels; j++) 
          {
            tnode = (NodeType *) calloc(1, sizeof(NodeType));
            
            if (tnode == NULL) 
              Error("Insufficient memory");
    
            tnode->mC.mType       = NT_PHONE | (node->mC.mType & NT_TRUE);
            tnode->mC.mpName      = pronun->model[j].mpName;
            tnode->mC.SetStart(node->mC.Start());
            tnode->mC.SetStop (node->mC.Stop());
            tnode->mC.SetPhoneAccuracy(1.0);
    
            if (j == 0) 
            {
              pronun_first = tnode;
            } 
            else 
            {
              if ((pronun_prev->rpLinks() = (LinkType *) malloc(sizeof(LinkType))) == NULL ||
                  (tnode->rpBackLinks()   = (LinkType *) malloc(sizeof(LinkType))) == NULL) 
              {
                Error("Insufficient memory");
              }
              
              tnode->rNBackLinks()              = 1;
              tnode->rpBackLinks()[0].SetNode(pronun_prev);
              tnode->rpBackLinks()[0].SetLmLike(0.0);
              tnode->rpBackLinks()[0].SetAcousticLike(0.0);
              pronun_prev->rNLinks()            = 1;
              pronun_prev->rpLinks()[0].SetNode(tnode);
              pronun_prev->rpLinks()[0].SetLmLike(0.0);
              pronun_prev->rpLinks()[0].SetAcousticLike(0.0);
              pronun_prev->mpNext             = tnode;
            }
            pronun_prev = tnode;
          }
          
          if (keep_word_nodes || j == 0) 
          {
            tnode = (NodeType *) calloc(1, sizeof(NodeType));
            
            if (tnode == NULL) Error("Insufficient memory");
    
            tnode->mC.mpName   = NULL;
            tnode->mC.mType    = NT_WORD | (node->mC.mType & NT_TRUE);
            tnode->mC.mpPronun = keep_word_nodes ? word->pronuns[i] : NULL;
            tnode->mC.SetStart (node->mC.Start());
            tnode->mC.SetStop  (node->mC.Stop());
    
            if (j == 0) {
              pronun_first = tnode;
            } else {
              if ((pronun_prev->rpLinks() = (LinkType *) malloc(sizeof(LinkType))) == NULL ||
                (tnode->rpBackLinks()   = (LinkType *) malloc(sizeof(LinkType))) == NULL) {
                Error("Insufficient memory");
              }
              tnode->rNBackLinks()          = 1;
              tnode->rpBackLinks()[0].SetNode(pronun_prev);
              tnode->rpBackLinks()[0].SetLmLike(0.0);
              tnode->rpBackLinks()[0].SetAcousticLike(0.0);
              pronun_prev->rNLinks()        = 1;
              pronun_prev->rpLinks()[0].SetNode(tnode);
              pronun_prev->rpLinks()[0].SetLmLike(0.0);
              pronun_prev->rpLinks()[0].SetAcousticLike(0.0);
              pronun_prev->mpNext          = tnode;
            }
            pronun_prev = tnode;
          }
          if ((pronun_prev->rpLinks() =
                (LinkType *) malloc(sizeof(LinkType) * node->NLinks()))==NULL ||
            (pronun_first->rpBackLinks() =
                (LinkType *) malloc(sizeof(LinkType) * node->rNBackLinks())) == NULL) {
            Error("Insufficient memory");
          }
          pronun_prev->rNLinks()      = node->NLinks();
          pronun_first->rNBackLinks() = node->rNBackLinks();
    
          for (j = 0; j < node->rNBackLinks(); j++) 
          {
            NodeType* backnode(node->rpBackLinks()[j].pNode());

            backnode->rpLinks()[backnode->NLinks()].SetNode(pronun_first);
            backnode->rpLinks()[backnode->NLinks()].SetLmLike(node->rpBackLinks()[j].LmLike());
            backnode->rpLinks()[backnode->NLinks()].SetAcousticLike(node->rpBackLinks()[j].AcousticLike());
            backnode->rNLinks()++;
            pronun_first->rpBackLinks()[j] = node->rpBackLinks()[j];
          }

          for (j=0; j < node->NLinks(); j++) 
          {
            NodeType *forwnode = node->rpLinks()[j].pNode();
            forwnode->rpBackLinks()[forwnode->rNBackLinks()].SetNode(pronun_prev);
            forwnode->rpBackLinks()[forwnode->rNBackLinks()].SetLmLike(node->rpLinks()[j].LmLike());
            forwnode->rpBackLinks()[forwnode->rNBackLinks()].SetAcousticLike(node->rpLinks()[j].AcousticLike());
            forwnode->rNBackLinks()++;
            pronun_prev->rpLinks()[j] = node->rpLinks()[j];
          }
          if (prev != NULL) prev->mpNext = pronun_first;
          prev = pronun_prev;
        }
        prev->mpNext = node->mpNext;
        delete node->mC.mpAlphaBeta;
        node->mC.mpAlphaBeta = NULL;
        free(node->rpLinks());
        free(node->rpBackLinks());
        free(node);
        node = prev;
      }
    }
  // ExpandWordNetworkByDictionary( ...
  //****************************************************************************


  //***************************************************************************
  //***************************************************************************
    void 
    DecoderNetwork::
    ExpandMonophonesToTriphones(MyHSearchData *nonCDphones, 
        MyHSearchData *CDphones)
    {
      NodeType *  p_node;
      int     did_we_clone;
      int     i;
      int     j;
      int     k;

      // Find all Tee model, Word, and Null nodes (except the pFirstNode and last Null node)
      // and clone those not having single input and output
      do 
      {
        ENTRY           e    = {0}; //{0} is just to make compiler happy
        ENTRY*          ep;
        NodeType*      prev = NULL;
        
        did_we_clone = 0;
        
        for (p_node = pFirst(); p_node != NULL; prev = p_node, p_node = p_node->mpNext) 
        {
          if ((p_node->NLinks() == 0 || p_node->rNBackLinks() == 0) ||
              (p_node->NLinks() == 1 && p_node->rNBackLinks() == 1)) 
          {
            continue;
          }
          
          if (p_node->mC.mType & NT_PHONE) 
          {
            e.key = p_node->mC.mpName;
            my_hsearch_r(e, FIND, &ep, nonCDphones);
            if (ep == NULL || !reinterpret_cast<size_t>(ep->data))
              continue; // NodeType is not a Tee model
          }
          
          did_we_clone = 1;
          assert(prev != NULL); //Otherwise pFirstNode node is not Null node
    
          // Remove links to current node form back-linked nodes and realloc
          // link arrays of back-linked nodes to hold node->NLinks() more links
          for (j=0; j < p_node->rNBackLinks(); j++) 
          {
            NodeType *backnode = p_node->rpBackLinks()[j].pNode();
            
            for (k=0; k<backnode->NLinks() && backnode->rpLinks()[k].pNode()!=p_node; k++)
            {}
            
            assert(k < backnode->NLinks());
            
            // Otherwise link to 'node' is missing from which backlink exists
            backnode->rpLinks()[k] = backnode->rpLinks()[backnode->NLinks()-1];
    
            backnode->rpLinks() = 
              (LinkType *) realloc((backnode->rpLinks()), 
                               (backnode->NLinks()-1+p_node->NLinks())*sizeof(LinkType));
            
            if (backnode->rpLinks() == NULL) 
              Error("Insufficient memory");
            
            backnode->rNLinks()--;
          }
          
          // Remove backlinks to current node form linked nodes and realloc
          // backlink arrays of linked nodes to hold node->rNBackLinks() more backlinks
          for (j=0; j < p_node->NLinks(); j++) 
          {
            NodeType *forwnode = p_node->rpLinks()[j].pNode();
            
            for (k=0;k<forwnode->rNBackLinks()&&forwnode->rpBackLinks()[k].pNode()!=p_node;k++)
            {}
            
            assert(k < forwnode->rNBackLinks());
            // Otherwise link to 'p_node' is missing from which backlink exists
            forwnode->rpBackLinks()[k] = forwnode->rpBackLinks()[forwnode->rNBackLinks()-1];
    
            forwnode->rpBackLinks() = 
              (LinkType *) realloc((forwnode->rpBackLinks()),
                               (forwnode->rNBackLinks()-1+p_node->rNBackLinks())*sizeof(LinkType));
                               
            if (forwnode->rpBackLinks() == NULL) 
              Error("Insufficient memory");
              
            forwnode->rNBackLinks()--;
          }
          // Alloc new p_node->rNLinks() * p_node->rNBackLinks() nodes and create new links
          // so that each backlinked node is conected with each linked node through
          // one new node.
          for (i=0; i < p_node->NLinks(); i++) 
          {
            for (j=0; j < p_node->rNBackLinks(); j++) 
            {
              NodeType *  tnode;
              LinkType    forwlink = p_node->rpLinks()[i];
              LinkType    backlink = p_node->rpBackLinks()[j];
    
              if ((tnode = (NodeType *) calloc(1, sizeof(NodeType))) == NULL)
                Error("Insufficient memory");
              
              tnode->mC.mpName = NULL;
              *tnode = *p_node;
    
              if ((tnode->rpLinks()   = (LinkType *) malloc(sizeof(LinkType))) == NULL ||
                (tnode->rpBackLinks() = (LinkType *) malloc(sizeof(LinkType))) == NULL) 
              {
                Error("Insufficient memory");
              }
              
              tnode->rNLinks()       = 1;
              tnode->rNBackLinks()   = 1;
              tnode->rpLinks()[0]     = forwlink;
              tnode->rpBackLinks()[0] = backlink;
              
              forwlink.pNode()->rpBackLinks()[forwlink.pNode()->rNBackLinks()].SetNode(tnode);
              forwlink.pNode()->rpBackLinks()[forwlink.pNode()->rNBackLinks()].SetLmLike(forwlink.LmLike());
              forwlink.pNode()->rpBackLinks()[forwlink.pNode()->rNBackLinks()].SetAcousticLike(forwlink.AcousticLike());

              forwlink.pNode()->rNBackLinks()++;

              backlink.pNode()->rpLinks()    [backlink.pNode()->NLinks()   ].SetNode(tnode);
              backlink.pNode()->rpLinks()    [backlink.pNode()->NLinks()   ].SetLmLike(backlink.LmLike());
              backlink.pNode()->rpLinks()    [backlink.pNode()->NLinks()   ].SetAcousticLike(backlink.AcousticLike());

              backlink.pNode()->rNLinks()++;

              prev->mpNext = tnode;
              prev = tnode;
            }
          }
          prev->mpNext = p_node->mpNext;
          delete p_node->mC.mpAlphaBeta;
          p_node->mC.mpAlphaBeta = NULL;
          free(p_node->rpLinks());
          free(p_node->rpBackLinks());
          free(p_node);
          p_node = prev;
        }
      } while (did_we_clone);
    
      // Assign to each node unique number, which will later allow to find groups of
      // expanded triphone nodes corresponding to original monophone nodes.
      int nbackmononodes;
      int nforwmononodes;
      int id = 0;
      
      for (p_node = pFirst(); p_node != NULL; p_node = p_node->mpNext) 
        p_node->mAux = id++;
    
      // Expand monophone nodes to triphone nodes
      NodeType *prev = NULL;
      for (p_node = pFirst(); p_node != NULL; prev = p_node, p_node = p_node->mpNext) 
      {
        ENTRY e = {0}; //{0} is just to make compiler happy
        ENTRY * ep;
    
        if ((p_node->mC.mType & NT_WORD) ||
            (p_node->NLinks() == 1 && p_node->rNBackLinks() == 1)) 
        {
          continue;
        }
        
        assert(p_node->mC.mType & NT_PHONE);
        e.key = p_node->mC.mpName;
        my_hsearch_r(e, FIND, &ep, nonCDphones);
        if (ep != NULL && reinterpret_cast<size_t>(ep->data)) continue; // NodeType is a Tee model
    
        assert(prev != NULL); //Otherwise first node is not Null node
    
        // Count groups of backlinked nodes corresponding to different monophones
        id = -1;
        nbackmononodes = 0;
        for (j=0; j < p_node->rNBackLinks(); j++) 
        {
          if (p_node->rpBackLinks()[j].pNode()->mAux != id) 
          {
            id = p_node->rpBackLinks()[j].pNode()->mAux;
            nbackmononodes++;
          }
        }
        // Count groups of linked nodes corresponding to different monophones
        id = -1;
        nforwmononodes = 0;
        for (j=0; j < p_node->NLinks(); j++) 
        {
          if (p_node->rpLinks()[j].pNode()->mAux != id) 
          {
            id = p_node->rpLinks()[j].pNode()->mAux;
            nforwmononodes++;
          }
        }
    
        // Remove links to current node form backlinked nodes and realloc
        // link arrays of backlinked nodes to hold nforwmononodes more links
        for (j=0; j < p_node->rNBackLinks(); j++) 
        {
          NodeType *backnode = p_node->rpBackLinks()[j].pNode();
          for (k=0; k<backnode->NLinks() && backnode->rpLinks()[k].pNode()!=p_node; k++);
          assert(k < backnode->NLinks());
          // Otherwise link to 'node' is missing from which backlink exists
          memmove(backnode->rpLinks()+k, backnode->rpLinks()+k+1,
                  (backnode->NLinks()-k-1) * sizeof(LinkType));
    
          backnode->rpLinks() = (LinkType *)
            realloc(backnode->rpLinks(),
                  (backnode->NLinks()-1+nforwmononodes)*sizeof(LinkType));
          if (backnode->rpLinks() == NULL) Error("Insufficient memory");
          backnode->rNLinks()--;
        }
    
        // Remove backlinks to current p_node form linked nodes and realloc
        // backlink arrays of linked nodes to hold nbackmononodes more backlinks
        for (j=0; j < p_node->NLinks(); j++) 
        {
          NodeType *forwnode = p_node->rpLinks()[j].pNode();
          for (k=0;k<forwnode->rNBackLinks()&&forwnode->rpBackLinks()[k].pNode()!=p_node;k++);
          assert(k < forwnode->rNBackLinks());
          // Otherwise link to 'p_node' is missing from which backlink exists
          memmove(forwnode->rpBackLinks()+k, forwnode->rpBackLinks()+k+1,
                  (forwnode->rNBackLinks()-k-1) * sizeof(LinkType));
    
          forwnode->rpBackLinks() = (LinkType *)
            realloc(forwnode->rpBackLinks(),
                  (forwnode->rNBackLinks()-1+nbackmononodes)*sizeof(LinkType));
          if (forwnode->rpBackLinks() == NULL) Error("Insufficient memory");
          forwnode->rNBackLinks()--;
        }
    
        // Alloc new nforwmononodes * nbackmononodes nodes and create new links
        // so that each backlinked node is conected through one new node with all
        // linked nodes belonging to one monophone group and vice versa each
        // linked node is conected through one new node with all backlinked nodes
        // belonging to one monophone group
        LinkType *forwmono_start, *forwmono_end = p_node->rpLinks();
        for (i=0; i < nforwmononodes; i++) 
        {
          for (forwmono_start = forwmono_end;
              forwmono_end < p_node->rpLinks()+p_node->NLinks() &&
              forwmono_start->pNode()->mAux == forwmono_end->pNode()->mAux;
              forwmono_end++)
          {}
    
          assert((i <  nforwmononodes-1 && forwmono_end <  p_node->rpLinks()+p_node->NLinks()) ||
                (i == nforwmononodes-1 && forwmono_end == p_node->rpLinks()+p_node->NLinks()));
    
          LinkType *tlink, *backmono_start, *backmono_end = p_node->rpBackLinks();
          
          for (j=0; j < nbackmononodes; j++) 
          {
            for (backmono_start = backmono_end;
              backmono_end < p_node->rpBackLinks()+p_node->rNBackLinks() &&
              backmono_start->pNode()->mAux == backmono_end->pNode()->mAux;
              backmono_end++)
            {}
    
            assert((j <  nbackmononodes-1 && backmono_end <  p_node->rpBackLinks()+p_node->rNBackLinks()) ||
                   (j == nbackmononodes-1 && backmono_end == p_node->rpBackLinks()+p_node->rNBackLinks()));
    
            NodeType * tnode;
            if ((tnode = (NodeType *) calloc(1, sizeof(NodeType))) == NULL)
              Error("Insufficient memory");
            
            *tnode = *p_node;
            tnode->rNLinks()       = forwmono_end-forwmono_start;
            tnode->rNBackLinks()   = backmono_end-backmono_start;
    
            if ((tnode->rpLinks() =
                (LinkType *) malloc(tnode->NLinks() * sizeof(LinkType))) == NULL ||
              (tnode->rpBackLinks() =
                (LinkType *) malloc(tnode->rNBackLinks() * sizeof(LinkType))) == NULL) 
            {
              Error("Insufficient memory");
            }
            
            for (tlink = forwmono_start; tlink < forwmono_end; tlink++) 
            {
              tnode->rpLinks()[tlink-forwmono_start] = *tlink;
              tlink->pNode()->rpBackLinks()[tlink->pNode()->rNBackLinks()].SetNode(tnode);
              tlink->pNode()->rpBackLinks()[tlink->pNode()->rNBackLinks()].SetLmLike(tlink->LmLike());
              tlink->pNode()->rpBackLinks()[tlink->pNode()->rNBackLinks()].SetAcousticLike(tlink->AcousticLike());

              tlink->pNode()->rNBackLinks()++;
            }
            
            for (tlink = backmono_start; tlink < backmono_end; tlink++) 
            {
              tnode->rpBackLinks()[tlink-backmono_start] = *tlink;
              tlink->pNode()->rpLinks()[tlink->pNode()->NLinks()].SetNode(tnode);
              tlink->pNode()->rpLinks()[tlink->pNode()->NLinks()].SetLmLike(tlink->LmLike());
              tlink->pNode()->rpLinks()[tlink->pNode()->NLinks()].SetAcousticLike(tlink->AcousticLike());

              tlink->pNode()->rNLinks()++;
            }
            
            prev->mpNext = tnode;
            prev = tnode;
          }
        }
        prev->mpNext = p_node->mpNext;
        delete p_node->mC.mpAlphaBeta;
        p_node->mC.mpAlphaBeta = NULL;
        free(p_node->rpLinks());
        free(p_node->rpBackLinks());
        free(p_node);
        p_node = prev;
      }

      // Give triphone names to phone nodes and create hash of these names
      for (p_node = pFirst(); p_node != NULL; p_node = p_node->mpNext) 
      {
        ENTRY     e       = {0}; //{0} is just to make compiler happy
        ENTRY*    ep      = NULL;
        NodeType*     lc      = NULL;
        NodeType*     rc      = NULL;
        char*     lcname  = NULL;
        char*     rcname  = NULL;
        char*     triname = NULL;
        int       lcnlen  = 0;
        int       rcnlen  = 0;
    
        if (!(p_node->mC.mType & NT_PHONE)) 
          continue;
    
        if (nonCDphones) 
        {
          e.key  = p_node->mC.mpName;
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
          for (lc = p_node;;) 
          {
            lc = lc->rNBackLinks() ? lc->rpBackLinks()[0].pNode() : NULL;
            
            if (lc == NULL)               break;
            if (!(lc->mC.mType & NT_PHONE))  continue;
            if (nonCDphones == NULL)      break;
            
            e.key  = lc->mC.mpName;
            my_hsearch_r(e, FIND, &ep, nonCDphones);
            if (ep == NULL || !reinterpret_cast<size_t>(ep->data)) break; // NodeType represents Tee model
          }
          
          for (rc = p_node;;) 
          {
            rc = rc->NLinks() ? rc->rpLinks()[0].pNode() : NULL;
            
            if (rc == NULL)               break;
            if (!(rc->mC.mType & NT_PHONE))  continue;
            if (nonCDphones == NULL)      break;
            
            e.key  = rc->mC.mpName;
            my_hsearch_r(e, FIND, &ep, nonCDphones);
            if (ep == NULL || !reinterpret_cast<size_t>(ep->data)) break; // NodeType represents Tee model
          }
        }
        
        lcnlen = -1;
        if (lc != NULL) 
        {
          lcname = strrchr(lc->mC.mpName, '-');
          
          if (lcname == NULL) 
            lcname = lc->mC.mpName;
          else 
            lcname++;
            
          lcnlen = strcspn(lcname, "+");
        }
        
        rcnlen = -1;
        if (rc != NULL) 
        {
          rcname = strrchr(rc->mC.mpName, '-');
          
          if (rcname == NULL) 
            rcname = rc->mC.mpName;
          else 
            rcname++;
            
          rcnlen = strcspn(rcname, "+");
        }
        
        triname = (char *) malloc(lcnlen+1+strlen(p_node->mC.mpName)+1+rcnlen+1);
        
        if (triname == NULL) 
          Error("Insufficient memory");
    
        triname[0] = '\0';
    
        if (lcnlen > 0) 
          strcat(strncat(triname, lcname, lcnlen), "-");
          
        strcat(triname, p_node->mC.mpName);
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
          
          p_node->mC.mpName = triname;
        } 
        else 
        {
          free(triname);
          p_node->mC.mpName = ep->key;
        }
      }
    }
  // ExpandMonophoneNetworkToTriphones(MyHSearchData *nonCDphones, ...
  //****************************************************************************

  
  //****************************************************************************
  //****************************************************************************
    void 
    DecoderNetwork::
    ExpansionsAndOptimizations(
      ExpansionOptions        expOptions,
      const STKNetworkOutputFormat&  rFormat,
      MyHSearchData *         wordHash,
      MyHSearchData *         nonCDphHash,
      MyHSearchData *         triphHash,
      FLOAT                   wordPenalty,
      FLOAT                   modelPenalty,
      FLOAT                   lmScale,
      FLOAT                   posteriorScale)
    {
      //NodeType*              p_node(pFirst());

      if (expOptions.mNoWordExpansion  && !expOptions.mCDPhoneExpansion &&
          expOptions.mNoOptimization   && !rFormat.mNoLMLikes &&
          !rFormat.mNoTimes        && !rFormat.mNoWordNodes &&
          !rFormat.mNoModelNodes   && !rFormat.mNoPronunVars && 
          !expOptions.mRemoveNulls) 
      {
        return;
      }
      
      SelfLinksToNullNodes();

      if (!expOptions.mNoWordExpansion) {
        if (!expOptions.mNoOptimization) {
          LatticeLocalOptimization(expOptions, wordPenalty, modelPenalty,
              lmScale, posteriorScale);
        }
        assert(wordHash != NULL);
        ExpandByDictionary(wordHash, !expOptions.mRemoveWordsNodes, 
            !expOptions.mRespectPronunVar);
      }

      if (expOptions.mCDPhoneExpansion) {
        if (!expOptions.mNoOptimization) {
          LatticeLocalOptimization(expOptions, wordPenalty, modelPenalty,
              lmScale, posteriorScale);
        }
        assert(triphHash != NULL && nonCDphHash != NULL);
        ExpandMonophonesToTriphones(nonCDphHash, triphHash);
      }

      DiscardUnwantedInfo(rFormat);    
      RemoveRedundantNullNodes(expOptions.mRemoveNulls == 1, wordPenalty, 
          modelPenalty, lmScale, posteriorScale);

      if (!expOptions.mNoOptimization) {
        LatticeLocalOptimization(expOptions, wordPenalty, modelPenalty,
            lmScale, posteriorScale);
        RemoveRedundantNullNodes(false, wordPenalty, modelPenalty, lmScale, 
            posteriorScale);
      }
    } 
  // void NetworkExpansionsAndOptimizations( )
  //****************************************************************************

  
  //****************************************************************************
  //****************************************************************************
    void 
    DecoderNetwork::
    DiscardUnwantedInfo(const STKNetworkOutputFormat& format)
    {
      // TODO : change to iterator
      iterator    p_node;
      int         i;
    
      for (p_node = begin(); p_node != end(); p_node++)  
      {
        if (format.mNoLMLikes && !format.mPosteriors) 
        {
          for (i=0; i < p_node->NLinks();     i++) 
            p_node->rpLinks()    [i].SetLmLike(0.0);
          
          for (i=0; i < p_node->rNBackLinks(); i++) 
            p_node->rpBackLinks()[i].SetLmLike(0.0);
        }

        if (format.mNoAcousticLikes && !format.mPosteriors) 
        {
          for (i=0; i < p_node->NLinks();     i++) 
            p_node->rpLinks()    [i].SetAcousticLike(0.0);

          for (i=0; i < p_node->rNBackLinks(); i++) 
            p_node->rpBackLinks()[i].SetAcousticLike(0.0);
        }

        if (format.mNoTimes) 
        {
          p_node->mC.SetStart(UNDEF_TIME);
          p_node->mC.SetStop(UNDEF_TIME);
        }
        if (format.mNoWordNodes && p_node->mC.mType & NT_WORD) {
          p_node->mC.mpPronun = NULL;
        }
        if (format.mNoModelNodes && (p_node->mC.mType&NT_MODEL || p_node->mC.mType&NT_PHONE)) {
          p_node->mC.mType = NT_WORD;
          p_node->mC.mpPronun = NULL;
        }
        if (format.mNoPronunVars && p_node->mC.mType & NT_WORD && p_node->mC.mpPronun != NULL) {
          p_node->mC.mpPronun = p_node->mC.mpPronun->mpWord->pronuns[0];
        }
      }
    }
    // DiscardUnwantedInfo(STKNetworkOutputFormat format)
    //***************************************************************************
    


  //***************************************************************************
  //***************************************************************************
    void
    DecoderNetwork:: 
    SortNodes()
    {
      int       i;
      int       j;
      NodeType *    chain;
      NodeType *    last;
      NodeType *    node;
    
      // Sort nodes for forward (Viterbi) propagation
    
      for (node = pFirst(); node != NULL; node = node->mpNext) 
      {
        node->mAux = node->rNBackLinks();
      }
    
      for (i = 0; i < pFirst()->rNLinks(); i++) 
      {
        pFirst()->rpLinks()[i].pNode()->mAux--;
      }
    
      last = pFirst();
      chain = pFirst()->mpNext;
    
      while (chain) 
      {
        bool    short_curcuit = true;
        NodeType ** curPtr = &chain;
        i = 0;
    
        while (*curPtr) 
        {
          if (((((*curPtr)->Content().mType & NT_MODEL) || ((*curPtr)->Content().mType & NT_PHONE)) && !((*curPtr)->Content().mType & NT_TEE))
            || (*curPtr)->mAux == 0) 
          {
            for (j = 0; j < (*curPtr)->rNLinks(); j++) 
            {
              (*curPtr)->rpLinks()[j].pNode()->mAux--;
            }
    
            last = (last->mpNext = *curPtr);
            last->mAux = i++;
            *curPtr = (*curPtr)->mpNext;
            short_curcuit = false;
          } 
          else 
          {
            curPtr = &(*curPtr)->mpNext;
          }
        }
    
        if (short_curcuit) 
        {
    //      fprintf(stderr, "Nodes in loop: ");
    //      for (curPtr = &chain; *curPtr; curPtr = &(*curPtr)->next)
    //        fprintf(stderr, "%d %d", *curPtr - mpNodes, (*curPtr)->mType);
    //      fprintf(stderr, "\n");
          Error("Loop of non-emiting nodes found in network");
        }
      }
    
      last->mpNext = NULL;
    
      // /// !!! What is this sorting links good for ???
      // for (node = pFirst(); node != NULL; node = node->mpNext) 
      // {
      //   if (node->rNLinks() > 1)
      //     qsort(node->rpLinks(), node->rNLinks(), sizeof(LinkType), cmplnk);
      // }
    
    // Sort nodes for backward propagation
    
      for (node = pFirst(); node != NULL; node = node->mpNext)
        node->mAux = node->rNLinks();
    
      for (i = 0; i < pLast()->rNBackLinks(); i++)
        pLast()->rpBackLinks()[i].pNode()->mAux--;
    
      last = pLast();
      chain = pLast()->mpBackNext;
      i = 0;
    
      while (chain) 
      {
        bool short_curcuit = true;
        NodeType **curPtr = &chain;
    
        while (*curPtr) 
        {
          if (((((*curPtr)->Content().mType & NT_MODEL) || ((*curPtr)->Content().mType & NT_PHONE)) && !((*curPtr)->Content().mType & NT_TEE))
            || (*curPtr)->mAux == 0) 
          {
            for (j = 0; j < (*curPtr)->rNBackLinks(); j++) 
            {
              (*curPtr)->rpBackLinks()[j].pNode()->mAux--;
            }
    
            last = (last->mpBackNext = *curPtr);
            last->mAux = i++;
            *curPtr = (*curPtr)->mpBackNext;
            short_curcuit = false;
          } 
          else 
          {
            curPtr = &(*curPtr)->mpBackNext;
          }
        }
    
        assert(!short_curcuit); // Shouldn't happen, since it didnot happen before
      }
    
      last->mpBackNext = NULL;
    
      // /// !!! What is this sorting links good for ???
      // for (node = pFirst(); node != NULL; node = node->mpNext) 
      // {
      //   if (node->rNBackLinks() > 1)
      //     qsort(node->rpBackLinks(), node->rNBackLinks(), sizeof(LinkType), cmplnk);
      // }
    } // Decoder::SortNodes();

  //***************************************************************************
  //***************************************************************************
    void 
    DecoderNetwork::
    TopologicalSort()
    {
      NodeType *    node;
      NodeType *    lastnode;
      int       i;
      int       j;
      int       unreachable = 0;

      // Sort nodes in topological order
      // printf("Sorting nodes...\n");
      pFirst()->mAux = 1;
      for (lastnode = node = pFirst(); node != NULL; node = node->mpNext) 
      {
        for (i=0; i < node->NLinks(); i++) 
        {
          NodeType *lnknode = node->rpLinks()[i].pNode();
          
          if (lnknode->mAux == 0) 
          {
            for (j=0; j<lnknode->rNBackLinks() && lnknode->rpBackLinks()[j].pNode()->mAux==1; j++)
            {}
            
            if (j == lnknode->rNBackLinks()) 
            {
              lastnode->mpNext = lnknode;
              lastnode  = lnknode;
              lnknode->mAux = 1;
              lnknode->mpNext = NULL;
            }
          }
        }
      }
      
      if (lastnode->NLinks() != 0) 
      {
        // There is a cycle in graph so we cannot sort nodes
        // topologicaly, so sort it at least somehow. :o|
        // Anyway this optimization algorithm is not optimal for graphs with cycles.
        for (node = pFirst(); node != NULL; node = node->mpBackNext) 
        {
          node->mAux = 0;
        }
          
        pFirst()->mAux = 1;
        for (lastnode = node = pFirst(); node != NULL; node = node->mpNext) 
        {
          for (i=0; i < node->NLinks(); i++) 
          {
            NodeType *lnknode = node->rpLinks()[i].pNode();
            if (lnknode->mAux == 0) 
            {
              lastnode->mpNext = lnknode;
              lastnode  = lnknode;
              lnknode->mAux = 1;
              lnknode->mpNext = NULL;
            }
          }
        }
        
        for (node=pFirst(); node->mpNext->NLinks() != 0; node=node->mpNext)
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
    
      for (node=pFirst(); node != NULL; node=node->mpBackNext) 
      {
        while (node->mpBackNext && node->mpBackNext->mAux == 0) 
        {
          NodeType *tnode = node->mpBackNext;
          node->mpBackNext = node->mpBackNext->mpBackNext;
          unreachable++;
          delete tnode->mC.mpAlphaBeta;
          tnode->mC.mpAlphaBeta = NULL;
          free(tnode->rpLinks());
          free(tnode->rpBackLinks());
          free(tnode);
        }
      }
    
      //  if (unreachable) Warning("Removing %d unreachable nodes", unreachable);
      if (unreachable) 
        Error("Networks contains unreachable nodes");
    
      pFirst()->mpBackNext = NULL;
      
      for (node=pFirst(); node->mpNext != NULL; node=node->mpNext) 
      {
        node->mpNext->mpBackNext = node;
      }
    }

  //***************************************************************************
  //***************************************************************************
    void 
    DecoderNetwork::
    LatticeLocalOptimization(const ExpansionOptions &expOptions,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale)
    {
      NodeType *    node;
      int       i;
      
      // For each node, sort links by pointer value to allow
      // for easy comparison whether two nodes have the same set of links
      for (node = pFirst(); node != NULL; node = node->mpNext)  
      {
        node->mAux = 0;
        node->mpBackNext = node->mpNext;
        qsort(node->rpLinks(), node->NLinks(), sizeof(LinkType), lnkcmp);
        qsort(node->rpBackLinks(), node->rNBackLinks(), sizeof(LinkType), lnkcmp);
      }
    
      TopologicalSort();
      
      for (i=1, node=pFirst(); node != NULL; node = node->mpNext, i++) 
      {
        node->mAux=i;
      }
      
      for (;;) 
      {
        if (expOptions.mTraceFlag & 2) 
        {
          for (i=0,node=pFirst(); node; node=node->mpNext,i++)
          {}
          
          TraceLog("Forward pass.... (number of nodes: %d)", i);
        }
        
        LatticeLocalOptimization_ForwardPass(expOptions, wordPenalty, 
            modelPenalty, lmScale, posteriorScale);
    
        if (expOptions.mTraceFlag & 2) 
        {
          for (i=0,node=pFirst(); node; node=node->mpNext,i++)
          {}
          
          TraceLog("Backward pass... (number of nodes: %d)", i);
        }
        
        if (!LatticeLocalOptimization_BackwardPass(expOptions, wordPenalty, 
              modelPenalty, lmScale, posteriorScale)) {
          break;
        }
      }
    }  
  // LatticeLocalOptimization(const ExpansionOptions &expOptions)
  //****************************************************************************


  //****************************************************************************
  //****************************************************************************
    int 
    DecoderNetwork::
    LatticeLocalOptimization_ForwardPass(const ExpansionOptions &expOptions,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale)
    {
      int     i; 
      int     j;
      int     k; 
      int     l; 
      int     m;
      int     rep;
      NodeType*   p_tnode;
      int     node_removed = 0;
      //FLOAT   t_acoustic_like;
      FLOAT   t_lm_like;

      
      //for (NodeType* p_node = pFirst(); p_node != NULL; p_node = p_node->mpNext) 
      for (iterator p_node = begin(); p_node != end(); ++p_node) 
      {
        if(!expOptions.mNoWeightPushing) {
          for (i = 0; i < p_node->NLinks(); i++) {
            p_tnode = p_node->rpLinks()[i].pNode();

            if (p_tnode->NLinks() == 0) 
              continue;

            // Weight pushing
            t_lm_like = p_tnode->rpBackLinks()[0].LmLike();
            //t_acoustic_like = p_tnode->rpBackLinks()[0].AcousticLike();

            for (l=1; l <  p_tnode->rNBackLinks(); l++) 
            {
              t_lm_like = HIGHER_OF(t_lm_like, p_tnode->rpBackLinks()[l].LmLike());
              //t_acoustic_like = HIGHER_OF(t_acoustic_like, p_tnode->rpBackLinks()[l].AcousticLike());
            }

            for (l=0; l < p_tnode->rNBackLinks(); l++) 
            {
              NodeType* backnode = p_tnode->rpBackLinks()[l].pNode();

              p_tnode->rpBackLinks()[l].AddLmLike(-t_lm_like);
              //p_tnode->rpBackLinks()[l].AddAcousticLike(-t_acoustic_like);

              for (k=0; k<backnode->NLinks() && backnode->rpLinks()[k].pNode()!=p_tnode; k++)
              {}

              assert(k < backnode->NLinks());

              backnode->rpLinks()[k].AddLmLike(-t_lm_like);
              //backnode->rpLinks()[k].AddAcousticLike(-t_acoustic_like);

#ifndef NDEBUG
              for (k++; k<backnode->NLinks() && backnode->rpLinks()[k].pNode()!=p_tnode; k++)
              {}
#endif
              assert(k == backnode->NLinks());
            }
            
            for (l=0; l < p_tnode->NLinks(); l++) 
            {
              NodeType* forwnode = p_tnode->rpLinks()[l].pNode();
  
              p_tnode->rpLinks()[l].AddLmLike(t_lm_like);
              //p_tnode->rpLinks()[l].AddAcousticLike(t_acoustic_like);
              
              for (k=0; k<forwnode->rNBackLinks() && forwnode->rpBackLinks()[k].pNode()!=p_tnode;k++)
              {}
              
              assert(k < forwnode->rNBackLinks());
  
              forwnode->rpBackLinks()[k].AddLmLike(t_lm_like);
              //forwnode->rpBackLinks()[k].AddAcousticLike(t_acoustic_like);
  
#ifndef NDEBUG
              for (k++; k<forwnode->rNBackLinks() && forwnode->rpBackLinks()[k].pNode()!=p_tnode;k++)
              {}
#endif
              assert(k == forwnode->rNBackLinks());
            }
          }
        }      
        //dnet(pFirstNode, 1, p_node);

        // For current node 'p_node', check for each possible pair of its successors
        // ('inode' and 'jnode') whether the pair may be merged to single node.
        for (i = 0; i < p_node->NLinks()-1; i++) 
        {
          for (j = i+1; j < p_node->NLinks(); j++) 
          {
            //NodeType* inode = p_node->rpLinks()[i].pNode();
          //NodeType* jnode = p_node->rpLinks()[j].pNode();
            iterator inode(p_node->rpLinks()[i].pNode());
            iterator jnode(p_node->rpLinks()[j].pNode());

            // Final node may be never merged.
            if (inode->NLinks() == 0 || jnode->NLinks() == 0) {
              continue;
            }

            // Two nodes ('inode' and 'jnode') may be mergeg if they are of the 
            // same type, name, ... with the same predecessors and with the same
            // weights on the links from predecesors.
            if ((inode->mC.mType & ~NT_TRUE) != (jnode->mC.mType & ~NT_TRUE)
            || ( inode->mC.mType & NT_PHONE && inode->mC.mpName   != jnode->mC.mpName)
            || ( inode->mC.mType & NT_WORD  && inode->mC.mpPronun != jnode->mC.mpPronun)
    //          &&  (inode->mpPronun == NULL ||
    //             jnode->mpPronun == NULL ||
    //             inode->mpPronun->mpWord       != jnode->mpPronun->mpWord ||
    //             inode->mpPronun->outSymbol  != jnode->mpPronun->outSymbol ||
    //             inode->mpPronun->variant_no != jnode->mpPronun->variant_no ||
    //             inode->mpPronun->prob       != jnode->mpPronun->prob)
            || (inode->NBackLinks() != jnode->NBackLinks())) 
            {
              continue;
            }
            
            if (expOptions.mStrictTiming && (inode->mC.Start() != jnode->mC.Start()
                            ||  inode->mC.Stop()  != jnode->mC.Stop())) 
            {
              continue;
            }

            // Weights on the links from predecesors does not have to be exactely
            // the same, but the must not differ more than by 
            // SIGNIFICANT_PROB_DIFFERENCE
            for (l=0; l < inode->NBackLinks(); l++) 
            {
              if (inode->rpBackLinks()[l].pNode() != jnode->rpBackLinks()[l].pNode()) 
                break;
              
              FLOAT ldiff =  inode->rpBackLinks()[l].LmLike() 
                - jnode->rpBackLinks()[l].LmLike();
              FLOAT adiff =  inode->rpBackLinks()[l].AcousticLike() 
                - jnode->rpBackLinks()[l].AcousticLike();

              // TODO: is the thing with acoustic like OK???
              if (ldiff < -SIGNIFICANT_PROB_DIFFERENCE 
              ||  ldiff >  SIGNIFICANT_PROB_DIFFERENCE
              ||  adiff < -SIGNIFICANT_PROB_DIFFERENCE 
              ||  adiff >  SIGNIFICANT_PROB_DIFFERENCE ) 
              {
                break;
              }
            }
            
            if (l < inode->NBackLinks()) 
              continue;
    
    /*        if (memcmp(inode->rpBackLinks(), jnode->rpBackLinks(),
                      inode->NBackLinks() * sizeof(Link<NODE_REGULAR, LINK_BASIC>))) {
              continue;
            }*/
              // inode and jnode are the same nodes with the same predeccessors
              // Remove jnode and add its links to inode
    
            assert(inode->NLinks() && jnode->NLinks());
    
            //TraceLog("Removing node: %s", 
            //    inode->mType & NT_PHONE ? inode->mpName : 
            //    inode->mType & NT_WORD  ? (inode->mpPronun ? 
            //         inode->mpPronun->mpWord->mpName : "!NULL") : "UNKNOWN_TYPE");


            // Remove links to jnode form predeccessors
            for (l=0; l < jnode->NBackLinks(); l++) 
            {
              NodeType* backnode = jnode->rpBackLinks()[l].pNode();

              for (k=0; k<backnode->NLinks() && backnode->rpLinks()[k].pNode()!=jnode.mpPtr;
                   k++)
              { }

              assert(k < backnode->NLinks());
              // Otherwise link to 'p_node' is missing from which backlink exists
              memmove(backnode->rpLinks()+k, backnode->rpLinks()+k+1,
                      (backnode->NLinks()-k-1) * sizeof(LinkType));

              backnode->rNLinks()--;
            }
    
            // Merge jnode's links with inode links
    
            //Count jnode's links not present among inode's links
            rep = l = k = 0;
            while (k < jnode->NLinks()) 
            {
              LinkType* ill = inode->rpLinks()+l;
              LinkType* jlk = jnode->rpLinks()+k;

              if (l == inode->NLinks() || ill->pNode() > jlk->pNode())
              {
                // k-th link of jnode will be included among inode's links.
                // Redirect corresponding baclink to inode
                for (m = 0; m < jlk->pNode()->NBackLinks()
                            && jlk->pNode()->rpBackLinks()[m].pNode() != jnode.mpPtr; m++)
                {}
                
                assert(m < jlk->pNode()->NBackLinks());
                jlk->pNode()->rpBackLinks()[m].SetNode(inode.mpPtr);
                qsort(jlk->pNode()->rpBackLinks(), jlk->pNode()->NBackLinks(),
                      sizeof(LinkType), lnkcmp);
                k++;
              } 
              else if (ill->pNode() == jlk->pNode()) 
              {
                // l-th link of inode and k-th link of jnode points to
                // the same node. Link from jnode is redundant.
                // Remove backlinks to jnode form jnode's succesors
                for (m = 0; m < jlk->pNode()->rNBackLinks()
                          && jlk->pNode()->rpBackLinks()[m].pNode() != jnode.mpPtr; m++);
                {}
                
                assert(m < jlk->pNode()->rNBackLinks());
                memmove(jlk->pNode()->rpBackLinks()+m, jlk->pNode()->rpBackLinks()+m+1,
                        (jlk->pNode()->rNBackLinks()-m-1) * sizeof(LinkType));

                jlk->pNode()->rNBackLinks()--;
    
                // update likelihood correctly
                if (ill->Like(lmScale) < jlk->Like(lmScale))
                {
                  ill->SetLmLike(jlk->LmLike());
                  ill->SetAcousticLike(jlk->AcousticLike());
                }

                jlk->SetNode(NULL); // Mark link to be removed

                rep++; 
                k++;
                l++;
              } 
              else 
              {
                l++;
              }
            }
            
            l = inode->NLinks();
            inode->rNLinks() += jnode->NLinks()-rep;
            inode->rpLinks() = (LinkType *) realloc(inode->rpLinks(),
                                            inode->NLinks() * sizeof(LinkType));
            
            if (inode->rpLinks() == NULL) 
              Error("Insufficient memory");
    
            for (k = 0; k < jnode->NLinks(); k++) 
            {
              if (jnode->rpLinks()[k].pNode() != NULL) 
                inode->rpLinks()[l++] = jnode->rpLinks()[k];
            }
            
            qsort(inode->rpLinks(), inode->NLinks(), sizeof(LinkType), lnkcmp);
    
            inode->mC.SetStart(inode->mC.Start() == UNDEF_TIME || jnode->mC.Start() == UNDEF_TIME
                            ? UNDEF_TIME 
                            : LOWER_OF(inode->mC.Start(), jnode->mC.Start()));
    
            inode->mC.SetStop (inode->mC.Stop() == UNDEF_TIME || jnode->mC.Stop() == UNDEF_TIME
                            ? UNDEF_TIME 
                            : HIGHER_OF(inode->mC.Stop(), jnode->mC.Stop()));
    
            if (inode->mAux > jnode->mAux) 
            {
            // Make sure that topological order of new inode's links
            // (inherited from jnode) is higher than inode's order.
            // In the 'next' list, move inode to jnode's lower position
              inode->mpBackNext->mpNext = inode->mpNext;
              inode->mpNext->mpBackNext = inode->mpBackNext;
              inode->mpNext             = jnode->mpNext;
              inode->mpBackNext         = jnode->mpBackNext;
              inode->mpBackNext->mpNext = inode.mpPtr;
              inode->mpNext->mpBackNext = inode.mpPtr;
              inode->mAux = jnode->mAux;
            } 
            else
            {
              jnode->mpNext->mpBackNext = jnode->mpBackNext;
              jnode->mpBackNext->mpNext = jnode->mpNext;
            }
            
            inode->mC.mType |= jnode->mC.mType & NT_TRUE;
            
            delete jnode->mC.mpAlphaBeta;
            jnode->mC.mpAlphaBeta = NULL;
            free(jnode->rpLinks());
            free(jnode->rpBackLinks());
            free(jnode.mpPtr);
//            jnode.mpPtr->mpNext = NULL;
//            jnode.mpPtr->mpBackNext = NULL;
//            erase(jnode);

            --j; // Process j-th node again
                 // there is new shifted node on this index
    
            node_removed = 1;
          }
        }
      }
      return node_removed;
    }
  //  LatticeLocalOptimization_ForwardPass(const ExpansionOptions &expOptions)
  //****************************************************************************

  //***************************************************************************
  //***************************************************************************
  struct CorrPhnRec 
  {
    DecoderNetwork::NodeType*   mpNode;
    TimingType              maxStopTimeTillNow;
    int                     mId;
  };

  //***************************************************************************
  //***************************************************************************
  static int 
  cmp_starts(const void *a, const void *b)
  {
    TimingType diff = ((CorrPhnRec *) a)->mpNode->mC.Start()
                    - ((CorrPhnRec *) b)->mpNode->mC.Start();
  
    if (diff != 0)
      return diff;
  
    return ((CorrPhnRec *) a)->mpNode->mC.Stop()
         - ((CorrPhnRec *) b)->mpNode->mC.Stop();
  }

  //***************************************************************************
  //***************************************************************************
  static int 
  cmp_overlapped(const void *key, const void *elem)
  {
    CorrPhnRec *corr_phn = (CorrPhnRec *) elem;
  
    if (((DecoderNetwork::NodeType *) key)->mC.Start() >=
          corr_phn->maxStopTimeTillNow) 
    {
      return 1;
    }
    if (corr_phn->mId == 0 || // first element in the array
      ((DecoderNetwork::NodeType *) key)->mC.Start() >= 
       (corr_phn-1)->maxStopTimeTillNow)
    {
      return  0;
    }
    else 
    {
      return -1;
    }
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
  void 
  DecoderNetwork::
  ComputePhoneCorrectnes(PhoneCorrectnessApproximationType approxType, MyHSearchData *silencePhones)
  {
//    NodeType*    node;
    CorrPhnRec*  corr_phn;
    CorrPhnRec*  overlapped;
    int          ncorr_phns = 0, i = 0;
    TimingType   maxStopTime;
    ENTRY        e; // = {0}; //{0} is just to make compiler happy
    ENTRY*       ep;

  
    iterator     p_node;

    for (p_node = begin(); p_node != end(); ++p_node) 
    {
      if (p_node->mC.mType & NT_PHONE && p_node->mC.mType & NT_TRUE) 
        ncorr_phns++;
    }

    if (ncorr_phns == 0) 
      Error("No correct phoneme node in network");
  
    corr_phn = new CorrPhnRec[ncorr_phns];

    if (corr_phn == NULL) 
      Error("Insufficient memory");
  
    for (p_node = begin(); p_node != end(); ++p_node)
    { 
      if (p_node->mC.mType & NT_PHONE && p_node->mC.mType & NT_TRUE) 
      {
        corr_phn[i++].mpNode = &(*p_node);
      }
    }
  
    qsort(corr_phn, ncorr_phns, sizeof(CorrPhnRec), cmp_starts);
  
    maxStopTime = corr_phn[0].mpNode->mC.Stop();

    for (i = 0; i < ncorr_phns; i++) 
    {
      corr_phn[i].mId = i;
      maxStopTime = HIGHER_OF(maxStopTime, corr_phn[i].mpNode->mC.Stop());
      corr_phn[i].maxStopTimeTillNow = maxStopTime;
    }

    //for (node = pFirstNode; node != NULL; node = node->mpNext) 
    for (p_node = begin(); p_node != end(); ++p_node)
    {
      if (!(p_node->mC.mType & NT_PHONE)) 
        continue;
  
      
      if (p_node->mC.Stop()  <= p_node->mC.Start())
      {    
        p_node->mC.SetPhoneAccuracy(0.0);
        continue;
      }

      e.key = p_node->mC.mpName;
      my_hsearch_r(e, FIND, &ep, silencePhones);
      if (ep != NULL) { // Silence phone
        //!!! List of ignored phonemes should be provided by some switch !!!
        if(approxType == MPE_FrameError) {
          // Silences are not considered as incorect
          p_node->mC.SetPhoneAccuracy((p_node->mC.Stop() - p_node->mC.Start()) / 100000.0);
        } else {      // MPE_ApproximateAccuracy, MFPE_FrameAccuracy
          p_node->mC.SetPhoneAccuracy(0.0);
        }
        continue;    
      }
      
      if(approxType == MPE_ApproximateAccuracy) {
        p_node->mC.SetPhoneAccuracy(-1.0);
      } else { // MPE_FrameAccuracy, MFPE_FrameError
        p_node->mC.SetPhoneAccuracy(0.0);
      }
        
      // Find first overlapping reference phone
      overlapped = (CorrPhnRec*) bsearch(&(*p_node), corr_phn, ncorr_phns,
                                         sizeof(CorrPhnRec), cmp_overlapped);
  
      if (overlapped) 
      {
        TimingType maxCorrectPhoneStopTimeTillNow = p_node->mC.Start();
        
        // Go over all reference phones that overlap with curent phone 'p_node'
        for (; overlapped < corr_phn + ncorr_phns &&
             overlapped->mpNode->mC.Start() < p_node->mC.Stop(); overlapped++) 
        {
          if (overlapped->mpNode->mC.Stop()  <= overlapped->mpNode->mC.Start() 
          ||  overlapped->mpNode->mC.Stop()  <= p_node->mC.Start()) {
            continue;
          }
          if(approxType == MPE_ApproximateAccuracy) {
            p_node->mC.SetPhoneAccuracy(
              HIGHER_OF(p_node->mC.PhoneAccuracy(), 
                        (SamePhoneme(overlapped->mpNode->mC.mpName, p_node->mC.mpName) + 1.0) 
                        * (LOWER_OF(overlapped->mpNode->mC.Stop(), p_node->mC.Stop()) 
                           -  HIGHER_OF(overlapped->mpNode->mC.Start(), p_node->mC.Start())) 
                        / (overlapped->mpNode->mC.Stop() - overlapped->mpNode->mC.Start()) - 1.0));
  
            if (p_node->mC.PhoneAccuracy() >= 1.0) {
              break;
            }
          } else { // MPE_FrameAccuracy, MFPE_FrameError
            if(SamePhoneme(overlapped->mpNode->mC.mpName, p_node->mC.mpName)) {
              FLOAT startingPhoneAccuracu = p_node->mC.PhoneAccuracy();
              
              p_node->mC.SetPhoneAccuracy(p_node->mC.PhoneAccuracy() + 
                  (LOWER_OF(overlapped->mpNode->mC.Stop(), p_node->mC.Stop()) - 
                   HIGHER_OF(overlapped->mpNode->mC.Start(), maxCorrectPhoneStopTimeTillNow)) / 100000.0);

              maxCorrectPhoneStopTimeTillNow = HIGHER_OF(maxCorrectPhoneStopTimeTillNow, 
                                                           overlapped->mpNode->mC.Stop());

              assert((p_node->mC.PhoneAccuracy() - startingPhoneAccuracu) <= (overlapped->mpNode->mC.Stop() - overlapped->mpNode->mC.Start())/ 100000.0);
              assert(p_node->mC.PhoneAccuracy() <= (p_node->mC.Stop() - p_node->mC.Start())/ 100000.0);
              
              if(maxCorrectPhoneStopTimeTillNow >= p_node->mC.Stop()) {
                break;
              }
            }
          }
        }
      }
    }
    delete []  corr_phn;
  }  

  //***************************************************************************
  //***************************************************************************
    int 
    DecoderNetwork::
    LatticeLocalOptimization_BackwardPass(const ExpansionOptions &expOptions,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale)
    {
      int     node_removed;

      Reverse();
      node_removed = LatticeLocalOptimization_ForwardPass(expOptions, 
          wordPenalty, modelPenalty, lmScale, posteriorScale);
      Reverse();

      return node_removed;
    }
  //  LatticeLocalOptimization_BackwardPass(const ExpansionOptions &expOptions)
  //***************************************************************************


/*
  void 
  FreeNetwork(DecoderNetwork::NodeType * pNode, bool compactRepresentation) 
  {
    if (!compactRepresentation)
    {
      DecoderNetwork::NodeType*  tnode;
      while (pNode) 
      {
        tnode = pNode->mpNext;
        delete pNode->mC.mpAlphaBeta;
        pNode->mC.mpAlphaBeta = NULL;
        free(pNode->rpLinks());
        free(pNode->rpBackLinks());
        free(pNode);
        pNode = tnode;
      }
    }
    else
    {
      NodeBasic<NodeBasicContent, LinkContent, LinkArray> *  tnode;
      for(tnode =  reinterpret_cast<NodeBasic<NodeBasicContent, LinkContent, LinkArray>* >(pNode); tnode->rNLinks() != 0; tnode++) 
        free(tnode->rpLinks());
      
      free(pNode);
    }
  }
*/

  

#ifndef NDEBUG
  // Debug function showing network using AT&T dot utility
  void dnet(DecoderNetwork::NodeType *net, int nAuxNodePtrs, ...)
  {
    static int dnetcnt=1;
    va_list ap;
    DecoderNetwork::NodeType* node;
    int i = 1;
  
    FILE *fp = popen("cat | (tf=`mktemp /tmp/netps.XXXXXX`;"
                    "dot -Tps > $tf; gv -scale -4 $tf; rm $tf)",
                    "w");
  //  FILE *fp = stdout;
  
    if (fp == NULL) return;
  
    for (node = net; node != NULL; node = node->mpNext) {
      node->mC.mEmittingStateId = i++;
    }
    fprintf(fp, "digraph \"dnet%d\" {\nrankdir=LR\n", dnetcnt++);
  
  
    for (node = net; node != NULL; node = node->mpNext) {
      fprintf(fp, "n%d [shape=%s,label=\"%d:%s", node->mC.mEmittingStateId,
              node->mC.mType & NT_WORD ? "box" : "ellipse", node->mC.mEmittingStateId,
              node->mC.mType & NT_WORD ? (node->mC.mpPronun ?
                                      node->mC.mpPronun->mpWord->mpName : "-"):
              node->mC.mType & NT_PHONE? node->mC.mpName :
              node->mC.mType & NT_MODEL? node->mC.mpHmm->mpMacro->mpName : "???");
  
      if (node->mC.mType & NT_WORD && node->mC.mpPronun != NULL) {
        if (node->mC.mpPronun != node->mC.mpPronun->mpWord->pronuns[0]) {
          fprintf(fp, ":%d", node->mC.mpPronun->variant_no);
        }
        fprintf(fp, "\\n");
  
        if (node->mC.mpPronun->outSymbol != node->mC.mpPronun->mpWord->mpName) {
          fprintf(fp, "[%s]", node->mC.mpPronun->outSymbol ?
                              node->mC.mpPronun->outSymbol : "");
        }
        if (node->mC.mpPronun->prob != 0.0) {
          fprintf(fp, " "FLOAT_FMT, node->mC.mpPronun->prob);
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
      for (i = 0; i < node->rNLinks(); i++) {
        fprintf(fp,"n%d -> n%d [color=blue,weight=1",
                node->mC.mEmittingStateId,node->rpLinks()[i].pNode()->mC.mEmittingStateId);
        if (node->rpLinks()[i].LmLike() != 0.0) {
          fprintf(fp,",label=\""FLOAT_FMT"\"", node->rpLinks()[i].LmLike());
        }
        fprintf(fp,"];\n");
      }
  //    for (i = 0; i < node->rNBackLinks(); i++) {
  //      fprintf(fp,"n%d -> n%d [color=red,weight=1",
  //              node->mEmittingStateId,node->rpBackLinks()[i].pNode()->mEmittingStateId);
  //      if (node->rpBackLinks()[i].mLmLike != 0.0) {
  //        fprintf(fp,",label=\""FLOAT_FMT"\"", node->rpBackLinks()[i].mLmLike);
  //      }
  //      fprintf(fp,"];\n");
  //    }
    }
    va_start(ap, nAuxNodePtrs);

    typedef  DecoderNetwork::NodeType  my_node;
    for (i = 0; i < nAuxNodePtrs; i++) 
    {
      my_node* ptr = va_arg(ap, my_node* );
      fprintf(fp, "AuxPtr%d [shape=plaintext];\nAuxPtr%d -> n%d\n",
              i, i, ptr->mC.mEmittingStateId);
    }
    va_end(ap);
  
    fprintf(fp, "}\n");
    pclose(fp);
  }
#endif

#include <iostream>
  
  //***************************************************************************
  //***************************************************************************
  // Remove null nones having less than three predecessors or less than three
  // successors
  int 
  DecoderNetwork::
  RemoveRedundantNullNodes(bool removeAllNullNodes,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale)
  {
    NodeType *    p_node;
    NodeType *    tnode;
    int       i;
    int       j;
    int       k;
    int       node_removed = 0;
  
    pFirst()->mpBackNext = NULL;

    for (p_node = pFirst(); p_node->mpNext != NULL; p_node = p_node->mpNext) 
    {
      p_node->mpNext->mpBackNext = p_node;
    }

    for (p_node = pFirst(); p_node != NULL; p_node = p_node->mpNext) 
    {
      if (p_node->mC.mType & NT_WORD && p_node->mC.mpPronun == NULL
      &&  p_node->rNLinks() != 0 && p_node->rNBackLinks() != 0
      &&  (removeAllNullNodes
          ||  p_node->rNLinks() == 1 
          ||  p_node->rNBackLinks() == 1 
          || (p_node->rNLinks() == 2 && p_node->rNBackLinks() == 2)))
      {
  
      node_removed = 1;
  
      // Remove links to current node form backlinked nodes and realloc
      // link arrays of backlinked nodes to hold p_node->rNLinks() more backlinks
        for (i = 0; i < p_node->rNBackLinks(); i++) 
        {
          NodeType* bakcnode = p_node->rpBackLinks()[i].pNode();

          for (j=0; j<bakcnode->NLinks() && bakcnode->rpLinks()[j].pNode()!=p_node; j++)
          {}

          assert(j < bakcnode->NLinks()); // Otherwise link to 'p_node' is missing
                                        // from which backlink exists
          bakcnode->rpLinks()[j] = bakcnode->rpLinks()[bakcnode->NLinks()-1];
  
          bakcnode->rpLinks() = (LinkType *) 
            realloc(bakcnode->rpLinks(), (bakcnode->NLinks() - 1 + 
                  p_node->NLinks()) * sizeof(LinkType));

          if (bakcnode->rpLinks() == NULL) 
            Error("Insufficient memory");

          bakcnode->rNLinks()--;// += word->npronuns-1;
        }
  
        // Remove backlinks to current node form linked nodes and realloc
        // backlink arrays of linked nodes to hold word->npronuns more backlinks
        for (i=0; i < p_node->NLinks(); i++) {
          NodeType *forwnode = p_node->rpLinks()[i].pNode();
          for (j=0;j<forwnode->rNBackLinks()&&forwnode->rpBackLinks()[j].pNode()!=p_node;j++);
          assert(j < forwnode->rNBackLinks());
          // Otherwise link to 'node' is missing from which backlink exists
          forwnode->rpBackLinks()[j] = forwnode->rpBackLinks()[forwnode->rNBackLinks()-1];
  
          forwnode->rpBackLinks() = (LinkType *)
            realloc(forwnode->rpBackLinks(),
                  (forwnode->rNBackLinks() - 1 + p_node->rNBackLinks()) * sizeof(LinkType));
          if (forwnode->rpBackLinks() == NULL) Error("Insufficient memory");
          forwnode->rNBackLinks()--;
        }
        for (j = 0; j < p_node->rNBackLinks(); j++) {
          NodeType *backnode = p_node->rpBackLinks()[j].pNode();
          int orig_nlinks = backnode->NLinks();
  
          for (i=0; i < p_node->NLinks(); i++) 
          {
            for(k = 0; k < orig_nlinks && backnode->rpLinks()[k].pNode() != p_node->rpLinks()[i].pNode(); k++);

            if(k < orig_nlinks) 
            {
              // Link which is to be created already exists. Its duplication must be avoided.
              if (backnode->rpLinks()[k].Like(lmScale) >  
                 p_node->rpLinks()[i].Like(lmScale) + p_node->rpBackLinks()[j].Like(lmScale))
              {
                backnode->rpLinks()[k].SetLmLike(backnode->rpLinks()[k].LmLike());
                backnode->rpLinks()[k].SetAcousticLike(backnode->rpLinks()[k].AcousticLike());
              }
              else
              {
                backnode->rpLinks()[k].SetLmLike(p_node->rpLinks()[i].LmLike() + p_node->rpBackLinks()[j].LmLike());
                backnode->rpLinks()[k].SetAcousticLike( p_node->rpLinks()[i].AcousticLike() + p_node->rpBackLinks()[j].AcousticLike());
              }
              /*
              backnode->rpLinks()[k].SetLmLike(HIGHER_OF(backnode->rpLinks()[k].LmLike(), 
                                             p_node->rpLinks()[i].LmLike() + p_node->rpBackLinks()[j].LmLike()));

              backnode->rpLinks()[k].SetAcousticLike(HIGHER_OF(backnode->rpLinks()[k].AcousticLike(), 
                                             p_node->rpLinks()[i].AcousticLike() + p_node->rpBackLinks()[j].AcousticLike()));
                                             */
            } 
            else 
            {
              backnode->rpLinks()[backnode->NLinks()].SetNode(p_node->rpLinks()[i].pNode());

              backnode->rpLinks()[backnode->NLinks()].SetLmLike(
                  p_node->rpLinks()[i].LmLike() + p_node->rpBackLinks()[j].LmLike());

              backnode->rpLinks()[backnode->NLinks()].SetAcousticLike(
                  p_node->rpLinks()[i].AcousticLike() + p_node->rpBackLinks()[j].AcousticLike());

              ++(backnode->rNLinks());
            }
          }
        }

        for (j = 0; j < p_node->NLinks(); j++) 
        {
          NodeType *forwnode = p_node->rpLinks()[j].pNode();
          int orig_nbacklinks = forwnode->rNBackLinks();
  
          for (i=0; i < p_node->rNBackLinks(); i++) 
          {
            for(k = 0; k < orig_nbacklinks && forwnode->rpBackLinks()[k].pNode() != p_node->rpBackLinks()[i].pNode(); k++);

            if (k < orig_nbacklinks) 
            {
              if (forwnode->rpBackLinks()[k].Like(lmScale) > 
                p_node->rpBackLinks()[i].Like(lmScale) + p_node->rpLinks()[j].Like(lmScale))
              {
                forwnode->rpBackLinks()[k].SetLmLike(forwnode->rpBackLinks()[k].LmLike()); 
                forwnode->rpBackLinks()[k].SetAcousticLike(forwnode->rpBackLinks()[k].AcousticLike()); 
              }
              else
              {
                forwnode->rpBackLinks()[k].SetLmLike(p_node->rpBackLinks()[i].LmLike() + p_node->rpLinks()[j].LmLike());
                forwnode->rpBackLinks()[k].SetAcousticLike(p_node->rpBackLinks()[i].AcousticLike() + p_node->rpLinks()[j].AcousticLike());
              }

              // Link which is to be created already exists. Its duplication must be avoided.
              /*
              forwnode->rpBackLinks()[k].SetLmLike(HIGHER_OF(forwnode->rpBackLinks()[k].LmLike(), 
                                                      p_node->rpBackLinks()[i].LmLike() + p_node->rpLinks()[j].LmLike()));

              forwnode->rpBackLinks()[k].SetAcousticLike(HIGHER_OF(forwnode->rpBackLinks()[k].AcousticLike(), 
                                                      p_node->rpBackLinks()[i].AcousticLike() + p_node->rpLinks()[j].AcousticLike()));
                                                      */
            } 
            else 
            {
              forwnode->rpBackLinks()[forwnode->rNBackLinks()].SetNode(p_node->rpBackLinks()[i].pNode());

              forwnode->rpBackLinks()[forwnode->rNBackLinks()].SetLmLike(
                  p_node->rpBackLinks()[i].LmLike() + p_node->rpLinks()[j].LmLike());

              forwnode->rpBackLinks()[forwnode->rNBackLinks()].SetAcousticLike(
                  p_node->rpBackLinks()[i].AcousticLike() + p_node->rpLinks()[j].AcousticLike());

              forwnode->rNBackLinks()++;
            }
          }
        }
        p_node->mpBackNext->mpNext = p_node->mpNext;
        p_node->mpNext->mpBackNext = p_node->mpBackNext;
        tnode = p_node;
        p_node = p_node->mpBackNext;
        delete tnode->mC.mpAlphaBeta;
        tnode->mC.mpAlphaBeta = NULL;
        free(tnode->rpLinks());
        free(tnode->rpBackLinks());
        free(tnode);
      }
    }
    return node_removed;
  }
  //  RemoveRedundantNullNodes(NodeType<NODE_REGULAR, LINK_REGULAR> *pFirstNode)
  //****************************************************************************


  //***************************************************************************
  //***************************************************************************
  void 
  DecoderNetwork::
  SelfLinksToNullNodes()
  {
    int   i;
    int   j;

    //NodeType* node;
    NodeType* tnode;
    iterator   p_node;
    
    for (p_node = begin(); p_node != end(); p_node++) 
    {
      NodeType* p_node_real_address = &(*p_node);

      for (i=0; i < p_node->NLinks(); i++) 
      {
        if (p_node->rpLinks()[i].pNode() == &(*p_node))
        {
          if ((tnode           = (NodeType *) calloc(1, sizeof(NodeType))) == NULL ||
              (tnode->rpLinks()     = (LinkType *) malloc(sizeof(LinkType))) == NULL ||
              (tnode->rpBackLinks() = (LinkType *) malloc(sizeof(LinkType))) == NULL) 
          {
            Error("Insufficient memory");
          }
  
          tnode->mC.mpName = NULL;
          p_node->rpLinks()[i].SetNode(tnode);
          
          for (j=0; j<p_node->rNBackLinks() && p_node->rpBackLinks()[j].pNode() != 
                    p_node_real_address ; j++)
          {}
          
          assert(j<p_node->rNBackLinks());
          
          p_node->rpBackLinks()[j].SetNode(tnode);
          p_node->rpBackLinks()[j].SetLmLike(0.0);
          p_node->rpBackLinks()[j].SetAcousticLike(0.0);
  
          tnode->mC.mType       = NT_WORD;
          tnode->mC.mpPronun     = NULL;
          tnode->rNLinks()     = 1;
          tnode->rNBackLinks() = 1;
          tnode->mC.SetStart(UNDEF_TIME);
          tnode->mC.SetStop(UNDEF_TIME);
  //        tnode->mpTokens     = NULL;
  //        tnode->mpExitToken  = NULL;
          tnode->rpLinks()[0].SetNode(p_node_real_address);
          tnode->rpLinks()[0].SetLmLike(0.0);
          tnode->rpLinks()[0].SetAcousticLike(0.0);
          tnode->rpBackLinks()[0].SetNode(p_node_real_address);
          tnode->rpBackLinks()[0].SetLmLike(p_node->rpLinks()[i].LmLike());
          tnode->rpBackLinks()[0].SetAcousticLike(p_node->rpLinks()[i].AcousticLike());
          tnode->mpNext = p_node->mpNext;
          p_node->mpNext = tnode;
        }
      }
    }
    /*
    for (node = pFirst(); node != NULL; node = node->mpNext) 
    {
      for (i=0; i < node->NLinks(); i++) 
      {
        if (node->rpLinks()[i].pNode() == node) 
        {
          if ((tnode           = (NodeType *) calloc(1, sizeof(NodeType))) == NULL ||
              (tnode->rpLinks()     = (LinkType *) malloc(sizeof(LinkType))) == NULL ||
              (tnode->rpBackLinks() = (LinkType *) malloc(sizeof(LinkType))) == NULL) 
          {
            Error("Insufficient memory");
          }
  
          tnode->mpName = NULL;
          node->rpLinks()[i].pNode() = tnode;
          
          for (j=0; j<node->rNBackLinks() && node->rpBackLinks()[j].pNode()!=node; j++)
          {}
          
          assert(j<node->rNBackLinks());
          
          node->rpBackLinks()[j].pNode() = tnode;
          node->rpBackLinks()[j].mLmLike = 0.0;
          node->rpBackLinks()[j].mAcousticLike = 0.0;
  
          tnode->mType       = NT_WORD;
          tnode->mpPronun     = NULL;
          tnode->rNLinks()     = 1;
          tnode->rNBackLinks() = 1;
          tnode->Start()      = UNDEF_TIME;
          tnode->mStop       = UNDEF_TIME;
  //        tnode->mpTokens     = NULL;
  //        tnode->mpExitToken  = NULL;
          tnode->rpLinks()[0].mpNode     = node;
          tnode->rpLinks()[0].mLmLike     = 0.0;
          tnode->rpLinks()[0].mAcousticLike     = 0.0;
          tnode->rpBackLinks()[0].mpNode = node;
          tnode->rpBackLinks()[0].mLmLike = node->rpLinks()[i].mLmLike;
          tnode->rpBackLinks()[0].mAcousticLike = node->rpLinks()[i].mAcousticLike;
          tnode->mpNext = node->mpNext;
          node->mpNext = tnode;
        }
      }
    }
    */
  }
  // SelfLinksToNullNodes()
  //****************************************************************************




}
// namespace STK
// ****************************************************************************

