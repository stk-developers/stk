#include <cstdlib>


namespace STK
{
  //***************************************************************************
  //***************************************************************************
  static int 
  lnkcmp(const void *a, const void *b)
  {
  //  return ((link_type *) a)->mpNode - ((link_type *) b)->mpNode;
  //  Did not work with gcc, probably bug in gcc pointer arithmetic
    return (char *)((Link<LINK_BASIC> *) a)->mpNode - (char *)((Link<LINK_BASIC> *) b)->mpNode;
  }

  
  //***************************************************************************
  //***************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    Release() 
    {
      if (IsEmpty()) 
        return;

      if (!mCompactRepresentation)
      {
        node_type*  p_tmp_node;
        node_type*  p_node(pFirst());
        
        while (p_node) 
        {
          p_tmp_node = p_node->mpNext;
          free(p_node->mpLinks);
          free(p_node->mpBackLinks);
          free(p_node);
          p_node = p_tmp_node;
        }
      }
      else
      {
        NodeBasic*  p_node(pFirst());
        NodeBasic*  p_tmp_node;

        for(p_tmp_node =  reinterpret_cast<NodeBasic*>(p_node); 
            p_tmp_node->mNLinks != 0; 
            p_tmp_node++) 
        {
          free(p_tmp_node->mpLinks);
        }
        
        free(p_node);
      }

      SetFirst(NULL);
    }


  //***************************************************************************
  //***************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    BuildFromLabels(const Label* pLabels, NodeType nodeType)
    {
      const Label*    p_lp;
      node_type*      p_first;
      node_type*      p_last = NULL;
      node_type*      p_node;
    
      // allocate the memory
      if ((p_first             = (node_type *) calloc(1, sizeof(node_type))) == NULL ||
          (p_last              = (node_type *) calloc(1, sizeof(node_type))) == NULL ||
          (p_first->mpLinks    = (link_type *) malloc(sizeof(link_type))) == NULL    ||
          (p_last->mpBackLinks = (link_type *) malloc(sizeof(link_type))) == NULL) 
      {
        Error("Insufficient memory");
      }

      p_first->mpName       = p_last->mpName      = NULL;     
      p_first->mType        = p_last->mType       = NT_WORD;
      p_first->mpPronun     = p_last->mpPronun    = NULL;
      p_first->mNLinks      = p_last->mNBackLinks = 1;
      p_first->mNBackLinks  = p_last->mNLinks     = 0;
      p_first->mpBackLinks  = p_last->mpLinks     = NULL;
      p_first->mStart       = p_last->mStart      = UNDEF_TIME;
      p_first->mStop        = p_last->mStop       = UNDEF_TIME;
    //  p_first->mpTokens        = p_last->mpTokens        = NULL;
    //  p_first->mpExitToken     = p_last->mpExitToken     = NULL;
    
      p_node = p_first;
      
      for (p_lp = pLabels; p_lp != NULL; p_lp = p_lp->mpNext) 
      {
        node_type*  p_tmp_node;
    
        if ((p_tmp_node              = (node_type *) calloc(1, sizeof(node_type))) == NULL ||
            (p_tmp_node->mpLinks     = (link_type *) malloc(sizeof(link_type))) == NULL    ||
            (p_tmp_node->mpBackLinks = (link_type *) malloc(sizeof(link_type))) == NULL) 
        {
          Error("Insufficient memory");
        }
        
        p_tmp_node->mpName = NULL;
        
        p_node->mpLinks[0].mpNode = p_tmp_node;
        p_node->mpLinks[0].mLmLike  = 0.0;
        p_node->mpLinks[0].mAcousticLike  = 0.0;
        
        switch (nodeType) 
        {
          case NT_WORD:  p_tmp_node->mpPronun = ((Word *) p_lp->mpData)->pronuns[0]; break;
          case NT_MODEL: p_tmp_node->mpHmm    =   (Hmm *) p_lp->mpData;              break;
          case NT_PHONE: p_tmp_node->mpName   =  (char *) p_lp->mpData;              break;
          default:       Error("Fatal: Invalid node type");
        }
        
        p_tmp_node->mType       = nodeType;
        p_tmp_node->mNLinks     = 1;
        p_tmp_node->mNBackLinks = 1;
        p_tmp_node->mStart      = p_lp->mStart;
        p_tmp_node->mStop       = p_lp->mStop;
        p_tmp_node->mpBackLinks[0].mpNode = p_node;
        p_tmp_node->mpBackLinks[0].mLmLike   = 0.0;
        p_tmp_node->mpBackLinks[0].mAcousticLike   = 0.0;
        p_node->mpNext = p_tmp_node;
        p_node = p_tmp_node;
      }
      
      p_node->mpNext = p_last;
      p_node->mpLinks[0].mpNode    = p_last;
      p_node->mpLinks[0].mLmLike    = 0.0;
      p_node->mpLinks[0].mAcousticLike    = 0.0;
      p_last->mpBackLinks[0].mpNode = p_node;
      p_last->mpBackLinks[0].mLmLike = 0.0;
      p_last->mpBackLinks[0].mAcousticLike = 0.0;
      p_last->mpNext = NULL;
    }
  // Network::MakeNetworkFromLabels(Label * pLabels, NodeType nodeType)
  //***************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    SelfLinksToNullNodes()
    {
      int   i;
      int   j;

      node_type* node;
      node_type* tnode;
      iterator   p_node;
      
      for (p_node = begin(); p_node != end(); p_node++) 
      {
        for (i=0; i < node->mNLinks; i++) 
        {
          if (p_node->mpLinks[i].mpNode == p_node) 
          {
            if ((tnode           = (node_type *) calloc(1, sizeof(node_type))) == NULL ||
                (tnode->mpLinks     = (link_type *) malloc(sizeof(link_type))) == NULL ||
                (tnode->mpBackLinks = (link_type *) malloc(sizeof(link_type))) == NULL) 
            {
              Error("Insufficient memory");
            }
    
            tnode->mpName = NULL;
            node->mpLinks[i].mpNode = tnode;
            
            for (j=0; j<node->mNBackLinks && node->mpBackLinks[j].mpNode!=node; j++)
            {}
            
            assert(j<node->mNBackLinks);
            
            node->mpBackLinks[j].mpNode = tnode;
            node->mpBackLinks[j].mLmLike = 0.0;
            node->mpBackLinks[j].mAcousticLike = 0.0;
    
            tnode->mType       = NT_WORD;
            tnode->mpPronun     = NULL;
            tnode->mNLinks     = 1;
            tnode->mNBackLinks = 1;
            tnode->mStart      = UNDEF_TIME;
            tnode->mStop       = UNDEF_TIME;
    //        tnode->mpTokens     = NULL;
    //        tnode->mpExitToken  = NULL;
            tnode->mpLinks[0].mpNode     = node;
            tnode->mpLinks[0].mLmLike     = 0.0;
            tnode->mpLinks[0].mAcousticLike     = 0.0;
            tnode->mpBackLinks[0].mpNode = node;
            tnode->mpBackLinks[0].mLmLike = node->mpLinks[i].mLmLike;
            tnode->mpBackLinks[0].mAcousticLike = node->mpLinks[i].mAcousticLike;
            tnode->mpNext = node->mpNext;
            node->mpNext = tnode;
          }
        }
      }
      /*
      for (node = pFirst(); node != NULL; node = node->mpNext) 
      {
        for (i=0; i < node->mNLinks; i++) 
        {
          if (node->mpLinks[i].mpNode == node) 
          {
            if ((tnode           = (node_type *) calloc(1, sizeof(node_type))) == NULL ||
                (tnode->mpLinks     = (link_type *) malloc(sizeof(link_type))) == NULL ||
                (tnode->mpBackLinks = (link_type *) malloc(sizeof(link_type))) == NULL) 
            {
              Error("Insufficient memory");
            }
    
            tnode->mpName = NULL;
            node->mpLinks[i].mpNode = tnode;
            
            for (j=0; j<node->mNBackLinks && node->mpBackLinks[j].mpNode!=node; j++)
            {}
            
            assert(j<node->mNBackLinks);
            
            node->mpBackLinks[j].mpNode = tnode;
            node->mpBackLinks[j].mLmLike = 0.0;
            node->mpBackLinks[j].mAcousticLike = 0.0;
    
            tnode->mType       = NT_WORD;
            tnode->mpPronun     = NULL;
            tnode->mNLinks     = 1;
            tnode->mNBackLinks = 1;
            tnode->mStart      = UNDEF_TIME;
            tnode->mStop       = UNDEF_TIME;
    //        tnode->mpTokens     = NULL;
    //        tnode->mpExitToken  = NULL;
            tnode->mpLinks[0].mpNode     = node;
            tnode->mpLinks[0].mLmLike     = 0.0;
            tnode->mpLinks[0].mAcousticLike     = 0.0;
            tnode->mpBackLinks[0].mpNode = node;
            tnode->mpBackLinks[0].mLmLike = node->mpLinks[i].mLmLike;
            tnode->mpBackLinks[0].mAcousticLike = node->mpLinks[i].mAcousticLike;
            tnode->mpNext = node->mpNext;
            node->mpNext = tnode;
          }
        }
      }
      */
    }
  // SelfLinksToNullNodes()
  //****************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    LatticeLocalOptimization(int strictTiming, int trace_flag)
    {
      node_type *    node;
      node_type *    lastnode;
      int            i;
      int            j;
      int            unreachable = 0;
      
      // For each node, sort links by pointer value to allow
      // for easy comparison whether two nodes have the same set of links
      for (node = pFirst(); node != NULL; node = node->mpNext)  
      {
        node->mAux = 0;
        node->mpBackNext = node->mpNext;
        qsort(node->mpLinks, node->mNLinks, sizeof(link_type), lnkcmp);
        qsort(node->mpBackLinks, node->mNBackLinks, sizeof(link_type), lnkcmp);
      }
    
      // Sort nodes in topological order
      // printf("Sorting nodes...\n");
      pFirst()->mAux = 1;
      for (lastnode = node = pFirst(); node != NULL; node = node->mpNext) 
      {
        for (i=0; i < node->mNLinks; i++) 
        {
          node_type *lnknode = node->mpLinks[i].mpNode;
          
          if (lnknode->mAux == 0) 
          {
            for (j=0; j<lnknode->mNBackLinks && lnknode->mpBackLinks[j].mpNode->mAux==1; j++)
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
        for (node = pFirst(); node != NULL; node = node->mpBackNext) 
        {
          node->mAux = 0;
        }
          
        pFirst()->mAux = 1;
        for (lastnode = node = pFirst(); node != NULL; node = node->mpNext) 
        {
          for (i=0; i < node->mNLinks; i++) 
          {
            node_type *lnknode = node->mpLinks[i].mpNode;
            if (lnknode->mAux == 0) 
            {
              lastnode->mpNext = lnknode;
              lastnode  = lnknode;
              lnknode->mAux = 1;
              lnknode->mpNext = NULL;
            }
          }
        }
        
        for (node=pFirst(); node->mpNext->mNLinks != 0; node=node->mpNext)
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
          node_type *tnode = node->mpBackNext;
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
    
      pFirst()->mpBackNext = NULL;
      
      for (node=pFirst(); node->mpNext != NULL; node=node->mpNext) 
      {
        node->mpNext->mpBackNext = node;
      }
      
      for (i=1, node=pFirst(); node != NULL; node = node->mpNext, i++) 
      {
        node->mAux=i;
      }
      
      for (;;) 
      {
        if (trace_flag & 2) 
        {
          for (i=0,node=pFirst(); node; node=node->mpNext,i++)
          {}
          
          TraceLog("Forward pass.... (number of nodes: %d)", i);
        }
        
        LatticeLocalOptimization_ForwardPass(pFirst(), strictTiming);
    
        if (trace_flag & 2) 
        {
          for (i=0,node=pFirst(); node; node=node->mpNext,i++)
          {}
          
          TraceLog("Backward pass... (number of nodes: %d)", i);
        }
        
        if (!LatticeLocalOptimization_BackwardPass(pFirst(), strictTiming)) 
          break;
      }
    }  
  // LatticeLocalOptimization(int strictTiming, int trace_flag)
  //****************************************************************************


  //****************************************************************************
  //****************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    ExpandByDictionary(
      MyHSearchData* pDict,
      bool keep_word_nodes,
      bool multiple_pronun)
    {
      node_type* node; 
      node_type* prev = NULL;
      int   i;
      int   j;
    
      Pronun    singlePronun;
      Pronun*   singlePronunPtr;
      Word      singlePronunWrd;

      singlePronunWrd.npronuns = 1;
      singlePronunWrd.pronuns  = &singlePronunPtr;
      singlePronunPtr = &singlePronun;
    
      assert(!IsEmpty() || pFirst()->mType & NT_WORD || pFirst()->mpPronun == NULL);
    
      for (node = pFirst(); node != NULL; prev = node, node = node->mpNext) 
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
          node_type *bakcnode = node->mpBackLinks[i].mpNode;

          for (j=0; j<bakcnode->mNLinks && bakcnode->mpLinks[j].mpNode!=node; j++)
          {}

          assert(j < bakcnode->mNLinks); // Otherwise link to 'node' is missing
                                        // from which backlink exists
          bakcnode->mpLinks[j] = bakcnode->mpLinks[bakcnode->mNLinks-1];
    
          bakcnode->mpLinks = (link_type *)
            realloc(bakcnode->mpLinks,
                  (bakcnode->mNLinks - 1 + word->npronuns) * sizeof(link_type));
          if (bakcnode->mpLinks == NULL) Error("Insufficient memory");
          bakcnode->mNLinks--;// += word->npronuns-1;
        }
    
        // Remove backlinks to current node form linked nodes and realloc
        // backlink arrays of linked nodes to hold word->npronuns more backlinks
        for (i=0; i < node->mNLinks; i++) 
        {
          node_type* forwnode = node->mpLinks[i].mpNode;

          for (j=0;j<forwnode->mNBackLinks&&forwnode->mpBackLinks[j].mpNode!=node;j++)
          {}

          assert(j < forwnode->mNBackLinks);
          // Otherwise link to 'node' is missing from which backlink exists
          forwnode->mpBackLinks[j] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
    
          forwnode->mpBackLinks = (link_type *)
            realloc(forwnode->mpBackLinks,
                  (forwnode->mNBackLinks - 1 + word->npronuns) * sizeof(link_type));

          if (forwnode->mpBackLinks == NULL) 
            Error("Insufficient memory");

          forwnode->mNBackLinks--;
        }

        for (i = 0; i < word->npronuns; i++) 
        {
          Pronun *pronun = word->pronuns[i];
          node_type *pronun_first = NULL, *pronun_prev = NULL, *tnode;
    
          for (j = 0; j < pronun->nmodels; j++) 
          {
            tnode = (node_type *) calloc(1, sizeof(node_type));
            
            if (tnode == NULL) 
              Error("Insufficient memory");
    
            tnode->mType       = NT_PHONE | (node->mType & NT_TRUE);
            tnode->mpName      = pronun->model[j].mpName;
            tnode->mStart      = node->mStart;
            tnode->mStop       = node->mStop;
            tnode->mPhoneAccuracy = 1.0;
    
            if (j == 0) 
            {
              pronun_first = tnode;
            } 
            else 
            {
              if ((pronun_prev->mpLinks = (link_type *) malloc(sizeof(link_type))) == NULL ||
                  (tnode->mpBackLinks   = (link_type *) malloc(sizeof(link_type))) == NULL) 
              {
                Error("Insufficient memory");
              }
              
              tnode->mNBackLinks              = 1;
              tnode->mpBackLinks[0].mpNode    = pronun_prev;
              tnode->mpBackLinks[0].mLmLike     = 0.0;
              tnode->mpBackLinks[0].mAcousticLike     = 0.0;
              pronun_prev->mNLinks            = 1;
              pronun_prev->mpLinks[0].mpNode  = tnode;
              pronun_prev->mpLinks[0].mLmLike   = 0.0;
              pronun_prev->mpLinks[0].mAcousticLike   = 0.0;
              pronun_prev->mpNext             = tnode;
            }
            pronun_prev = tnode;
          }
          
          if (keep_word_nodes || j == 0) 
          {
            tnode = (node_type *) calloc(1, sizeof(node_type));
            
            if (tnode == NULL) Error("Insufficient memory");
    
            tnode->mpName   = NULL;
            tnode->mType    = NT_WORD | (node->mType & NT_TRUE);
            tnode->mpPronun = keep_word_nodes ? word->pronuns[i] : NULL;
            tnode->mStart   = node->mStart;
            tnode->mStop    = node->mStop;
    
            if (j == 0) {
              pronun_first = tnode;
            } else {
              if ((pronun_prev->mpLinks = (link_type *) malloc(sizeof(link_type))) == NULL ||
                (tnode->mpBackLinks   = (link_type *) malloc(sizeof(link_type))) == NULL) {
                Error("Insufficient memory");
              }
              tnode->mNBackLinks          = 1;
              tnode->mpBackLinks[0].mpNode   = pronun_prev;
              tnode->mpBackLinks[0].mLmLike   = 0.0;
              tnode->mpBackLinks[0].mAcousticLike   = 0.0;
              pronun_prev->mNLinks        = 1;
              pronun_prev->mpLinks[0].mpNode = tnode;
              pronun_prev->mpLinks[0].mLmLike = 0.0;
              pronun_prev->mpLinks[0].mAcousticLike = 0.0;
              pronun_prev->mpNext          = tnode;
            }
            pronun_prev = tnode;
          }
          if ((pronun_prev->mpLinks =
                (link_type *) malloc(sizeof(link_type) * node->mNLinks))==NULL ||
            (pronun_first->mpBackLinks =
                (link_type *) malloc(sizeof(link_type) * node->mNBackLinks)) == NULL) {
            Error("Insufficient memory");
          }
          pronun_prev->mNLinks      = node->mNLinks;
          pronun_first->mNBackLinks = node->mNBackLinks;
    
          for (j = 0; j < node->mNBackLinks; j++) {
            node_type *backnode = node->mpBackLinks[j].mpNode;
            backnode->mpLinks[backnode->mNLinks  ].mpNode = pronun_first;
            backnode->mpLinks[backnode->mNLinks++].mLmLike = node->mpBackLinks[j].mLmLike;
            backnode->mpLinks[backnode->mNLinks++].mAcousticLike = node->mpBackLinks[j].mAcousticLike;
            pronun_first->mpBackLinks[j] = node->mpBackLinks[j];
          }
          for (j=0; j < node->mNLinks; j++) {
            node_type *forwnode = node->mpLinks[j].mpNode;
            forwnode->mpBackLinks[forwnode->mNBackLinks  ].mpNode = pronun_prev;
            forwnode->mpBackLinks[forwnode->mNBackLinks++].mLmLike = node->mpLinks[j].mLmLike;
            forwnode->mpBackLinks[forwnode->mNBackLinks++].mAcousticLike = node->mpLinks[j].mAcousticLike;
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
  // ExpandWordNetworkByDictionary( ...
  //****************************************************************************


  //***************************************************************************
  //***************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    ExpandMonophonesToTriphones(MyHSearchData *nonCDphones, 
        MyHSearchData *CDphones)
    {
      node_type *  p_node;
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
        node_type*      prev = NULL;
        
        did_we_clone = 0;
        
        for (p_node = pFirst(); p_node != NULL; prev = p_node, p_node = p_node->mpNext) 
        {
          if ((p_node->mNLinks == 0 || p_node->mNBackLinks == 0) ||
              (p_node->mNLinks == 1 && p_node->mNBackLinks == 1)) 
          {
            continue;
          }
          
          if (p_node->mType & NT_PHONE) 
          {
            e.key = p_node->mpName;
            my_hsearch_r(e, FIND, &ep, nonCDphones);
            if (ep == NULL || !reinterpret_cast<size_t>(ep->data))
              continue; // Node is not a Tee model
          }
          
          did_we_clone = 1;
          assert(prev != NULL); //Otherwise pFirstNode node is not Null node
    
          // Remove links to current node form back-linked nodes and realloc
          // link arrays of back-linked nodes to hold node->mNLinks more links
          for (j=0; j < p_node->mNBackLinks; j++) 
          {
            node_type *backnode = p_node->mpBackLinks[j].mpNode;
            
            for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].mpNode!=p_node; k++)
            {}
            
            assert(k < backnode->mNLinks);
            
            // Otherwise link to 'node' is missing from which backlink exists
            backnode->mpLinks[k] = backnode->mpLinks[backnode->mNLinks-1];
    
            backnode->mpLinks = 
              (link_type *) realloc((backnode->mpLinks), 
                               (backnode->mNLinks-1+p_node->mNLinks)*sizeof(link_type));
            
            if (backnode->mpLinks == NULL) 
              Error("Insufficient memory");
            
            backnode->mNLinks--;
          }
          
          // Remove backlinks to current node form linked nodes and realloc
          // backlink arrays of linked nodes to hold node->mNBackLinks more backlinks
          for (j=0; j < p_node->mNLinks; j++) 
          {
            node_type *forwnode = p_node->mpLinks[j].mpNode;
            
            for (k=0;k<forwnode->mNBackLinks&&forwnode->mpBackLinks[k].mpNode!=p_node;k++)
            {}
            
            assert(k < forwnode->mNBackLinks);
            // Otherwise link to 'p_node' is missing from which backlink exists
            forwnode->mpBackLinks[k] = forwnode->mpBackLinks[forwnode->mNBackLinks-1];
    
            forwnode->mpBackLinks = 
              (link_type *) realloc((forwnode->mpBackLinks),
                               (forwnode->mNBackLinks-1+p_node->mNBackLinks)*sizeof(link_type));
                               
            if (forwnode->mpBackLinks == NULL) 
              Error("Insufficient memory");
              
            forwnode->mNBackLinks--;
          }
          // Alloc new p_node->mNLinks * p_node->mNBackLinks nodes and create new links
          // so that each backlinked node is conected with each linked node through
          // one new node.
          for (i=0; i < p_node->mNLinks; i++) 
          {
            for (j=0; j < p_node->mNBackLinks; j++) 
            {
              node_type *  tnode;
              link_type    forwlink = p_node->mpLinks[i];
              link_type    backlink = p_node->mpBackLinks[j];
    
              if ((tnode = (node_type *) calloc(1, sizeof(node_type))) == NULL)
                Error("Insufficient memory");
              
              tnode->mpName = NULL;
              *tnode = *p_node;
    
              if ((tnode->mpLinks     = (link_type *) malloc(sizeof(link_type))) == NULL ||
                (tnode->mpBackLinks = (link_type *) malloc(sizeof(link_type))) == NULL) 
              {
                Error("Insufficient memory");
              }
              
              tnode->mNLinks       = 1;
              tnode->mNBackLinks   = 1;
              tnode->mpLinks[0]     = forwlink;
              tnode->mpBackLinks[0] = backlink;
              
              forwlink.mpNode->mpBackLinks[forwlink.mpNode->mNBackLinks  ].mpNode = tnode;
              forwlink.mpNode->mpBackLinks[forwlink.mpNode->mNBackLinks++].mLmLike  = forwlink.mLmLike;
              forwlink.mpNode->mpBackLinks[forwlink.mpNode->mNBackLinks++].mAcousticLike  = forwlink.mAcousticLike;
              backlink.mpNode->mpLinks    [backlink.mpNode->mNLinks      ].mpNode = tnode;
              backlink.mpNode->mpLinks    [backlink.mpNode->mNLinks++    ].mLmLike  = backlink.mLmLike;
              backlink.mpNode->mpLinks    [backlink.mpNode->mNLinks++    ].mAcousticLike  = backlink.mAcousticLike;
              
              prev->mpNext = tnode;
              prev = tnode;
            }
          }
          prev->mpNext = p_node->mpNext;
          free(p_node->mpLinks);
          free(p_node->mpBackLinks);
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
      node_type *prev = NULL;
      for (p_node = pFirst(); p_node != NULL; prev = p_node, p_node = p_node->mpNext) 
      {
        ENTRY e = {0}; //{0} is just to make compiler happy
        ENTRY * ep;
    
        if ((p_node->mType & NT_WORD) ||
            (p_node->mNLinks == 1 && p_node->mNBackLinks == 1)) 
        {
          continue;
        }
        
        assert(p_node->mType & NT_PHONE);
        e.key = p_node->mpName;
        my_hsearch_r(e, FIND, &ep, nonCDphones);
        if (ep != NULL && reinterpret_cast<size_t>(ep->data)) continue; // Node is a Tee model
    
        assert(prev != NULL); //Otherwise first node is not Null node
    
        // Count groups of backlinked nodes corresponding to different monophones
        id = -1;
        nbackmononodes = 0;
        for (j=0; j < p_node->mNBackLinks; j++) 
        {
          if (p_node->mpBackLinks[j].mpNode->mAux != id) 
          {
            id = p_node->mpBackLinks[j].mpNode->mAux;
            nbackmononodes++;
          }
        }
        // Count groups of linked nodes corresponding to different monophones
        id = -1;
        nforwmononodes = 0;
        for (j=0; j < p_node->mNLinks; j++) 
        {
          if (p_node->mpLinks[j].mpNode->mAux != id) 
          {
            id = p_node->mpLinks[j].mpNode->mAux;
            nforwmononodes++;
          }
        }
    
        // Remove links to current node form backlinked nodes and realloc
        // link arrays of backlinked nodes to hold nforwmononodes more links
        for (j=0; j < p_node->mNBackLinks; j++) 
        {
          node_type *backnode = p_node->mpBackLinks[j].mpNode;
          for (k=0; k<backnode->mNLinks && backnode->mpLinks[k].mpNode!=p_node; k++);
          assert(k < backnode->mNLinks);
          // Otherwise link to 'node' is missing from which backlink exists
          memmove(backnode->mpLinks+k, backnode->mpLinks+k+1,
                  (backnode->mNLinks-k-1) * sizeof(link_type));
    
          backnode->mpLinks = (link_type *)
            realloc(backnode->mpLinks,
                  (backnode->mNLinks-1+nforwmononodes)*sizeof(link_type));
          if (backnode->mpLinks == NULL) Error("Insufficient memory");
          backnode->mNLinks--;
        }
    
        // Remove backlinks to current p_node form linked nodes and realloc
        // backlink arrays of linked nodes to hold nbackmononodes more backlinks
        for (j=0; j < p_node->mNLinks; j++) 
        {
          node_type *forwnode = p_node->mpLinks[j].mpNode;
          for (k=0;k<forwnode->mNBackLinks&&forwnode->mpBackLinks[k].mpNode!=p_node;k++);
          assert(k < forwnode->mNBackLinks);
          // Otherwise link to 'p_node' is missing from which backlink exists
          memmove(forwnode->mpBackLinks+k, forwnode->mpBackLinks+k+1,
                  (forwnode->mNBackLinks-k-1) * sizeof(link_type));
    
          forwnode->mpBackLinks = (link_type *)
            realloc(forwnode->mpBackLinks,
                  (forwnode->mNBackLinks-1+nbackmononodes)*sizeof(link_type));
          if (forwnode->mpBackLinks == NULL) Error("Insufficient memory");
          forwnode->mNBackLinks--;
        }
    
        // Alloc new nforwmononodes * nbackmononodes nodes and create new links
        // so that each backlinked node is conected through one new node with all
        // linked nodes belonging to one monophone group and vice versa each
        // linked node is conected through one new node with all backlinked nodes
        // belonging to one monophone group
        link_type *forwmono_start, *forwmono_end = p_node->mpLinks;
        for (i=0; i < nforwmononodes; i++) 
        {
          for (forwmono_start = forwmono_end;
              forwmono_end < p_node->mpLinks+p_node->mNLinks &&
              forwmono_start->mpNode->mAux == forwmono_end->mpNode->mAux;
              forwmono_end++)
          {}
    
          assert((i <  nforwmononodes-1 && forwmono_end <  p_node->mpLinks+p_node->mNLinks) ||
                (i == nforwmononodes-1 && forwmono_end == p_node->mpLinks+p_node->mNLinks));
    
          link_type *tlink, *backmono_start, *backmono_end = p_node->mpBackLinks;
          
          for (j=0; j < nbackmononodes; j++) 
          {
            for (backmono_start = backmono_end;
              backmono_end < p_node->mpBackLinks+p_node->mNBackLinks &&
              backmono_start->mpNode->mAux == backmono_end->mpNode->mAux;
              backmono_end++)
            {}
    
            assert((j <  nbackmononodes-1 && backmono_end <  p_node->mpBackLinks+p_node->mNBackLinks) ||
                   (j == nbackmononodes-1 && backmono_end == p_node->mpBackLinks+p_node->mNBackLinks));
    
            node_type * tnode;
            if ((tnode = (node_type *) calloc(1, sizeof(node_type))) == NULL)
              Error("Insufficient memory");
            
            *tnode = *p_node;
            tnode->mNLinks       = forwmono_end-forwmono_start;
            tnode->mNBackLinks   = backmono_end-backmono_start;
    
            if ((tnode->mpLinks =
                (link_type *) malloc(tnode->mNLinks * sizeof(link_type))) == NULL ||
              (tnode->mpBackLinks =
                (link_type *) malloc(tnode->mNBackLinks * sizeof(link_type))) == NULL) 
            {
              Error("Insufficient memory");
            }
            
            for (tlink = forwmono_start; tlink < forwmono_end; tlink++) 
            {
              tnode->mpLinks[tlink-forwmono_start] = *tlink;
              tlink->mpNode->mpBackLinks[tlink->mpNode->mNBackLinks  ].mpNode = tnode;
              tlink->mpNode->mpBackLinks[tlink->mpNode->mNBackLinks++].mLmLike = tlink->mLmLike;
              tlink->mpNode->mpBackLinks[tlink->mpNode->mNBackLinks++].mAcousticLike = tlink->mAcousticLike;
            }
            
            for (tlink = backmono_start; tlink < backmono_end; tlink++) 
            {
              tnode->mpBackLinks[tlink-backmono_start] = *tlink;
              tlink->mpNode->mpLinks[tlink->mpNode->mNLinks  ].mpNode = tnode;
              tlink->mpNode->mpLinks[tlink->mpNode->mNLinks++].mLmLike = tlink->mLmLike;
              tlink->mpNode->mpLinks[tlink->mpNode->mNLinks++].mAcousticLike = tlink->mAcousticLike;
            }
            
            prev->mpNext = tnode;
            prev = tnode;
          }
        }
        prev->mpNext = p_node->mpNext;
        free(p_node->mpLinks);
        free(p_node->mpBackLinks);
        free(p_node);
        p_node = prev;
      }
    
      // Give triphone names to phone nodes and create hash of these names
      for (p_node = pFirst(); p_node != NULL; p_node = p_node->mpNext) 
      {
        ENTRY     e       = {0}; //{0} is just to make compiler happy
        ENTRY*    ep      = NULL;
        node_type*     lc      = NULL;
        node_type*     rc      = NULL;
        char*     lcname  = NULL;
        char*     rcname  = NULL;
        char*     triname = NULL;
        int       lcnlen  = 0;
        int       rcnlen  = 0;
    
        if (!(p_node->mType & NT_PHONE)) 
          continue;
    
        if (nonCDphones) 
        {
          e.key  = p_node->mpName;
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
            lc = lc->mNBackLinks ? lc->mpBackLinks[0].mpNode : NULL;
            
            if (lc == NULL)               break;
            if (!(lc->mType & NT_PHONE))  continue;
            if (nonCDphones == NULL)      break;
            
            e.key  = lc->mpName;
            my_hsearch_r(e, FIND, &ep, nonCDphones);
            if (ep == NULL || !reinterpret_cast<size_t>(ep->data)) break; // Node represents Tee model
          }
          
          for (rc = p_node;;) 
          {
            rc = rc->mNLinks ? rc->mpLinks[0].mpNode : NULL;
            
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
        
        triname = (char *) malloc(lcnlen+1+strlen(p_node->mpName)+1+rcnlen+1);
        
        if (triname == NULL) 
          Error("Insufficient memory");
    
        triname[0] = '\0';
    
        if (lcnlen > 0) 
          strcat(strncat(triname, lcname, lcnlen), "-");
          
        strcat(triname, p_node->mpName);
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
          
          p_node->mpName = triname;
        } 
        else 
        {
          free(triname);
          p_node->mpName = ep->key;
        }
      }
    }
  // ExpandMonophoneNetworkToTriphones(MyHSearchData *nonCDphones, ...
  //****************************************************************************

  
  //****************************************************************************
  //****************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    ExpansionsAndOptimizations(
      ExpansionOptions        expOptions,
      const STKNetworkOutputFormat&  rFormat,
      MyHSearchData *         wordHash,
      MyHSearchData *         nonCDphHash,
      MyHSearchData *         triphHash)
    {
      node_type*              p_node(pFirst());

      if (expOptions.mNoWordExpansion  && !expOptions.mCDPhoneExpansion &&
          expOptions.mNoOptimization   && !rFormat.mNoLMLikes &&
          !rFormat.mNoTimes        && !rFormat.mNoWordNodes &&
          !rFormat.mNoModelNodes   && !rFormat.mNoPronunVars) 
      {
        return;
      }
      
      SelfLinksToNullNodes();

      if (!expOptions.mNoWordExpansion) {
        if (!expOptions.mNoOptimization) {
          LatticeLocalOptimization(expOptions.mStrictTiming, expOptions.mTraceFlag);
        }
        assert(wordHash != NULL);
        ExpandByDictionary(wordHash, !expOptions.mRemoveWordsNodes, 
            !expOptions.mRespectPronunVar);
      }

      if (expOptions.mCDPhoneExpansion) {
        if (!expOptions.mNoOptimization) {
          LatticeLocalOptimization(expOptions.mStrictTiming, expOptions.mTraceFlag);
        }
        assert(triphHash != NULL && nonCDphHash != NULL);
        ExpandMonophonesToTriphones(nonCDphHash, triphHash);
      }

      DiscardUnwantedInfo(rFormat);
    
      if (!expOptions.mNoOptimization) {
        LatticeLocalOptimization(expOptions.mStrictTiming, expOptions.mTraceFlag);
      }

      RemoveRedundantNullNodes(p_node);
    } 
  // void NetworkExpansionsAndOptimizations( )
  //****************************************************************************

  
  //****************************************************************************
  //****************************************************************************
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    void 
    Network<_NodeType, _LinkType, _NetworkType>::
    DiscardUnwantedInfo(const STKNetworkOutputFormat& format)
    {
      // TODO : change to iterator
      node_type*  p_node;
      int         i;
    
      for (p_node = pFirst(); p_node != NULL; p_node = p_node->mpNext)  
      {
        if (format.mNoLMLikes) 
        {
          for (i=0; i < p_node->mNLinks;     i++) 
          { 
            p_node->mpLinks    [i].mLmLike = 0.0;
            p_node->mpBackLinks[i].mLmLike = 0.0;
          }
        }

        if (format.mNoAcousticLikes) 
        {
          for (i=0; i < p_node->mNBackLinks; i++) 
          {
            p_node->mpBackLinks[i].mAcousticLike = 0.0;
            p_node->mpLinks    [i].mAcousticLike = 0.0;
          }
        }

        if (format.mNoTimes) {
          p_node->mStop = p_node->mStart = UNDEF_TIME;
        }
        if (format.mNoWordNodes && p_node->mType & NT_WORD) {
          p_node->mpPronun = NULL;
        }
        if (format.mNoModelNodes && (p_node->mType&NT_MODEL || p_node->mType&NT_PHONE)) {
          p_node->mType = NT_WORD;
          p_node->mpPronun = NULL;
        }
        if (format.mNoPronunVars && p_node->mType & NT_WORD && p_node->mpPronun != NULL) {
          p_node->mpPronun = p_node->mpPronun->mpWord->pronuns[0];
        }
      }
    }
    // DiscardUnwantedInfo(STKNetworkOutputFormat format)
    //***************************************************************************
    


} // namespace STK

