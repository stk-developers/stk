
/** @file Net_IO.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

namespace STK
{
  //****************************************************************************
  //****************************************************************************
  template <class _NetworkType>
    void
    ReadSTKNetwork(
      FILE*                     lfp,
      _NetworkType&             rNetwork,
      struct MyHSearchData *    word_hash,
      struct MyHSearchData *    phone_hash,
      int                       notInDict,
      LabelFormat               labelFormat,
      long                      sampPeriod,
      const char *              file_name,
      const char *              in_MLF,
      bool                      compactRepresentation)
    {
      typename _NetworkType::node_type*         p_node;
      const NetworkStorageType  network_storage = _NetworkType::StorageType;

      Node*            node;
      Node*            enode = NULL;
      Node*            first = NULL;
      Node*            last  = NULL;
      Node*            fnode = NULL;
      Node*            lnode;
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
      
      NodeBasic* first_basic;
    
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
          } 
          else if (!*valptr || isspace(*valptr)) 
          { // label definition (field without '=' )
            valptr = chptr;
            chptr = "";
          } 
          else 
          {
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
              if (network_storage == NETWORK_COMPACT)
              {
                 node = reinterpret_cast<Node *>(&first_basic[getInteger(valptr, &chptr, file_name, line_no)]);
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

          if (state == HEADER_DEF) 
          { // label definition
            if (!strcmp(chptr, "S") || !strcmp(chptr, "SUBLAT")) {
              Error("%s not supported (%s:%d)", chptr, file_name, line_no);
            } 
            else if (!strcmp(chptr, "N") || !strcmp(chptr, "NODES")) 
            {
              nnodes  = getInteger(valptr, &chptr, file_name, line_no);

              if (network_storage == NETWORK_COMPACT)
              {
                if (first != NULL)
                  Error("Redefinition of N= (NODES=) is not allowed in CSTK format (%s:%d)", file_name, line_no);
                  
                first = reinterpret_cast<Node *>(first_basic = new NodeBasic[nnodes]);
              }
              else if (node_hash.mTabSize == 0 && !my_hcreate_r(nnodes, &node_hash)) 
              {
                Error("Insufficient memory");
              }
            } 
            else 
            { // Skip unknown header term
              if (getHTKstr(valptr, &chptr)) 
              {
                Error("%s (%s:%d)", chptr, file_name, line_no);
              }
            }
            continue;
          }

          if (state == NODE_DEF) 
          {
            if ((!strcmp(chptr, "time") || !strcmp(chptr, "t")) 
            && !(labelFormat.TIMES_OFF) && !(network_storage == NETWORK_COMPACT)) 
            {
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
            } else if (!strcmp(chptr, "p") && !(network_storage == NETWORK_COMPACT)) {
              node->mPhoneAccuracy = getFloat(valptr, &chptr, file_name, line_no);
            } else if (!strcmp(chptr, "flag") || !strcmp(chptr, "f")) {
              if (getHTKstr(valptr, &chptr)) {
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
                // Count number of link definitions on the rest of the line and prealocate memory for links
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
                
                node->mpLinks = (Link<LINK_BASIC> *) realloc(node->mpLinks, (node->mNLinks + nl) * sizeof(Link<LINK_BASIC>));
                if (node->mpLinks == NULL) Error("Insufficient memory");              
              }
            }
          }
          if (state == ARC_DEF) {
            if (!*chptr || !strcmp(chptr, "END") || !strcmp(chptr, "E")) {
              if (network_storage == NETWORK_COMPACT)
              {
                 enode = reinterpret_cast<Node *>(&first_basic[getInteger(valptr, &chptr, file_name, line_no)]);
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
              node->mpLinks[nl-1].mpNode = enode;
              node->mpLinks[nl-1].mLmLike = 0.0;
              node->mpLinks[nl-1].mAcousticLike = 0.0;
    
              if (!(network_storage == NETWORK_COMPACT))
                ++enode->mNBackLinks;  
            } 
            
            else if (!strcmp(chptr, "language") || !strcmp(chptr, "l")) 
            {
              FLOAT lm_like = getFloat(valptr, &chptr, file_name, line_no);
              
              // Set LM score to link starting in node. This link can possibly
              // lead to a phone node already inserted (div=) between 'node' and'enode'
              node->mpLinks[node->mNLinks-1].mLmLike = lm_like;
            } 
            
            else if (!strcmp(chptr, "acoustic") || !strcmp(chptr, "a")) 
            {
              // Set acoustic score to link starting in node. This link can possibly
              // lead to a phone node already inserted (div=) between 'node' and'enode'
              FLOAT acoustic_like = getFloat(valptr, &chptr, file_name, line_no);
              node->mpLinks[node->mNLinks-1].mAcousticLike = acoustic_like;
            } 

            else if (!strcmp(chptr, "div") || !strcmp(chptr, "d")) 
            {
              
              if (network_storage == NETWORK_COMPACT)
                Error("d= or div= is not allowed in CSTK format (%s:%d)", file_name, line_no);
              
              ENTRY e, *ep;
              char  name[1024];
              float time;
              int   n;
              Node*  last = node;
              FLOAT lm_like  = node->mpLinks[node->mNLinks-1].mLmLike;
    
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
    
                if ((tnode            = (Node *) calloc(1, sizeof(Node))) == NULL
                || (tnode->mpLinks     = (Link<LINK_BASIC> *) malloc(sizeof(Link<LINK_BASIC>))) == NULL
                ) {
                  Error("Insufficient memory");
                }
    
                //Use special type to mark nodes inserted by d=..., they will need
                //special treatment. Later, they will become ordinary NT_PHONE nodes
                tnode->mType      = NT_PHONE | NT_MODEL;
                tnode->mNLinks    = tnode->mNBackLinks = 1;
                tnode->mpBackNext  = enode->mpBackNext;
                enode->mpBackNext  = tnode;
                tnode->mPhoneAccuracy = 1.0;
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
                last->mpLinks[last->mNLinks-1].mLmLike = 0.0;
                last->mpLinks[last->mNLinks-1].mAcousticLike = 0.0;
                last = tnode;
              }
              if (strcmp(phn_marks,":")) {
                Error("Invalid specification of phone marks (d=) (%s:%d)",
                      file_name, line_no);
              }
              last->mpLinks[last->mNLinks-1].mpNode = enode;
              last->mpLinks[last->mNLinks-1].mLmLike = lm_like;
            } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }
          }
        }
      }
      
      if (network_storage == NETWORK_COMPACT)
      {
        if (nnodes == 0)
          Error("No node defined in the network file (%s)", file_name);
          
        for(i = 0; i < nnodes; i++)
        {
          if (first_basic[i].mType == NT_UNDEF)
            Error("Node %d not defined in network file (%s)", i, file_name);
          
          for (int j = 0; j < first_basic[i].mNLinks; j++)
            if (first_basic[i].mpLinks[j].mpNode == first_basic)
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
        if (lnode) lnode->mpNext = NULL;
      
        for (node = lnode; node != NULL; node = node->mpBackNext)
        {
          node->mpBackLinks = (Link<LINK_BASIC> *) malloc(node->mNBackLinks * sizeof(Link<LINK_BASIC>));
          if (node->mpBackLinks == NULL) Error("Insufficient memory");
          
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
          for (i = 0; i < node->mNLinks; i++) {
            Node *forwnode = node->mpLinks[i].mpNode;
            forwnode->mpBackLinks[forwnode->mNBackLinks].mpNode = node;
            forwnode->mpBackLinks[forwnode->mNBackLinks++].mLmLike = node->mpLinks[i].mLmLike;
          }
        }
      
        for (node = lnode; node != NULL; fnode = node, node = node->mpBackNext) {
          //If only stop time is specified, set start time to lowest predecessor stop time
          if (node->mStart == UNDEF_TIME && node->mStop != UNDEF_TIME) {
            int i;
            for (i = 0; i < node->mNBackLinks; i++) {
              Node *backnode = node->mpBackLinks[i].mpNode;
              // skip nodes inserted by d=...
              while (backnode->mType == (NT_PHONE | NT_MODEL)) {
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
          //For model nodes defined by d=... (NT_PHONE | NT_MODEL), node->mStart contains
          //only phone durations. Absolute times must be computed derived starting from
          //the end time of the node to which arc with d=... definition points.
          if (node->mType == (NT_PHONE | NT_MODEL)) {
            assert(node->mNLinks == 1);
            node->mStop = node->mpLinks[0].mpNode->mType == (NT_PHONE | NT_MODEL)
                        && node->mStart != UNDEF_TIME
                        ? node->mpLinks[0].mpNode->mStart : node->mpLinks[0].mpNode->mStop;
            node->mStart = node->mStart != UNDEF_TIME && node->mStop != UNDEF_TIME
                          ? node->mStop - node->mStart : node->mpLinks[0].mpNode->mStart;
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
          if (node->mStart != UNDEF_TIME) {
            node->mStart = (node->mStart - labelFormat.left_extent) / sampPeriod;
          }
          if (node->mStop  != UNDEF_TIME) {
            node->mStop  = (node->mStop + labelFormat.right_extent) / sampPeriod;
          }
        }

        if (first->mpPronun != NULL) 
        {
          node = (Node *) calloc(1, sizeof(Node));
          if (node == NULL) Error("Insufficient memory");
          node->mpNext       = first;
          node->mpBackNext   = NULL;
          first->mpBackNext  = node;
          node->mType       = NT_WORD;
          node->mpPronun     = NULL;
          node->mStart      = UNDEF_TIME;
          node->mStop       = UNDEF_TIME;
          node->mNBackLinks = 0;
          node->mpBackLinks  = NULL;
          node->mNLinks     = 1;
          node->mpLinks      = (Link<LINK_BASIC>*) malloc(sizeof(Link<LINK_BASIC>));
          if (node->mpLinks == NULL) Error("Insufficient memory");
          node->mpLinks[0].mLmLike = 0.0;
          node->mpLinks[0].mAcousticLike = 0.0;
          node->mpLinks[0].mpNode = first;
          first->mNBackLinks = 1;
          first->mpBackLinks  = (Link<LINK_BASIC>*) malloc(sizeof(Link<LINK_BASIC>));
          if (first->mpBackLinks == NULL) Error("Insufficient memory");
          first->mpBackLinks[0].mLmLike = 0.0;
          first->mpBackLinks[0].mAcousticLike = 0.0;
          first->mpBackLinks[0].mpNode = node;
          first = node;
        }

        if (last->mpPronun != NULL) 
        {
          node = (Node *) calloc(1, sizeof(Node));
          if (node == NULL) Error("Insufficient memory");
          last->mpNext      = node;
          node->mpNext      = NULL;
          node->mpBackNext  = last;
          node->mType      = NT_WORD;
          node->mpPronun    = NULL;
          node->mStart     = UNDEF_TIME;
          node->mStop      = UNDEF_TIME;
          node->mNLinks    = 0;
          node->mpLinks     = NULL;
          last->mNLinks = 1;
          last->mpLinks  = (Link<LINK_BASIC>*) malloc(sizeof(Link<LINK_BASIC>));
          if (last->mpLinks == NULL) Error("Insufficient memory");
          last->mpLinks[0].mLmLike = 0.0;
          last->mpLinks[0].mAcousticLike = 0.0;
          last->mpLinks[0].mpNode = node;
          node->mNBackLinks = 1;
          node->mpBackLinks  = (Link<LINK_BASIC>*) malloc(sizeof(Link<LINK_BASIC>));
          if (node->mpBackLinks == NULL) Error("Insufficient memory");
          node->mpBackLinks[0].mLmLike = 0.0;
          node->mpBackLinks[0].mAcousticLike = 0.0;
          node->mpBackLinks[0].mpNode = last;
        }
      }

      // return first;
      rNetwork.SetFirst(first);
    }           
    // Node *ReadSTKNetwork(...)
    //**************************************************************************


  //****************************************************************************
  //****************************************************************************
  template <class _NetworkType>
    void 
    WriteSTKNetwork(
      FILE*                   pFp,
      _NetworkType&           rNetwork,
      STKNetworkOutputFormat  format,
      long                    sampPeriod,
      const char*             net_file,
      const char*             out_MNF)
    {
      int                                     n;
      int                                     l=0;
      typename _NetworkType::iterator         p_node;

    
      // use the mAux field to index the nodes
      for (n = 0, p_node = rNetwork.begin(); p_node != rNetwork.end(); p_node++, n++)
      {
        p_node->mAux = n;
        l += p_node->mNLinks;
      }
      
      fprintf(pFp,"N=%d L=%d\n", n, l);

      for (p_node = rNetwork.begin(); p_node != rNetwork.end(); p_node++)
      {
        int j;
    
        if (format.mAllFieldNames) fputs("I=", pFp);
        if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
        else                      fprintf(pFp,"%d", p_node->mAux);
    
        if (!format.mNoTimes && p_node->mStop != UNDEF_TIME) {
          fputs(" t=", pFp);
    
          if (p_node->mStart != UNDEF_TIME && format.mStartTimes) {
            fprintf(pFp,"%g,", p_node->mStart * 1.0e-7 * sampPeriod);
          }
          fprintf(  pFp,"%g",  p_node->mStop  * 1.0e-7 * sampPeriod);
        }

        if (!(p_node->mType & NT_WORD && p_node->mpPronun == NULL)
          || !format.mNoDefaults) 
        {
          putc(' ', pFp);
          putc(p_node->mType & NT_WORD   ? 'W' :
               p_node->mType & NT_SUBNET ? 'S' :
                                         'M', pFp); // NT_MODEL, NT_PHONE
          putc('=', pFp);
          fprintHTKstr(pFp, p_node->mType & NT_MODEL   ? p_node->mpHmm->mpMacro->mpName   :
                            p_node->mType & NT_WORD    ? (!p_node->mpPronun ? "!NULL" :
                                                      p_node->mpPronun->mpWord->mpName) :
                                                      p_node->mpName); // NT_PHONE (NT_SUBNET)
        }

        if (!format.mNoPronunVars && p_node->mType & NT_WORD
        && p_node->mpPronun != NULL && p_node->mpPronun->mpWord->npronuns > 1
        && (p_node->mpPronun->variant_no > 1 || !format.mNoDefaults))
        {
          fprintf(pFp," v=%d", p_node->mpPronun->variant_no);
        }

        if (p_node->mType & NT_TRUE || p_node->mType & NT_STICKY) 
        {
          fputs(" f=", pFp);
          if (p_node->mType & NT_TRUE)   putc('T', pFp);
          if (p_node->mType & NT_STICKY) putc('K', pFp);
        }

        if (p_node->mType & NT_PHONE && p_node->mPhoneAccuracy != 1.0) {
          fprintf(pFp," p="FLOAT_FMT, p_node->mPhoneAccuracy);
        }

        if (!format.mArcDefsToEnd) 
        {
          if (format.mAllFieldNames) fprintf(pFp," J=%d", p_node->mNLinks);
    
          for (j = 0; j < p_node->mNLinks; j ++) 
          {
            putc(' ', pFp);
            if (format.mAllFieldNames) fputs("E=", pFp);
            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].mpNode->mAux);
            else                     fprintf(pFp,"%d", p_node->mpLinks[j].mpNode->mAux);

            // output language probability
            if ((!close_enough(p_node->mpLinks[j].mLmLike, 0.0, 10)) 
            &&  (!format.mNoLMLikes))
            {
              fprintf(pFp," l="FLOAT_FMT, p_node->mpLinks[j].mLmLike);
            }

            // output acoustic probability
            if ((p_node->mpLinks[j].mAcousticLike != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a="FLOAT_FMT, p_node->mpLinks[j].mAcousticLike);
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
        for (p_node = rNetwork.begin(); p_node != rNetwork.end(); p_node++) 
        {
          int j;
    
          for (j = 0; j < p_node->mNLinks; j ++) 
          {
            if (format.mAllFieldNames) 
              fprintf(pFp, format.mArcDefsWithJ ? "J=%d S=" : "I=", l++);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
            else                     fprintf(pFp,"%d", p_node->mAux);
            putc(' ', pFp); // space = ' ';
            if (format.mAllFieldNames) fputs("E=", pFp);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].mpNode->mAux);
            else                      fprintf(pFp,"%d", p_node->mpLinks[j].mpNode->mAux);

            // output language probability
            if ((!close_enough(p_node->mpLinks[j].mLmLike, 0.0, 10)) 
            && (!format.mNoLMLikes))
            {
              fprintf(pFp," l="FLOAT_FMT, p_node->mpLinks[j].mLmLike);
            }

            // output acoustic probability
            if ((p_node->mpLinks[j].mAcousticLike != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a="FLOAT_FMT, p_node->mpLinks[j].mAcousticLike);
            }

            fputs("\n", pFp);

            if (ferror(pFp)) 
              Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
            
          }
        }
      }
    }
    // WriteSTKNetwork(...
    //***************************************************************************

  /*
  //****************************************************************************
  //****************************************************************************
  template <class _NetworkType>
    void 
    WriteSTKNetwork(
      FILE*                   pFp,
      _NetworkType&           rNetwork,
      STKNetworkOutputFormat  format,
      long                    sampPeriod,
      const char*             net_file,
      const char*             out_MNF)
    {
      int                         n;
      int                         l=0;
      typename _NetworkType::node_type*         p_node;

    
      for (n = 0, p_node = rNetwork.pFirst(); p_node != NULL; p_node = p_node->mpNext, n++)  
      {
        p_node->mAux = n;
        l += p_node->mNLinks;
      }
      
      fprintf(pFp,"N=%d L=%d\n", n, l);
      for (p_node = rNetwork.pFirst(); p_node != NULL; p_node = p_node->mpNext)
      {
        int j;
    
        if (format.mAllFieldNames) fputs("I=", pFp);
        if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
        else                      fprintf(pFp,"%d", p_node->mAux);
    
        if (!format.mNoTimes && p_node->mStop != UNDEF_TIME) {
          fputs(" t=", pFp);
    
          if (p_node->mStart != UNDEF_TIME && format.mStartTimes) {
            fprintf(pFp,"%g,", p_node->mStart * 1.0e-7 * sampPeriod);
          }
          fprintf(  pFp,"%g",  p_node->mStop  * 1.0e-7 * sampPeriod);
        }

        if (!(p_node->mType & NT_WORD && p_node->mpPronun == NULL)
          || !format.mNoDefaults) 
        {
          putc(' ', pFp);
          putc(p_node->mType & NT_WORD   ? 'W' :
               p_node->mType & NT_SUBNET ? 'S' :
                                         'M', pFp); // NT_MODEL, NT_PHONE
          putc('=', pFp);
          fprintHTKstr(pFp, p_node->mType & NT_MODEL   ? p_node->mpHmm->mpMacro->mpName   :
                            p_node->mType & NT_WORD    ? (!p_node->mpPronun ? "!NULL" :
                                                      p_node->mpPronun->mpWord->mpName) :
                                                      p_node->mpName); // NT_PHONE (NT_SUBNET)
        }

        if (!format.mNoPronunVars && p_node->mType & NT_WORD
        && p_node->mpPronun != NULL && p_node->mpPronun->mpWord->npronuns > 1
        && (p_node->mpPronun->variant_no > 1 || !format.mNoDefaults))
        {
          fprintf(pFp," v=%d", p_node->mpPronun->variant_no);
        }

        if (p_node->mType & NT_TRUE || p_node->mType & NT_STICKY) 
        {
          fputs(" f=", pFp);
          if (p_node->mType & NT_TRUE)   putc('T', pFp);
          if (p_node->mType & NT_STICKY) putc('K', pFp);
        }

        if (p_node->mType & NT_PHONE && p_node->mPhoneAccuracy != 1.0) {
          fprintf(pFp," p="FLOAT_FMT, p_node->mPhoneAccuracy);
        }

        if (!format.mArcDefsToEnd) 
        {
          if (format.mAllFieldNames) fprintf(pFp," J=%d", p_node->mNLinks);
    
          for (j = 0; j < p_node->mNLinks; j ++) 
          {
            putc(' ', pFp);
            if (format.mAllFieldNames) fputs("E=", pFp);
            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].mpNode->mAux);
            else                     fprintf(pFp,"%d", p_node->mpLinks[j].mpNode->mAux);

            // output language probability
            if ((!close_enough(p_node->mpLinks[j].mLmLike, 0.0, 10)) 
            &&  (!format.mNoLMLikes))
            {
              fprintf(pFp," l="FLOAT_FMT, p_node->mpLinks[j].mLmLike);
            }

            // output acoustic probability
            if ((p_node->mpLinks[j].mAcousticLike != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a="FLOAT_FMT, p_node->mpLinks[j].mAcousticLike);
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
        for (p_node = rNetwork.pFirst(); p_node != NULL; p_node = p_node->mpNext) 
        {
          int j;
    
          for (j = 0; j < p_node->mNLinks; j ++) 
          {
            if (format.mAllFieldNames) 
              fprintf(pFp, format.mArcDefsWithJ ? "J=%d S=" : "I=", l++);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
            else                     fprintf(pFp,"%d", p_node->mAux);
            putc(' ', pFp); // space = ' ';
            if (format.mAllFieldNames) fputs("E=", pFp);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].mpNode->mAux);
            else                      fprintf(pFp,"%d", p_node->mpLinks[j].mpNode->mAux);

            // output language probability
            if ((!close_enough(p_node->mpLinks[j].mLmLike, 0.0, 10)) 
            && (!format.mNoLMLikes))
            {
              fprintf(pFp," l="FLOAT_FMT, p_node->mpLinks[j].mLmLike);
            }

            // output acoustic probability
            if ((p_node->mpLinks[j].mAcousticLike != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a="FLOAT_FMT, p_node->mpLinks[j].mAcousticLike);
            }

            fputs("\n", pFp);

            if (ferror(pFp)) 
              Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
            
          }
        }
      }
    }
    // WriteSTKNetwork(...
    //***************************************************************************
  */

} // namespace STK


