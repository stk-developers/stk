
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
      struct MyHSearchData *    word_hash,
      struct MyHSearchData *    phone_hash,
      int                       notInDict,
      LabelFormat               labelFormat,
      long                      sampPeriod,
      const char *              file_name,
      const char *              in_MLF,
      bool                      compactRepresentation,
      _NetworkType&             rNetwork)
    {
      // to make it easier, we define local typedefs
      typedef typename _NetworkType::NodeType    _node_type;
      typedef typename _NetworkType::LinkType    _link_type;
      
      const NetworkStorageType network_storage = _NetworkType::NetworkType;


      _node_type*      p_node;
      _node_type*      enode = NULL;
      _node_type*      first = NULL;
      _node_type*      last  = NULL;
      _node_type*      fnode = NULL;
      _node_type*      lnode;
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
      
      NodeBasic<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* first_basic;
    
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
        p_node = NULL;

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

          if (state == LINE_START && !strcmp(chptr, "J")) 
          {
            getInteger(valptr, &chptr, file_name, line_no);
            state = AFTER_J;
            continue;
          }

          if (state == AFTER_J) 
          {
            if (*chptr && strcmp(chptr,"START") && strcmp(chptr,"S")) {
              Error("Term 'J=' must be followed by term 'S=' (%s:%d)", file_name, line_no);
            }
            state = LINE_START;
            chptr="";
          }
          
          if (state == LINE_START) 
          {
            if (!*chptr || !strcmp(chptr, "I")) {
              if (compactRepresentation)
              {
                p_node = reinterpret_cast<_node_type *>(&first_basic[getInteger(valptr, &chptr, file_name, line_no)]);
                if (p_node->mType == NT_UNDEF)
                  p_node->mType = NT_WORD;
              }
              else
              {
                if (getHTKstr(valptr, &chptr)) {
                  Error("%s (%s:%d)", chptr, file_name, line_no);
                }
                p_node = find_or_create_node(&node_hash, valptr, &last);
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

              if (compactRepresentation)
              {
                if (first != NULL)
                  Error("Redefinition of N= (NODES=) is not allowed in CSTK format (%s:%d)", file_name, line_no);
                  
                first = reinterpret_cast<_node_type *>(first_basic = new NodeBasic<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>[nnodes]);
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
            && !(labelFormat.TIMES_OFF) && !(compactRepresentation)) 
            {
              char *colonptr=valptr;
              while (*colonptr && !isspace(*colonptr) && *colonptr != ',') colonptr++;
    
              if (*colonptr == ',') {
                if (colonptr != valptr) {
                  *colonptr = ' ';
                  p_node->SetStart(100 * (long long) (0.5 + 1e5 *
                                getFloat(valptr, &chptr, file_name, line_no)));
                }
                valptr = colonptr+1;
              }
              p_node->SetStop(100 * (long long) (0.5 + 1e5 *
                          getFloat(valptr, &chptr, file_name, line_no)));
            } else if (!strcmp(chptr, "var") || !strcmp(chptr, "v")) {
              pron_var = getInteger(valptr, &chptr, file_name, line_no);
              if (pron_var < 1) {
                Error("Invalid pronunciation variant (%s:%d)", file_name, line_no);
              }
            } else if (!strcmp(chptr, "p") && !(compactRepresentation)) {
              p_node->mPhoneAccuracy = getFloat(valptr, &chptr, file_name, line_no);
            } else if (!strcmp(chptr, "flag") || !strcmp(chptr, "f")) {
              if (getHTKstr(valptr, &chptr)) {
                Error("%s (%s:%d)", chptr, file_name, line_no);
              }
              for (; *valptr; valptr++) {
                switch (toupper(*valptr)) {
                  case 'K':
                  case 'F':  p_node->mType |= NT_STICKY; break;
                  case 'T':  p_node->mType |= NT_TRUE;   break;
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
              p_node->mType &= ~(NT_MODEL | NT_PHONE);
              p_node->mType |= NT_WORD;
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
              p_node->mpName = (char *) ep->data;
              p_node->mType &= ~NT_WORD;
              p_node->mType |= NT_PHONE;
            } else if (*chptr=='\0' || !strcmp(chptr,"END") || !strcmp(chptr,"E")) {
              state = ARC_DEF;
            } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }          
            if (state == ARC_DEF || *chptr == '\0') {
              // Node definition is over. For NT_WORD, select right pronun according to
              // word and pron_var; and continue with parsing the arc definition below
              if (p_node->mType & NT_WORD && word != NULL) {
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
                p_node->mpPronun = word->pronuns[pron_var-1];
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
                
                p_node->mpLinks = static_cast<_link_type * >
                  (realloc(p_node->mpLinks, (p_node->mNLinks + nl) * sizeof(_link_type)));

                if (p_node->mpLinks == NULL) 
                  Error("Insufficient memory");              

                for (size_t i_count(p_node->mNLinks); i_count < p_node->mNLinks + nl; i_count++)
                {
                  p_node->mpLinks[i_count].Init();
                }
              }
            }
          }

          if (state == ARC_DEF) 
          {
            if (!*chptr || !strcmp(chptr, "END") || !strcmp(chptr, "E")) {
              if (compactRepresentation)
              {
                 enode = reinterpret_cast<_node_type *>(&first_basic[getInteger(valptr, &chptr, file_name, line_no)]);
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
    
              int nl = ++p_node->mNLinks;
              
              // Links are counted and p_node->mpLinks is properly realocated 
              // at the end of the node definition above
              p_node->mpLinks[nl-1].SetNode(enode);
              p_node->mpLinks[nl-1].SetLmLike(0.0);
              p_node->mpLinks[nl-1].SetAcousticLike(0.0);
    
              // TODO: don't make it dependent on whether the network is compact
              // or not but on the fact, that the nodes are capable of storing
              // backlinks
              if (!(compactRepresentation))
                ++enode->mNBackLinks;  
            } 
            
            else if (!strcmp(chptr, "language") || !strcmp(chptr, "l")) 
            {
              FLOAT lm_like = getFloat(valptr, &chptr, file_name, line_no);
              
              // Set LM score to link starting in node. This link can possibly
              // lead to a phone node already inserted (div=) between 'p_node' and'enode'
              p_node->mpLinks[p_node->mNLinks-1].SetLmLike(lm_like);
            } 
            
            else if (!strcmp(chptr, "acoustic") || !strcmp(chptr, "a")) 
            {
              // Set acoustic score to link starting in p_node. This link can possibly
              // lead to a phone node already inserted (div=) between 'p_node' and'enode'
              FLOAT acoustic_like = getFloat(valptr, &chptr, file_name, line_no);
              p_node->mpLinks[p_node->mNLinks-1].SetAcousticLike(acoustic_like);
            } 

            else if (!strcmp(chptr, "div") || !strcmp(chptr, "d")) 
            {
              if (compactRepresentation)
                Error("d= or div= is not allowed in CSTK format (%s:%d)", file_name, line_no);
              
              ENTRY   e;
              ENTRY*  ep;
              char    name[1024];
              float   time;
              int     n;
              _node_type*   last = p_node;
              FLOAT   lm_like  = p_node->mpLinks[p_node->mNLinks-1].LmLike();
    
              if (p_node->mpLinks[p_node->mNLinks-1].pNode() != enode) 
              {
                Error("Redefinition of  (%s:%d)", chptr, file_name, line_no);
              }
              
              if (getHTKstr(phn_marks=valptr, &chptr)) 
              {
                Error("%s (%s:%d)", chptr, file_name, line_no);
              }
              
              time = -FLT_MAX;
              
              while (sscanf(phn_marks, ":%[^,:]%n,%f%n", name, &n, &time, &n) > 0) 
              {
                _node_type* tnode;
                phn_marks+=n;
    
                if ((tnode          = (_node_type *) calloc(1, sizeof(_node_type))) == NULL
                || (tnode->mpLinks  = (_link_type *) malloc(sizeof(_link_type))) == NULL) 
                {
                  Error("Insufficient memory");
                }
    
                //Use special type to mark nodes inserted by d=..., they will need
                //special treatment. Later, they will become ordinary NT_PHONE nodes
                tnode->mType      = NT_PHONE | NT_MODEL;
                tnode->mNLinks    = tnode->mNBackLinks = 1;
                tnode->mpLinks[0].Init();
                tnode->mpBackNext  = enode->mpBackNext;
                enode->mpBackNext  = tnode;
                tnode->mPhoneAccuracy = 1.0;
                //Store phone durations now. Will be replaced by absolute times below.
                tnode->SetStart(time != -FLT_MAX ? 100 * (long long) (0.5 + 1e5 * time) : UNDEF_TIME);
                tnode->SetStop(UNDEF_TIME);
                e.key  = name;
                e.data = NULL;
                my_hsearch_r(e, FIND, &ep, phone_hash);
    
                if (ep == NULL) 
                {
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
            } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
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
        if (lnode) lnode->mpNext = NULL;
      
        for (p_node = lnode; p_node != NULL; p_node = p_node->mpBackNext)
        {
          p_node->mpBackLinks = (_link_type *) malloc(p_node->mNBackLinks * sizeof(_link_type));
          if (p_node->mpBackLinks == NULL) Error("Insufficient memory");
          
          for (size_t i_count(0); i_count < p_node->mNBackLinks; i_count++)
          { 
            p_node->mpBackLinks[i_count].Init(); 
          }

          if (p_node->mpBackNext) 
            p_node->mpBackNext->mpNext = p_node;
      
          if (p_node->mNLinks == 0) 
          {
            if (last)
              Error("Network has multiple nodes with no successors (%s)", file_name);
            last = p_node;
          }
          
          if (p_node->mNBackLinks == 0) 
          {
            if (first)
              Error("Network has multiple nodes with no predecessor (%s)", file_name);
            first = p_node;
          }
        
          p_node->mNBackLinks = 0;
        }
        if (!first || !last) {
          Error("Network contain no start node or no final node (%s)", file_name);
        }

        for (p_node = lnode; p_node != NULL; p_node = p_node->mpBackNext)
        {
          int i;
          for (i = 0; i < p_node->mNLinks; i++) {
            _node_type* p_forwnode = p_node->mpLinks[i].pNode();
            p_forwnode->mpBackLinks[p_forwnode->mNBackLinks].SetNode(p_node);
            p_forwnode->mpBackLinks[p_forwnode->mNBackLinks++].SetLmLike(p_node->mpLinks[i].LmLike());
          }
        }
      
        for (p_node = lnode; p_node != NULL; fnode = p_node, p_node = p_node->mpBackNext) 
        {
          //If only stop time is specified, set start time to lowest predecessor stop time
          if (p_node->Start() == UNDEF_TIME && p_node->Stop() != UNDEF_TIME) 
          {
            int i;
            for (i = 0; i < p_node->mNBackLinks; i++) {
              _node_type *backnode = p_node->mpBackLinks[i].pNode();
              // skip nodes inserted by d=...
              while (backnode->mType == (NT_PHONE | NT_MODEL)) {
                assert(backnode->mNBackLinks == 1);
                backnode = backnode->mpBackLinks[0].pNode();
              }
              if (backnode->Stop() != UNDEF_TIME) {
                p_node->SetStart(p_node->Start() == UNDEF_TIME
                              ?          backnode->Stop()
                              : LOWER_OF(backnode->Stop(), p_node->Start()));
              }
            }

            if (p_node->Start() == UNDEF_TIME) 
              p_node->SetStart(0);
          }
          //For model nodes defined by d=... (NT_PHONE | NT_MODEL), p_node->Start() contains
          //only phone durations. Absolute times must be computed derived starting from
          //the end time of the node to which arc with d=... definition points.
          if (p_node->mType == (NT_PHONE | NT_MODEL)) {
            assert(p_node->mNLinks == 1);
            p_node->SetStop(p_node->mpLinks[0].pNode()->mType == (NT_PHONE | NT_MODEL)
                        && p_node->Start() != UNDEF_TIME
                        ? p_node->mpLinks[0].pNode()->Start() : p_node->mpLinks[0].pNode()->Stop());
            p_node->SetStart(p_node->Start() != UNDEF_TIME && p_node->Stop() != UNDEF_TIME
                          ? p_node->Stop() - p_node->Start() : p_node->mpLinks[0].pNode()->Start());
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

        if (last != lnode) 
        {
          if (last->mpNext)     last->mpNext->mpBackNext = last->mpBackNext;
          if (last->mpBackNext) last->mpBackNext->mpNext = last->mpNext;
          last->mpNext = NULL;
          lnode->mpNext = last;
          last->mpBackNext = lnode;
        }
        
        for (p_node = first; p_node != NULL; p_node = p_node->mpNext) 
        {
          if (p_node->mType == (NT_PHONE | NT_MODEL)) 
          {
            p_node->mType = NT_PHONE;
          }
          
          if (p_node->Start() != UNDEF_TIME) 
          {
            p_node->SetStart((p_node->Start() - labelFormat.left_extent) / sampPeriod);
          }
          
          if (p_node->Stop()  != UNDEF_TIME) 
          {
            p_node->SetStop((p_node->Stop() + labelFormat.right_extent) / sampPeriod);
          }
        }

        if (first->mpPronun != NULL) 
        {
          p_node = (_node_type *) calloc(1, sizeof(_node_type));
          if (p_node == NULL) Error("Insufficient memory");
          p_node->mpNext       = first;
          p_node->mpBackNext   = NULL;
          first->mpBackNext  = p_node;
          p_node->mType       = NT_WORD;
          p_node->mpPronun     = NULL;
          p_node->SetStart(UNDEF_TIME);
          p_node->SetStop(UNDEF_TIME);
          p_node->mNBackLinks = 0;
          p_node->mpBackLinks  = NULL;
          p_node->mNLinks     = 1;
          p_node->mpLinks      = (_link_type*) malloc(sizeof(_link_type));
          if (p_node->mpLinks == NULL) Error("Insufficient memory");
          p_node->mpLinks[0].Init();
          p_node->mpLinks[0].SetNode(first);
          first->mNBackLinks = 1;
          first->mpBackLinks  = (_link_type*) malloc(sizeof(_link_type));
          if (first->mpBackLinks == NULL) Error("Insufficient memory");
          first->mpBackLinks[0].Init();
          first->mpBackLinks[0].SetNode(p_node);
          first = p_node;
        }

        if (last->mpPronun != NULL) 
        {
          p_node = (_node_type *) calloc(1, sizeof(_node_type));
          if (p_node == NULL) Error("Insufficient memory");
          last->mpNext      = p_node;
          p_node->mpNext      = NULL;
          p_node->mpBackNext  = last;
          p_node->mType      = NT_WORD;
          p_node->mpPronun    = NULL;
          p_node->SetStart(UNDEF_TIME);
          p_node->SetStop(UNDEF_TIME);
          p_node->mNLinks    = 0;
          p_node->mpLinks     = NULL;
          last->mNLinks = 1;
          last->mpLinks  = (_link_type*) malloc(sizeof(_link_type));
          if (last->mpLinks == NULL) Error("Insufficient memory");
          last->mpLinks[0].Init();
          last->mpLinks[0].SetNode(p_node);
          p_node->mNBackLinks = 1;
          p_node->mpBackLinks  = (_link_type*) malloc(sizeof(_link_type));
          if (p_node->mpBackLinks == NULL) Error("Insufficient memory");
          p_node->mpBackLinks[0].Init();
          p_node->mpBackLinks[0].SetNode(last);
          last = p_node;
        }
      }

      // return first;
      rNetwork.SetFirst(first);
      rNetwork.SetLast(last);
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
      const char*             pNetFile,
      const char*             out_MNF,
      const FLOAT&            wPenalty,
      const FLOAT&            lmScale)
    {
      int                                     n;
      int                                     l=0;
      typename _NetworkType::iterator         p_node;
      float                                   lm_scale(lmScale);

    
      // use the mAux field to index the nodes
      for (n = 0, p_node = rNetwork.begin(); p_node != rNetwork.end(); p_node++, n++)
      {
        p_node->mAux = n;
        l += p_node->mNLinks;
      }
      
      fprintf(pFp,"N=%d L=%d\n", n, l);

      for (p_node = rNetwork.begin(); p_node != rNetwork.end(); p_node++)
      {
        if (p_node->NSuccessors() < 1)
          continue;

        int j;
    
        if (format.mAllFieldNames) fputs("I=", pFp);
        if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
        else                      fprintf(pFp,"%d", p_node->mAux);
    
        if (!format.mNoTimes && p_node->Stop() != UNDEF_TIME) 
        {
          fputs(" t=", pFp);
    
          if (p_node->Start() != UNDEF_TIME && format.mStartTimes) {
            fprintf(pFp,"%g,", p_node->Start() * 1.0e-7 * sampPeriod);
          }
          fprintf(  pFp,"%g",  p_node->Stop()  * 1.0e-7 * sampPeriod);
        }

        if (!(p_node->mType & NT_WORD && p_node->mpPronun == NULL)
          || !format.mNoDefaults) 
        {
          putc(' ', pFp);
          putc(p_node->mType & NT_WORD   ? 'W' :
               p_node->mType & NT_SUBNET ? 'S' :
                                           'M', pFp); // NT_MODEL, NT_PHONE
          putc('=', pFp);
          fprintHTKstr(pFp, p_node->mType & NT_MODEL  ? p_node->mpHmm->mpMacro->mpName   :
                            p_node->mType & NT_WORD   ? (!p_node->mpPronun ? "!NULL" :
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
            if (p_node->mpLinks[j].PointsNowhere())
              continue;

            putc(' ', pFp);
            if (format.mAllFieldNames) fputs("E=", pFp);
            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].pNode()->mAux);
            else                      fprintf(pFp,"%d", p_node->mpLinks[j].pNode()->mAux);

            // output language probability
            if ((!close_enough(p_node->mpLinks[j].LmLike(), 0.0, 10)) 
            &&  (!format.mNoLMLikes))
            {
              fprintf(pFp," l="FLOAT_FMT, p_node->mpLinks[j].LmLike() / lm_scale);
            }

            // output acoustic probability
            if ((p_node->mpLinks[j].AcousticLike() != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a="FLOAT_FMT, p_node->mpLinks[j].AcousticLike());
            }
          }
        }

        fputs("\n", pFp);

        if (ferror(pFp)) {
          Error("Cannot write to output network file %s", out_MNF ? out_MNF : pNetFile);
        }
      }
    
      if (format.mArcDefsToEnd) 
      {
        l = 0;
        for (p_node = rNetwork.begin(); p_node != rNetwork.end(); p_node++) 
        {
          if (p_node->NSuccessors() < 1)
            continue;

          int j;
    
          for (j = 0; j < p_node->mNLinks; j ++) 
          {
            if (p_node->mpLinks[j].PointsNowhere())
              continue;

            if (format.mAllFieldNames) 
              fprintf(pFp, format.mArcDefsWithJ ? "J=%d S=" : "I=", l++);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
            else                     fprintf(pFp,"%d", p_node->mAux);
            putc(' ', pFp); // space = ' ';
            if (format.mAllFieldNames) fputs("E=", pFp);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].pNode()->mAux);
            else                      fprintf(pFp,"%d", p_node->mpLinks[j].pNode()->mAux);

            // output language probability
            if ((!close_enough(p_node->mpLinks[j].LmLike(), 0.0, 10)) 
            && (!format.mNoLMLikes))
            {
              fprintf(pFp," l="FLOAT_FMT, p_node->mpLinks[j].LmLike() / lm_scale);
            }

            // output acoustic probability
            if ((p_node->mpLinks[j].AcousticLike() != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a="FLOAT_FMT, p_node->mpLinks[j].AcousticLike());
            }

            fputs("\n", pFp);

            if (ferror(pFp)) 
              Error("Cannot write to output network file %s", out_MNF ? out_MNF : pNetFile);
            
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
      const char*             pNetFile,
      const char*             out_MNF)
    {
      int                         n;
      int                         l=0;
      typename _NetworkType::NodeType*         p_node;

    
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
    
        if (!format.mNoTimes && p_node->Stop() != UNDEF_TIME) {
          fputs(" t=", pFp);
    
          if (p_node->Start() != UNDEF_TIME && format.mStartTimes) {
            fprintf(pFp,"%g,", p_node->Start() * 1.0e-7 * sampPeriod);
          }
          fprintf(  pFp,"%g",  p_node->Stop()  * 1.0e-7 * sampPeriod);
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
            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].pNode()->mAux);
            else                     fprintf(pFp,"%d", p_node->mpLinks[j].pNode()->mAux);

            // output language probability
            if ((!close_enough(p_node->mpLinks[j].mLmLike, 0.0, 10)) 
            &&  (!format.mNoLMLikes))
            {
              fprintf(pFp," l="FLOAT_FMT, p_node->mpLinks[j].mLmLike);
            }

            // output acoustic probability
            if ((p_node->mpLinks[j].AcousticLike() != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a="FLOAT_FMT, p_node->mpLinks[j].mAcousticLike);
            }
          }
        }

        fputs("\n", pFp);

        if (ferror(pFp)) {
          Error("Cannot write to output network file %s", out_MNF ? out_MNF : pNetFile);
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

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mpLinks[j].pNode()->mAux);
            else                      fprintf(pFp,"%d", p_node->mpLinks[j].pNode()->mAux);

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
              Error("Cannot write to output network file %s", out_MNF ? out_MNF : pNetFile);
            
          }
        }
      }
    }
    // WriteSTKNetwork(...
    //***************************************************************************
  */

} // namespace STK


