
/** @file DecoderNetwork_IO.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#include <sstream>
#include "mymath.h"

#define INIT_NODE_HASH_SIZE 1000

namespace STK
{
  //***************************************************************************
  //***************************************************************************
  template <typename _NodeType>
    _NodeType*
    find_or_create_node(struct MyHSearchData *node_hash, const char *node_id, _NodeType** last)
    // Auxiliary function used by ReadHTKLattice_new. (Optionally initialize
    // uninitialized node_hash) Search for node record at key node_id. If found,
    // pointer to this record is returned, otherwise new node record is allocated
    // with type set to NT_UNDEF and entered to has at key node_id and pointer
    // to this new record is returned.
    {
      _NodeType*   p_node;
      ENTRY   e = {0}; //{0} is just to make compiler happy
      ENTRY*  ep;
    
      if (node_hash->mTabSize == 0 && !my_hcreate_r(INIT_NODE_HASH_SIZE, node_hash))
        Error("Insufficient memory");
      
      e.key = const_cast<char *>(node_id);
      my_hsearch_r(e, FIND, &ep, node_hash);
    
      if (ep != NULL) 
        return (_NodeType*) ep->data;  
    
      p_node = (_NodeType*) calloc(1, sizeof(_NodeType));
      
      if (p_node == NULL) 
        Error("Insufficient memory");
    
      p_node->rNLinks()     = 0;
      p_node->rpLinks()     = NULL;
      p_node->rNBackLinks() = 0;
      p_node->rpBackLinks() = NULL;
      p_node->mC.mType       = NT_WORD;
      p_node->mC.SetStart(UNDEF_TIME);
      p_node->mC.SetStop(UNDEF_TIME);
      p_node->mC.mpPronun    = NULL;
      p_node->mpBackNext  = *last;
      p_node->mC.SetPhoneAccuracy(1.0);
      
      *last  = p_node;
      e.key  = strdup(node_id);
      e.data = p_node;
    
      if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, node_hash))
        Error("Insufficient memory");
      
      return p_node;
    } // find_or_create_node(...)


  //****************************************************************************
  //****************************************************************************
  template <class _NetworkType>
    void
    ReadSTKNetwork(
      FILE*                     lfp,
      struct MyHSearchData *    word_hash,
      struct MyHSearchData *    phone_hash,
      int                       notInDict,
//      LabelFormat               labelFormat,
      STKNetworkOutputFormat    labelFormat,
      long                      sampPeriod,
      const char *              file_name,
      const char *              in_MLF,
      bool                      compactRepresentation,
      _NetworkType&             rNetwork)
    {
      // to make it easier, we define local typedefs
      typedef typename _NetworkType::NodeType    _node_type;
      typedef typename _NetworkType::LinkType    _link_type;


      _node_type*      p_node;
      _node_type*      enode = NULL;
      _node_type*      first = NULL;
      _node_type*      last  = NULL;
      _node_type*      fnode = NULL;
      _node_type*      lnode;
      char*            line;
      int              line_no   =  0;
      int              nnodes    =  0;
      char             null_char = '\0';
      char*            chptr;
      char*            valptr;
      char*            phn_marks = NULL;
      Word*            word     = NULL;
      int              i;
      int              pron_var  = 1;
      enum {LINE_START, AFTER_J, HEADER_DEF, ARC_DEF, NODE_DEF} state;
      MyHSearchData node_hash = {0};
      struct ReadlineData   rld       = {0};
      
      for (;;) 
      {
        do 
        {
          if ((chptr = line = readline(lfp, &rld)) == NULL) 
            break;
            
          if (chptr[0] == '.' && (chptr[1] == '\0' || IsSpace(chptr[1]))) {
            chptr = NULL;
            break;
          }
          line_no++;
          while (IsSpace(*chptr)) chptr++;
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
          else if (!*valptr || IsSpace(*valptr)) 
          { // label definition (field without '=' )
            valptr = chptr;
            chptr = &null_char; // = "" (to make compiler happy)
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
              Error("Term 'J=' must be followed by term 'S=' (%s:%d)", 
                  file_name, line_no);
            }
            state = LINE_START;
            chptr = &null_char; // = "" (to make compiler happy)
          }
          
          if (state == LINE_START) 
          {
            if (!*chptr || !strcmp(chptr, "I")) {
              if (compactRepresentation)
              {
                p_node = (&first[getInteger(valptr, &chptr, file_name, line_no)]);

                //p_node = reinterpret_cast<_node_type *>
                  //(&first_basic[getInteger(valptr, &chptr, file_name, line_no)]);

                if (p_node->mC.mType == NT_UNDEF)
                  p_node->mC.mType = NT_WORD;
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
                  Error("Redefinition of N= (NODES=) is not allowed in "
                        "CSTK format (%s:%d)", 
                        file_name, line_no);

                // TODO: !!! first was allocated by malloc, so elements were not initialized
                // e.g. "if (p_node->mC.mType == NT_UNDEF)" above relied on uninitialized mType.
                // calloc initialize all the allocated memory by zero, however, this shuld not be
                // final solution!!!
                
                if (NULL == (first = (_node_type *) calloc(sizeof(_node_type), nnodes)))
                  Error("Insufficient memory");

                last  = &first[nnodes - 1];
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
            && !(labelFormat.mNoTimes) && !(compactRepresentation)) 
            {
              char *colonptr=valptr;

              while (*colonptr && !IsSpace(*colonptr) && *colonptr != ',')
                colonptr++;
    
              if (*colonptr == ',') {
                if (colonptr != valptr) {
                  *colonptr = ' ';
                  p_node->mC.SetStart(100 * (long long) (0.5 + 1e5 *
                                getFloat(valptr, &chptr, file_name, line_no)));
                }
                valptr = colonptr+1;
              }
              p_node->mC.SetStop(100 * (long long) (0.5 + 1e5 *
                          getFloat(valptr, &chptr, file_name, line_no)));
            } else if (!strcmp(chptr, "var") || !strcmp(chptr, "v")) {
              pron_var = getInteger(valptr, &chptr, file_name, line_no);
              if (pron_var < 1) {
                Error("Invalid pronunciation variant (%s:%d)", file_name, line_no);
              }
            } else if (!strcmp(chptr, "p") && !(compactRepresentation)) {
              p_node->mC.SetPhoneAccuracy(getFloat(valptr, &chptr, file_name, line_no));
            } else if (!strcmp(chptr, "flag") || !strcmp(chptr, "f")) {
              if (getHTKstr(valptr, &chptr)) {
                Error("%s (%s:%d)", chptr, file_name, line_no);
              }
              for (; *valptr; valptr++) {
                switch (toupper(*valptr)) {
                  case 'K':
                  case 'F':  p_node->mC.mType |= NT_STICKY; break;
                  case 'T':  p_node->mC.mType |= NT_TRUE;   break;
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
              p_node->mC.mType &= ~(NT_MODEL | NT_PHONE);
              p_node->mC.mType |= NT_WORD;
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
              p_node->mC.mpName = (char *) ep->data;
              p_node->mC.mType &= ~NT_WORD;
              p_node->mC.mType |= NT_PHONE;
            } else if (*chptr=='\0' || !strcmp(chptr,"END") || !strcmp(chptr,"E")) {
              state = ARC_DEF;
            } else if (getHTKstr(valptr, &chptr)) { // Skip unknown term
              Error("%s (%s:%d)", chptr, file_name, line_no);
            }          
            if (state == ARC_DEF || *chptr == '\0') {
              // Node definition is over. For NT_WORD, select right pronun according to
              // word and pron_var; and continue with parsing the arc definition below
              if (p_node->mC.mType & NT_WORD && word != NULL) {
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
                p_node->mC.mpPronun = word->pronuns[pron_var-1];
              }

              if (state == ARC_DEF) 
              {
                // Count number of link definitions on the rest of the line and
                // prealocate memory for links
                const char *pCh;
                int nl = 1;
                
                if (skipHTKstr(valptr, &pCh)) 
                  Error("%s (%s:%d)", pCh, file_name, line_no);
                
                while(*pCh != '\0') 
                {
                  if (!strncmp("END=", pCh, 4) || !strncmp("E=", pCh, 2)) 
                    nl++;
                
                  while(isalnum(*pCh)) pCh++;
                  
                  if (*pCh == '=')
                  {
                    if (skipHTKstr(pCh, &pCh)) 
                      Error("%s (%s:%d)", pCh, file_name, line_no);      
                      
                    continue;
                  }
                  
                  nl++;

                  if(!IsSpace(*pCh) && *pCh != '\0')
                    Error("Invalid character '%c' (%s:%d, char %d)",
                          *pCh, file_name, line_no, pCh-line+1);
                  
                  while(IsSpace(*pCh)) pCh++;
                }
                
                p_node->rpLinks() = static_cast<_link_type * >
                  (realloc(p_node->rpLinks(), (p_node->rNLinks() + nl) 
                           * sizeof(_link_type)));

                if (p_node->rpLinks() == NULL) 
                  Error("Insufficient memory");              

                for (int i_count(p_node->rNLinks()); i_count < p_node->rNLinks()
                      + nl; i_count++)
                {
                  p_node->rpLinks()[i_count].Init();
                }
              }
            }
          }

          if (state == ARC_DEF) 
          {
            if (!*chptr || !strcmp(chptr, "END") || !strcmp(chptr, "E")) {
              if (compactRepresentation)
              {
                 enode = &first[getInteger(valptr, &chptr, file_name, line_no)];
                 //enode = reinterpret_cast<_node_type *>
                 //  (&first_basic[getInteger(valptr, &chptr, file_name, line_no)]);

                 if (enode->mC.mType == NT_UNDEF)
                   enode->mC.mType = NT_WORD;
              }
              else
              {
                if (getHTKstr(valptr, &chptr)) {
                  Error("%s (%s:%d)", chptr, file_name, line_no);
                }
                enode = find_or_create_node(&node_hash, valptr, &last);
              }
    
              int nl = ++p_node->rNLinks();
              
              // Links are counted and p_node->rpLinks() is properly realocated 
              // at the end of the node definition above
              p_node->rpLinks()[nl-1].SetNode(enode);
              p_node->rpLinks()[nl-1].SetLmLike(0.0);
              p_node->rpLinks()[nl-1].SetAcousticLike(0.0);
    
              // TODO: don't make it dependent on whether the network is compact
              // or not but on the fact, that the nodes are capable of storing
              // backlinks
              if (!(compactRepresentation))
                ++enode->rNBackLinks();  
            } 
            else if (!strcmp(chptr, "language") || !strcmp(chptr, "l")) 
            {              
              // Set LM score to link starting in node. This link can possibly
              // lead to a phone node already inserted (div=) between 'p_node'
              // and'enode'
              FLOAT lm_like = getFloat(valptr, &chptr, file_name, line_no);
	      if (!labelFormat.mNoLMLikes)
                p_node->rpLinks()[p_node->rNLinks()-1].SetLmLike(lm_like);
            } 
            else if (!strcmp(chptr, "acoustic") || !strcmp(chptr, "a")) 
            {
              // Set acoustic score to link starting in p_node. This link can
              // possibly lead to a phone node already inserted (div=) between
              // 'p_node' and'enode'
              FLOAT acoustic_like = getFloat(valptr, &chptr, file_name, line_no);
	      if (!labelFormat.mNoAcousticLikes)
                p_node->rpLinks()[p_node->rNLinks()-1].SetAcousticLike(acoustic_like);
            } 
	    else if (!strcmp(chptr, "div") || !strcmp(chptr, "d")) 
	    {
              if (compactRepresentation)
                Error("d= or div= is not allowed in CSTK format (%s:%d)", 
                    file_name, line_no);
              
              ENTRY   e;
              ENTRY*  ep;
              char    name[1024];
              float   time;
              int     n;
              _node_type*   last = p_node;
              FLOAT   lm_like  = p_node->rpLinks()[p_node->rNLinks()-1].LmLike();
    
              if (p_node->rpLinks()[p_node->rNLinks()-1].pNode() != enode) 
              {
                Error("Redefinition of  (%s:%d)", chptr, file_name, line_no);
              }
              
              if (getHTKstr(phn_marks=valptr, &chptr)) 
              {
                Error("%s (%s:%d)", chptr, file_name, line_no);
              }
              
              time = -FLT_MAX;
              
              while (sscanf(phn_marks, ":%[^,:]%n,%n", name, &n, &n) > 0)
              {
                _node_type* tnode;
                phn_marks+=n;
                Str2Number(phn_marks, &time, &phn_marks);
    
                if ((tnode = (_node_type *) calloc(1, sizeof(_node_type))) == NULL
                || (tnode->rpLinks()  = (_link_type *) malloc(sizeof(_link_type))) 
                     == NULL) 
                {
                  Error("Insufficient memory");
                }
    
                //Use special type to mark nodes inserted by d=..., they will need
                //special treatment. Later, they will become ordinary NT_PHONE nodes
                tnode->mC.mType      = NT_PHONE | NT_MODEL;
                tnode->rNLinks()    = tnode->rNBackLinks() = 1;
                tnode->rpLinks()[0].Init();
                tnode->mpBackNext  = enode->mpBackNext;
                enode->mpBackNext  = tnode;
                tnode->mC.SetPhoneAccuracy(1.0);
                //Store phone durations now. Will be replaced by absolute times below.
                tnode->mC.SetStart(time != -FLT_MAX ? 100 * (long long) (0.5 + 1e5 * time) : UNDEF_TIME);
                tnode->mC.SetStop(UNDEF_TIME);
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
                                                   
                tnode->mC.mpName = (char *) ep->data;
                last->rpLinks()[last->rNLinks()-1].SetNode(tnode);
                last->rpLinks()[last->rNLinks()-1].SetLmLike(0.0);
                last->rpLinks()[last->rNLinks()-1].SetAcousticLike(0.0);
                last = tnode;
              }
              if (strcmp(phn_marks,":")) {
                Error("Invalid specification of phone marks (d=) (%s:%d)",
                      file_name, line_no);
              }
              last->rpLinks()[last->rNLinks()-1].SetNode(enode);
              last->rpLinks()[last->rNLinks()-1].SetLmLike(lm_like);
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
          if (first[i].mC.mType == NT_UNDEF)
            Error("Node %d not defined in network file (%s)", i, file_name);
          
          for (int j = 0; j < first[i].rNLinks(); j++)
            if (first[i].rpLinks()[j].pNode() == first)
              Error("Node 0 must be the initial network node (%s)", file_name);
        }
        
        if (first[nnodes-1].rNLinks() != 0)
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
          p_node->rpBackLinks() = (_link_type *) 
            malloc(p_node->rNBackLinks() * sizeof(_link_type));

          if (p_node->rpBackLinks() == NULL) 
            Error("Insufficient memory");
          
          for (int i_count(0); i_count < p_node->rNBackLinks(); i_count++)
          { 
            p_node->rpBackLinks()[i_count].Init(); 
          }

          if (p_node->mpBackNext) 
            p_node->mpBackNext->mpNext = p_node;
      
          if (p_node->rNLinks() == 0) 
          {
            if (last)
              Error("Network has multiple nodes with no successors (%s)", file_name);
            last = p_node;
          }
          
          if (p_node->rNBackLinks() == 0) 
          {
            if (first)
              Error("Network has multiple nodes with no predecessor (%s)", file_name);
            first = p_node;
          }
        
          p_node->rNBackLinks() = 0;
        }
        if (!first || !last) {
          Error("Network contain no start node or no final node (%s)", file_name);
        }

        for (p_node = lnode; p_node != NULL; p_node = p_node->mpBackNext)
        {
          int i;
          for (i = 0; i < p_node->rNLinks(); i++) {
            _node_type* p_forwnode = p_node->rpLinks()[i].pNode();
            p_forwnode->rpBackLinks()[p_forwnode->rNBackLinks()].SetNode(p_node);
            p_forwnode->rpBackLinks()[p_forwnode->rNBackLinks()].SetLmLike(p_node->rpLinks()[i].LmLike());
            p_forwnode->rpBackLinks()[p_forwnode->rNBackLinks()].SetAcousticLike(p_node->rpLinks()[i].AcousticLike());
            ++p_forwnode->rNBackLinks();
          }
        }
      
        for (p_node = lnode; p_node != NULL; fnode = p_node, p_node = p_node->mpBackNext) 
        {
          //If only stop time is specified, set start time to lowest predecessor stop time
          if (p_node->mC.Start() == UNDEF_TIME && p_node->mC.Stop() != UNDEF_TIME) 
          {
            int i;
            for (i = 0; i < p_node->rNBackLinks(); i++) {
              _node_type *backnode = p_node->rpBackLinks()[i].pNode();
              // skip nodes inserted by d=...
              while (backnode->mC.mType == (NT_PHONE | NT_MODEL)) {
                assert(backnode->rNBackLinks() == 1);
                backnode = backnode->rpBackLinks()[0].pNode();
              }
              if (backnode->mC.Stop() != UNDEF_TIME) {
                p_node->mC.SetStart(p_node->mC.Start() == UNDEF_TIME
                              ?          backnode->mC.Stop()
                              : LOWER_OF(backnode->mC.Stop(), p_node->mC.Start()));
              }
            }

            if (p_node->mC.Start() == UNDEF_TIME) 
              p_node->mC.SetStart(0);
          }
          //For model nodes defined by d=... (NT_PHONE | NT_MODEL), p_node->Start() contains
          //only phone durations. Absolute times must be computed derived starting from
          //the end time of the node to which arc with d=... definition points.
          if (p_node->mC.mType == (NT_PHONE | NT_MODEL)) {
            assert(p_node->rNLinks() == 1);
            p_node->mC.SetStop(p_node->rpLinks()[0].pNode()->mC.mType == (NT_PHONE | NT_MODEL)
                        && p_node->mC.Start() != UNDEF_TIME
                        ? p_node->rpLinks()[0].pNode()->mC.Start() : p_node->rpLinks()[0].pNode()->mC.Stop());
            p_node->mC.SetStart(p_node->mC.Start() != UNDEF_TIME && p_node->mC.Stop() != UNDEF_TIME
                          ? p_node->mC.Stop() - p_node->mC.Start() : p_node->rpLinks()[0].pNode()->mC.Start());
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
          if (p_node->mC.mType == (NT_PHONE | NT_MODEL)) 
          {
            p_node->mC.mType = NT_PHONE;
          }
          
          if (p_node->mC.Start() != UNDEF_TIME) 
          {
            p_node->mC.SetStart((p_node->mC.Start() 
	                       + 100 * (long long) (0.5 + 1e5 * labelFormat.mStartTimeShift))
			       / sampPeriod);
          }
          
          if (p_node->mC.Stop()  != UNDEF_TIME) 
          {
            p_node->mC.SetStop((p_node->mC.Stop() 
	                      + 100 * (long long) (0.5 + 1e5 * labelFormat.mEndTimeShift)) 
			      / sampPeriod);
          }
//printf(">%d %d\n", (int) p_node->mC.Start(), (int) p_node->mC.Stop());
        }

        if (first->mC.mpPronun != NULL) 
        {
          p_node = (_node_type *) calloc(1, sizeof(_node_type));
          if (p_node == NULL) Error("Insufficient memory");
          p_node->mpNext       = first;
          p_node->mpBackNext   = NULL;
          first->mpBackNext  = p_node;
          p_node->mC.mType       = NT_WORD;
          p_node->mC.mpPronun     = NULL;
          p_node->mC.SetStart(UNDEF_TIME);
          p_node->mC.SetStop(UNDEF_TIME);
          p_node->rNBackLinks() = 0;
          p_node->rpBackLinks()  = NULL;
          p_node->rNLinks()     = 1;
          p_node->rpLinks()      = (_link_type*) malloc(sizeof(_link_type));
          if (p_node->rpLinks() == NULL) Error("Insufficient memory");
          p_node->rpLinks()[0].Init();
          p_node->rpLinks()[0].SetNode(first);
          first->rNBackLinks() = 1;
          first->rpBackLinks()  = (_link_type*) malloc(sizeof(_link_type));
          if (first->rpBackLinks() == NULL) Error("Insufficient memory");
          first->rpBackLinks()[0].Init();
          first->rpBackLinks()[0].SetNode(p_node);
          first = p_node;
        }

        if (last->mC.mpPronun != NULL) 
        {
          p_node = (_node_type *) calloc(1, sizeof(_node_type));
          if (p_node == NULL) Error("Insufficient memory");
          last->mpNext      = p_node;
          p_node->mpNext      = NULL;
          p_node->mpBackNext  = last;
          p_node->mC.mType      = NT_WORD;
          p_node->mC.mpPronun    = NULL;
          p_node->mC.SetStart(UNDEF_TIME);
          p_node->mC.SetStop(UNDEF_TIME);
          p_node->rNLinks()    = 0;
          p_node->rpLinks()     = NULL;
          last->rNLinks() = 1;
          last->rpLinks()  = (_link_type*) malloc(sizeof(_link_type));
          if (last->rpLinks() == NULL) Error("Insufficient memory");
          last->rpLinks()[0].Init();
          last->rpLinks()[0].SetNode(p_node);
          p_node->rNBackLinks() = 1;
          p_node->rpBackLinks()  = (_link_type*) malloc(sizeof(_link_type));
          if (p_node->rpBackLinks() == NULL) Error("Insufficient memory");
          p_node->rpBackLinks()[0].Init();
          p_node->rpBackLinks()[0].SetNode(last);
          last = p_node;
        }
      }

      // take over the original list
      rNetwork.splice(first, last);
    }           
    // NodeType *ReadSTKNetwork(...)
    //**************************************************************************



  //****************************************************************************
  //****************************************************************************
  template <class _NetworkType>
    void
    ReadSTKNetworkInOldFormat(
      FILE*                     lfp,
      MyHSearchData*            word_hash,
      MyHSearchData*            phone_hash,
      STKNetworkOutputFormat labelFormat,
//      LabelFormat               labelFormat,
      long                      sampPeriod,
      const char*               file_name,
      const char*               in_MLF,
      _NetworkType&             rNetwork)
    {
      // to make it easier, we define local typedefs
      typedef          _NetworkType              _network_type;
      typedef typename _NetworkType::NodeType    _node_type;
      typedef typename _NetworkType::LinkType    _link_type;
      

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
      _node_type*     node;
      _node_type**    nodes;
    
      std::stringstream ss; 
      
      RemoveCommentLines(lfp);
    
      if (fscanf(lfp," %1023[^0-9]", wordOrModelName) == 1) 
      {
        for ( i =0; i < strlen(wordOrModelName); i++) 
        {
          wordOrModelName[i] = toupper(wordOrModelName[i]);
        }
        while (--i>=0 && (wordOrModelName[i] == '='||IsSpace(wordOrModelName[i]))) 
        {
          wordOrModelName[i] = '\0';
        }
      }

      int t1, t2;    
      if ((strcmp(wordOrModelName, "NUMNODES:") &&
          strcmp(wordOrModelName, "NUMBEROFNODES")) ||        //Obsolete NumerOfNodes
        fscanf(lfp," %d NumberOfArcs=%d", &t1, &t2)<1)
      {       //Obsolete NumerOfArcs
        Error("Syntax error in file %s\nKeyword NumNodes: is missing", file_name);
      }
      numOfNodes = t1;
      
      if ((nodes = (_node_type **) calloc(numOfNodes, sizeof(_node_type *))) == NULL) 
      {
        Error("Insufficient memory");
      }

      for (i=0; i < numOfNodes; i++) 
      {
        if ((nodes[i] = (_node_type *) calloc(1, sizeof(_node_type))) == NULL) {
          Error("Insufficient memory");
        }
        nodes[i]->mC.mType = NT_UNDEF;
      }

      for (i=0; i < numOfNodes; i++) 
      {
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
    
        if (node->mC.mType != NT_UNDEF)
          Error("Redefinition of node %d in file %s", nodeId, file_name);
        
        if (toupper(nodeType) != nodeType) {
          nodeType = toupper(nodeType);
          node->mC.mType = 0;
        } else {
          node->mC.mType = NT_TRUE;
        }
        if (nodeType != 'M' && nodeType != 'W' && nodeType != 'N' &&
          nodeType != 'S' && nodeType != 'K' && nodeType != 'F') {
          Error("Invalid definition of node %d in file %s.\n"
                "Supported values for node type are: M - model, "
                "W - word, N - null, S - subnet, K - keyword, F - filler",
                nodeId, file_name);
        }
        node->mC.SetStart(UNDEF_TIME); 
        node->mC.SetStop (UNDEF_TIME);
        node->mC.SetPhoneAccuracy(1.0);
    //    node->mAux = *totalNumOfNodes;
    //    ++*totalNumOfNodes;
    
        if (nodeType == 'S') 
        {
          FILE *snfp;
          _node_type*    subnetFirst = NULL;
          _network_type  my_net(subnetFirst);
    
    //      --*totalNumOfNodes; // Subnet node doesn't count
    
          if ((snfp = fopen(wordOrModelName, "rt")) == NULL) {
            Error("Cannot open network file: %s", wordOrModelName);
          }

          // TODO: This needs an operation of comosing networks...
          ReadSTKNetworkInOldFormat(snfp, word_hash, phone_hash, labelFormat, 
              sampPeriod, wordOrModelName, NULL, my_net);

          my_net.pFirst()->rNBackLinks() = node->rNBackLinks();
          *node = *(my_net.pFirst());
          free(my_net.pFirst());

          fclose(snfp);
        } 
        else if (nodeType == 'M') 
        {
          ENTRY e, *ep;

          node->mC.mType |= NT_PHONE;
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
          node->mC.mpName = (char *) ep->data;
          fscanf(lfp, " {%lf}", &pronunProb);
          // We are not interested in PhoneAccuracy
        } else {
          node->mC.mType |= NT_WORD;
    
          if (nodeType == 'K' || nodeType == 'F') {
          node->mC.mType |= NT_STICKY;
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
    
            while (IsSpace(ch = fgetc(lfp)));
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
                      "Missing ']' after output symbol definition", nodeId,
                      file_name);
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
            node->mC.mpPronun = word ? word->pronuns[pronunVar] : NULL;
          } else {
            node->mC.mpPronun = NULL;
          }
        }

        if (nodeType != 'S') 
        {
          if (fscanf(lfp, " (%lld %lld)", &start, &stop)==2 
          && !(labelFormat.mNoTimes)) 
          {
            long center_shift = labelFormat.mCentreTimes ? sampPeriod / 2 : 0;
            node->mC.SetStart((start - center_shift 
	                     + 100 * (long long) (0.5 + 1e5 * labelFormat.mStartTimeShift))  
                             / sampPeriod);
            node->mC.SetStop((stop  + center_shift 
	                    + 100 * (long long) (0.5 + 1e5 * labelFormat.mEndTimeShift)) 
                            / sampPeriod);
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
          if ((node->rpLinks() = (_link_type *) malloc(numOfLinks * sizeof(_link_type))) == NULL) {
            Error("Insufficient memory");
          }
        } else {
          if (nodeType == 'M') {
            Error("Invalid definition of node %d in file %s.\n"
                  "Model node must have at least one link", nodeId, file_name);
          }
          node->rpLinks() = NULL;
        }
        node->rNLinks() = numOfLinks;
    
        for (j=0; j < static_cast<size_t>(numOfLinks); j++) {
          if (fscanf(lfp, "%d ", &linkId) != 1) {
            Error("Invalid definition of node %d in file %s.\n"
                  "Link Id is expected in link list", nodeId, file_name);
          }
          if (static_cast<size_t>(linkId) >= numOfNodes) {
            Error("Invalid definition of node %d in file %s.\n"
                  "Link Id is bigger than number of nodes", nodeId, file_name);
          }
          if (fscanf(lfp, "{%lf} ", &linkLike) != 1) {
            linkLike = 0.0;
          }
          node->rpLinks()[j].SetNode(nodes[linkId]);
          node->rpLinks()[j].SetLmLike(linkLike);
          ++nodes[linkId]->rNBackLinks();
        }
      }
      for (i = 1; i < numOfNodes-1; i++) {
        if (nodes[i]->rNLinks() == 0) {
          if (nodes[numOfNodes-1]->rNLinks() == 0) {
            Error("Network contains multiple nodes with no successors (%s)",
                  file_name);
          }
          node = nodes[numOfNodes-1];
          nodes[numOfNodes-1] = nodes[i];
          nodes[i] = node;
        }
        if (nodes[i]->rNBackLinks() == 0) {
          if (nodes[0]->rNBackLinks() == 0) {
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
      if (nodes[0]->rNBackLinks() != 0 || nodes[numOfNodes-1]->rNLinks() != 0) {
        Error("Network contain no start node or no final node (%s)", file_name);
      }
      if (!(nodes[0]           ->mC.mType & NT_WORD) || nodes[0]           ->mC.mpPronun != NULL ||
        !(nodes[numOfNodes-1]->mC.mType & NT_WORD) || nodes[numOfNodes-1]->mC.mpPronun != NULL) {
        Error("Start node and final node must be Null nodes (%s)", file_name);
      }
      for (i = 0; i < numOfNodes-1; i++) {
        nodes[i]->mpNext = nodes[i+1];
      }
    
      // create back links
      for (i = 0; i < numOfNodes; i++) {
        if (!nodes[i]->rpBackLinks()) // Could be allready alocated for subnetwork
          nodes[i]->rpBackLinks() = (_link_type *) malloc(nodes[i]->rNBackLinks() * sizeof(_link_type));
        if (nodes[i]->rpBackLinks() == NULL) Error("Insufficient memory");
        nodes[i]->rNBackLinks() = 0;
      }
      for (i = 0; i < numOfNodes; i++) {
        for (j=0; j < static_cast<size_t>(nodes[i]->rNLinks()); j++) {
          _node_type *forwNode = nodes[i]->rpLinks()[j].pNode();

          forwNode->rpBackLinks()[forwNode->rNBackLinks()].SetNode(nodes[i]);
          forwNode->rpBackLinks()[forwNode->rNBackLinks()].SetLmLike(nodes[i]->rpLinks()[j].LmLike());
          forwNode->rpBackLinks()[forwNode->rNBackLinks()].SetAcousticLike(nodes[i]->rpLinks()[j].AcousticLike());
          forwNode->rNBackLinks()++;
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
    
        if ((chptr[0] != '.' || (chptr[1] != '\0' && !IsSpace(chptr[1])))) {
          Error("Missing '.' at the end of network '%s' in NMF '%s'",
                file_name, in_MLF);
        }
      }

      
    }
    // ReadSTKNetworkInOldFormat(
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
      const FLOAT&            wordPenalty,
      const FLOAT&            modelPenalty,
      const FLOAT&            lmScale,
      const FLOAT             posteriorScale)
    {
      int                                     n;
      int                                     l=0;
      typename _NetworkType::Iterator         p_node;
      typedef typename _NetworkType::NodeType    _node_type;
//      float                                   lm_scale(lmScale);
      FLOAT                                   tot_log_like;
    
      // use the mAux field to index the nodes
      for (n = 0, p_node = rNetwork.Begin(); p_node != rNetwork.End(); ++p_node)
      {
        typename _NetworkType::iterator         tmp_node(p_node);
        ++tmp_node;

//        if (tmp_node != rNetwork.End() && p_node->NSuccessors() < 1)
//        {
//          continue;
//        }

        p_node->mAux = n++;
        l += p_node->NSuccessors();
      }
      
      if (format.mPosteriors) {
	assert(NULL != rNetwork.Begin()->mC.mpAlphaBeta);
//        tot_log_like = HIGHER_OF(rNetwork.End()->mC.mpAlphaBeta->mAlpha, 
//                                 rNetwork.Begin()->mC.mpAlphaBeta->mBeta);
	tot_log_like = rNetwork.Begin()->mC.mpAlphaBeta->mBeta;
      }
      
      fprintf(pFp,"N=%d L=%d\n", n, l);

      for (p_node = rNetwork.Begin(); p_node != rNetwork.End(); ++p_node)
      {
        typename _NetworkType::iterator         tmp_node(p_node);
        ++tmp_node;
        
//        if (tmp_node != rNetwork.End() && p_node->NSuccessors() < 1)
//        {
//          continue;
//        }

        int j;
    
        if (format.mAllFieldNames) fputs("I=", pFp);
        if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
        else                      fprintf(pFp,"%d", p_node->mAux);
    
    
// Derive times from neighbours if missing
        if(!format.mNoTimes && p_node->mC.Stop() == UNDEF_TIME && format.mDeriveTimes) {
	  if (p_node->mC.Stop() == UNDEF_TIME) {
            int i;
	    p_node->mC.SetStart(UNDEF_TIME);
            for (i = 0; i < p_node->rNBackLinks(); i++) {
              _node_type *backnode = p_node->rpBackLinks()[i].pNode();
              if (backnode->mC.Stop() != UNDEF_TIME) {
                p_node->mC.SetStart(p_node->mC.Start() == UNDEF_TIME
                                    ?          backnode->mC.Stop()
                                    : LOWER_OF(backnode->mC.Stop(), p_node->mC.Start()));
              }
            }
            if (p_node->mC.Start() != UNDEF_TIME) {

/*  	      /// Hack jak svina platny neb pro JHU workshop. Tohe musi prijit pryc
	      if(p_node->rNLinks() == 1
	         && p_node->rpLinks()[0].pNode()->mC.mType & NT_WORD
	         && p_node->rpLinks()[0].pNode()->mC.mpPronun
	         && strcmp(p_node->rpLinks()[0].pNode()->mC.mpPronun->mpWord->mpName, "sp")
	         && strcmp(p_node->rpLinks()[0].pNode()->mC.mpPronun->mpWord->mpName, "sil")
		 && p_node->rpLinks()[0].pNode()->mC.Start() == UNDEF_TIME
		 && p_node->rpLinks()[0].pNode()->mC.Stop()  == UNDEF_TIME
	         && p_node->rpLinks()[0].pNode()->rNLinks() == 1
	         && p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.mType & NT_WORD
	         && p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.mpPronun
	         && !strcmp(p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.mpPronun->mpWord->mpName, "sp")
		 && p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->rNLinks() > 0
		 && !strcmp(p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.mpPronun->mpWord->mpName, "sp")
	         ) {
		 
		Pronun* pron = p_node->rpLinks()[0].pNode()->mC.mpPronun;
	        p_node->rpLinks()[0].pNode()->mC.mpPronun = p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.mpPronun;
		p_node->rpLinks()[0].pNode()->mC.mpPronun = pron;
		p_node->rpLinks()[0].pNode()->mC.SetStart(p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.Start());
		p_node->rpLinks()[0].pNode()->mC.SetStop (p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.Stop());
		p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.SetStart(UNDEF_TIME);
		p_node->rpLinks()[0].pNode()->rpLinks()[0].pNode()->mC.SetStop(UNDEF_TIME);         
	      }
*/
	      for (i = 0; i < p_node->rNLinks(); i++) {
                _node_type *forwnode = p_node->rpLinks()[i].pNode();
                if (forwnode->mC.Start() != UNDEF_TIME) {
                  p_node->mC.SetStop(p_node->mC.Stop() == UNDEF_TIME
                                      ?           forwnode->mC.Start()
                                      : HIGHER_OF(forwnode->mC.Start(), p_node->mC.Stop()));
                }
              }
	    }
	    
          }
	}
//

        if (!format.mNoTimes && p_node->mC.Stop() != UNDEF_TIME) {
          fputs(" t=", pFp);
    
          if (p_node->mC.Start() != UNDEF_TIME && format.mStartTimes) {
            WriteNumber(pFp, p_node->mC.Start() * 1.0e-7 * sampPeriod);
            fputc(',', pFp);
          }
          WriteNumber(pFp, p_node->mC.Stop() * 1.0e-7 * sampPeriod);
        }

        if (!(p_node->mC.mType & NT_WORD && p_node->mC.mpPronun == NULL)
          || !format.mNoDefaults) 
        {
          putc(' ', pFp);
          putc(p_node->mC.mType & NT_WORD   ? 'W' :
               p_node->mC.mType & NT_SUBNET ? 'S' :
                                           'M', pFp); // NT_MODEL, NT_PHONE
          putc('=', pFp);
          fprintHTKstr(pFp, p_node->mC.mType & NT_MODEL  ? p_node->mC.mpHmm->mpMacro->mpName   :
                            p_node->mC.mType & NT_WORD   ? (!p_node->mC.mpPronun ? "!NULL" :
                                                      p_node->mC.mpPronun->mpWord->mpName) :
                                                      p_node->mC.mpName); // NT_PHONE (NT_SUBNET)
        }

        if (!format.mNoPronunVars && p_node->mC.mType & NT_WORD
        && p_node->mC.mpPronun != NULL && p_node->mC.mpPronun->mpWord->npronuns > 1
        && (p_node->mC.mpPronun->variant_no > 1 || !format.mNoDefaults)) {
          fprintf(pFp," v=%d", p_node->mC.mpPronun->variant_no); 
        } 
        
        if (p_node->mC.mType & NT_TRUE || p_node->mC.mType & NT_STICKY) { 
          fputs(" f=", pFp);
          if (p_node->mC.mType & NT_TRUE)   putc('T', pFp);
          if (p_node->mC.mType & NT_STICKY) putc('K', pFp);
        }
        
        if (format.mPosteriors) {
          assert(NULL != p_node->mC.mpAlphaBeta);
          fprintf(pFp," P=");
          WriteNumber(pFp, my_exp(p_node->mC.mpAlphaBeta->mAlpha + 
                p_node->mC.mpAlphaBeta->mBeta - tot_log_like));
        }
	
        if (p_node->mC.mType & NT_PHONE && p_node->mC.PhoneAccuracy() != 1.0) {
          fprintf(pFp," p=");
          WriteNumber(pFp, p_node->mC.PhoneAccuracy());
        }

        if (!format.mArcDefsToEnd) 
        {
          if (format.mAllFieldNames) {
            fprintf(pFp," J=%d", p_node->rNLinks());
          }
    
          for (j = 0; j < p_node->rNLinks(); j ++) 
          {
            if (p_node->rpLinks()[j].PointsNowhere())
              continue;

            putc(' ', pFp);
            if (format.mAllFieldNames) fputs("E=", pFp);
            if (format.mBase62Labels) fprintBase62(pFp, p_node->rpLinks()[j].pNode()->mAux);
            else                      fprintf(pFp,"%d", p_node->rpLinks()[j].pNode()->mAux);

            FLOAT lm_score = p_node->rpLinks()[j].LmLike();
	    FLOAT ac_score = p_node->rpLinks()[j].AcousticLike();
	    
	    if (format.mPosteriors) {
	      typename _NetworkType::NodeType* p_end_node = p_node->rpLinks()[j].pNode();
	      FLOAT score = p_node->mC.mpAlphaBeta->mAlpha + p_end_node->mC.mpAlphaBeta->mBeta 
                            + (lm_score * lmScale + ac_score) * posteriorScale;

              if (p_end_node->mC.mType & NT_MODEL) {
                score += modelPenalty * posteriorScale;
	      } else if (p_end_node->mC.mType & NT_WORD &&  p_end_node->mC.mpPronun != NULL) {
                score += wordPenalty * posteriorScale;
	      }
              fprintf(pFp," P=");
              WriteNumber(pFp, my_exp(score - tot_log_like));
            }

// THIS WOULD UN-SCALE THE SCORES, BUT THIS NOW HAPPENS JUST AFTER LATTICE
// GENERATION
//            if      (p_node->rpLinks()[j].pNode()->mC.mType & NT_MODEL)
//              lm_score = (p_node->rpLinks()[j].LmLike() - modelPenalty) / lm_scale; 
//            else if (p_node->rpLinks()[j].pNode()->mC.mType & NT_WORD 
//            &&       p_node->rpLinks()[j].pNode()->mC.mpPronun != NULL)
//              lm_score = (p_node->rpLinks()[j].LmLike() - wordPenalty) / lm_scale; 
//            else   
//              lm_score = p_node->rpLinks()[j].LmLike() / lm_scale; 

            // output language probability
            if ((!close_enough(lm_score, 0.0, 10)) 
            &&  (!format.mNoLMLikes))
            {
              fprintf(pFp," l=");
              WriteNumber(pFp, lm_score);
            }

            // output acoustic probability
            if ((!close_enough(p_node->rpLinks()[j].AcousticLike(), 0.0, 10)) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a=");
              WriteNumber(pFp, p_node->rpLinks()[j].AcousticLike());
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
        for (p_node = rNetwork.Begin(); p_node != rNetwork.End(); p_node++) 
        {
          if (p_node->NSuccessors() < 1) {
            continue;
          }

          for (int j = 0; j < p_node->rNLinks(); j ++) 
          {
            if (p_node->rpLinks()[j].PointsNowhere())
              continue;

            if (format.mAllFieldNames) 
              fprintf(pFp, format.mArcDefsWithJ ? "J=%d S=" : "I=", l++);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
            else                      fprintf(pFp,"%d", p_node->mAux);
            putc(' ', pFp); // space = ' ';
            if (format.mAllFieldNames) fputs("E=", pFp);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->rpLinks()[j].pNode()->mAux);
            else                      fprintf(pFp,"%d", p_node->rpLinks()[j].pNode()->mAux);
	    
	    FLOAT lm_score = p_node->rpLinks()[j].LmLike();
	    FLOAT ac_score = p_node->rpLinks()[j].AcousticLike();
	    
	    if (format.mPosteriors) {
	      typename _NetworkType::NodeType* p_end_node = p_node->rpLinks()[j].pNode();
	      FLOAT score = p_node->mC.mpAlphaBeta->mAlpha + p_end_node->mC.mpAlphaBeta->mBeta 
                            + (lm_score * lmScale + ac_score) * posteriorScale;

              if (p_end_node->mC.mType & NT_MODEL) {
                score += modelPenalty * posteriorScale;
	      } else if (p_end_node->mC.mType & NT_WORD &&  p_end_node->mC.mpPronun != NULL) {
                score += wordPenalty * posteriorScale;
	      }
              fprintf(pFp," P=");
              WriteNumber(pFp, my_exp(score - tot_log_like));
            }

// THIS WOULD UN-SCALE THE SCORES, BUT THIS NOW HAPPENS JUST AFTER LATTICE
// GENERATION
//            if      (p_node->rpLinks()[j].pNode()->mC.mType & NT_MODEL)
//              lm_score = (p_node->rpLinks()[j].LmLike() - modelPenalty) / lm_scale; 
//
//            else if (p_node->rpLinks()[j].pNode()->mC.mType & NT_WORD 
//            &&       p_node->rpLinks()[j].pNode()->mC.mpPronun != NULL)
//              lm_score = (p_node->rpLinks()[j].LmLike() - wordPenalty) / lm_scale; 
//
//            else   
//              lm_score = p_node->rpLinks()[j].LmLike() / lm_scale; 

            // output language probability
            if ((!close_enough(lm_score, 0.0, 10)) 
            && (!format.mNoLMLikes))
            {
              fprintf(pFp," l=");
              WriteNumber(pFp, lm_score);
            }

            // output acoustic probability
            if ((!close_enough(p_node->rpLinks()[j].AcousticLike(), 0.0, 10)) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a=");
              WriteNumber(pFp, ac_score);
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

  //***************************************************************************
  //***************************************************************************
  template <class _NetworkType>
    void 
    WriteSTKNetworkInOldFormat(
      FILE*                   pFp,
      _NetworkType&           rNetwork,
      STKNetworkOutputFormat labelFormat,
//      LabelFormat labelFormat,
      long        sampPeriod,
      const char  *net_file,
      const char  *out_MNF)
    {
      // to make it easier, we define local typedefs
      typedef          _NetworkType              _network_type;
      typedef typename _NetworkType::NodeType    _node_type;
      typedef typename _NetworkType::LinkType    _link_type;


      int   i;

      typename _NetworkType::iterator         p_node;
    
      for (i = 0, p_node = rNetwork.Begin(); p_node != rNetwork.End(); 
           ++p_node, ++i)  
      {
        p_node->mAux = i;
      }

      fprintf(pFp,"NUMNODES: %d\n", i);

      for (i = 0, p_node = rNetwork.Begin(); p_node != rNetwork.End(); 
          ++p_node, ++i) 
      {
        int j;
        int type = p_node->mC.mType & NT_MODEL       ? 'M'  :
                   p_node->mC.mType & NT_PHONE       ? 'M'  :
                   p_node->mC.mType & NT_SUBNET      ? 'S'  :
                   p_node->mC.mType & NT_WORD        ?
                    (p_node->mC.mpPronun == NULL     ?
                       (p_node->mC.mType & NT_STICKY ? 'F'  :
                                                    'N') :
                       (p_node->mC.mType & NT_STICKY ? 'K'  :
                                                    'W')):
                                                    '?';
        if (!(p_node->mC.mType & NT_TRUE)) 
        {
          type = tolower(type);
        }

        fprintf(pFp,"%d\t%c %s",
                i, type,
                p_node->mC.mType & NT_MODEL   ? p_node->mC.mpHmm->mpMacro->mpName :
                p_node->mC.mType & NT_PHONE   ? p_node->mC.mpName :
                p_node->mC.mType & NT_SUBNET  ? p_node->mC.mpName :
                p_node->mC.mType & NT_WORD    ?
                  (p_node->mC.mpPronun == NULL ? "-" :
                                          p_node->mC.mpPronun->mpWord->mpName):
                                          "?");
        if (p_node->mC.mType & NT_WORD && p_node->mC.mpPronun) 
        {
          if (p_node->mC.mpPronun->mpWord->mpName != p_node->mC.mpPronun->outSymbol) 
          {
            fprintf(pFp," [%s]", p_node->mC.mpPronun->outSymbol);
          }

          if (p_node->mC.mpPronun->prob != 0.0 
          ||  p_node->mC.mpPronun->mpWord->npronuns > 1) 
          {
            fprintf(pFp," {%d ", p_node->mC.mpPronun->variant_no);
            WriteNumber(pFp, p_node->mC.mpPronun->prob);
            fputc('}', pFp);
          }
        }

        if (p_node->mC.mType & NT_PHONE && p_node->mC.PhoneAccuracy() != 1.0) 
        {
          fprintf(pFp," {");
          WriteNumber(pFp, p_node->mC.PhoneAccuracy());
          fputc('}', pFp);
        }

        if (!(labelFormat.mNoTimes) &&
          p_node->mC.Start() != UNDEF_TIME && p_node->mC.Stop() != UNDEF_TIME) 
        {
          int ctm = labelFormat.mCentreTimes;
          fprintf   (pFp," (");
          fprintf_ll(pFp, sampPeriod * (2 * p_node->mC.Start() + ctm) / 2 
              + 100 * (long long) (0.5 + 1e5 * labelFormat.mStartTimeShift));
          fprintf   (pFp," ");
          fprintf_ll(pFp, sampPeriod * (2 * p_node->mC.Stop() - ctm)  / 2 
              + 100 * (long long) (0.5 + 1e5 * labelFormat.mEndTimeShift));
          fprintf   (pFp,")");
        }

        fprintf(pFp,"\t%d", p_node->rNLinks());
        for (j = 0; j < p_node->rNLinks(); j ++) 
        {
          fprintf(pFp," %d",p_node->rpLinks()[j].pNode()->mAux);
          if (p_node->rpLinks()[j].LmLike() != 0.0) 
          {
            fprintf(pFp," {");
            WriteNumber(pFp, p_node->rpLinks()[j].LmLike());
            fputc('}', pFp);
          }
        }

        fputs("\n", pFp);

        if (ferror(pFp)) {
          Error("Cannot write to output network file %s", out_MNF ? out_MNF : net_file);
        }
      }
    }


#if 0
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
        l += p_node->rNLinks();
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
            WriteNumber(pFp, p_node->Start() * 1.0e-7 * sampPeriod);
            fputc(',', pFp);
          }
          WriteNumber(pFp, p_node->Stop() * 1.0e-7 * sampPeriod);
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
          fprintf(pFp," p=");
          WriteNumber(pFp, p_node->mPhoneAccuracy);
        }

        if (!format.mArcDefsToEnd) 
        {
          if (format.mAllFieldNames) fprintf(pFp," J=%d", p_node->rNLinks());
    
          for (j = 0; j < p_node->rNLinks(); j ++) 
          {
            putc(' ', pFp);
            if (format.mAllFieldNames) fputs("E=", pFp);
            if (format.mBase62Labels) fprintBase62(pFp, p_node->rpLinks()[j].pNode()->mAux);
            else                     fprintf(pFp,"%d", p_node->rpLinks()[j].pNode()->mAux);

            // output language probability
            if ((!close_enough(p_node->rpLinks()[j].mLmLike, 0.0, 10)) 
            &&  (!format.mNoLMLikes))
            {
              fprintf(pFp," l=");
              WriteNumber(pFp, p_node->rpLinks()[j].mLmLike);
            }

            // output acoustic probability
            if ((p_node->rpLinks()[j].AcousticLike() != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a=");
              WriteNumber(pFp, p_node->rpLinks()[j].mAcousticLike);
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
    
          for (j = 0; j < p_node->rNLinks(); j ++) 
          {
            if (format.mAllFieldNames) 
              fprintf(pFp, format.mArcDefsWithJ ? "J=%d S=" : "I=", l++);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->mAux);
            else                     fprintf(pFp,"%d", p_node->mAux);
            putc(' ', pFp); // space = ' ';
            if (format.mAllFieldNames) fputs("E=", pFp);

            if (format.mBase62Labels) fprintBase62(pFp, p_node->rpLinks()[j].pNode()->mAux);
            else                      fprintf(pFp,"%d", p_node->rpLinks()[j].pNode()->mAux);

            // output language probability
            if ((!close_enough(p_node->rpLinks()[j].mLmLike, 0.0, 10)) 
            && (!format.mNoLMLikes))
            {
              fprintf(pFp," l=");
              WriteNumber(pFp, p_node->rpLinks()[j].mLmLike);
            }

            // output acoustic probability
            if ((p_node->rpLinks()[j].mAcousticLike != 0.0) 
            && !(format.mNoAcousticLikes))
            {
              fprintf(pFp," a=");
              WriteNumber(pFp, p_node->rpLinks()[j].mAcousticLike);
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
#endif

} // namespace STK


