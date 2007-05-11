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
  // void 
  // FreeNetwork(Node * pNode) 
  // {
  //   Node *  tnode;
  //   
  //   while (pNode) 
  //   {
  //     tnode = pNode->mpNext;
  //     free(pNode->mpLinks);
  //     free(pNode->mpBackLinks);
  //     free(pNode);
  //     pNode = tnode;
  //   }
  // }
  // 

} // namespace STK
