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

#include "DecoderNetwork.h"

#include <sstream>


// CODE
//
namespace STK
{
  
  //###########################################################################
  //###########################################################################
  // GENERAL FUNCTIONS
  //###########################################################################
  //###########################################################################
  
  
  
  //###########################################################################
  //###########################################################################
  // OUTPUT FUNCTIONS
  //###########################################################################
  //###########################################################################
  
  //***************************************************************************
  //***************************************************************************
  int fprintBase62(FILE *fp, int v)
  {
    int         i = 0;
    char        str[16];
    const char* tab = 
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

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
    float d;
    Str2Number(str, &d, endPtr);

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
  
  
} // namespace STK

