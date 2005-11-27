//
// C++ Interface: %{MODULE}
//
// Description: 
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef STK_Error_h
#define STK_Error_h

#include <stdexcept>


#define Error(...) _Error_(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define Warning(...) _Warning_(__func__, __FILE__, __LINE__, __VA_ARGS__)


namespace STK
{
  /**
   *  @brief Error throwing function
   */
  void _Error_(const char *func, const char *file, int line, char *msg, ...);
  
  /**
   *  @brief Warning handling function
   */
  void _Warning_(const char *func, const char *file, int line, char *msg, ...);
    
  /**
   *  @brief Warning handlingfunction
   */
  void TraceLog(char *msg, ...);

}; // namespace STK

//#define STK_Error_h
#endif
