//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "Error.h"
#include <stdarg.h>
#include <stdio.h>


namespace STK
{
  void _Error_(const char *func, const char *file, int line, char *msg, ...) {
    va_list ap;
    va_start(ap, msg);
    fflush(stdout);
    fprintf(stderr, "ERROR (%s:%d) %s: ", file, line, func);
    vfprintf(stderr, msg, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(1);
  }
  
  void _Warning_(const char *func, const char *file, int line, char *msg, ...) {
    va_list ap;
    va_start(ap, msg);
    fflush(stdout);
    fprintf(stderr, "WARNING (%s:%d) %s: ", file, line, func);
    vfprintf(stderr, msg, ap);
    fputc('\n', stderr);
    va_end(ap);
  }
  
  void TraceLog(char *msg, ...) {
    va_list ap;
    va_start(ap, msg);
    fflush(stderr);
    vfprintf(stdout, msg, ap);
    fputc('\n', stdout);
    va_end(ap);
  }

} // namespace STK
