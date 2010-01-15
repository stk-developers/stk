#ifndef STRTOD_LC_H
#define STRTOD_LC_H

/* locale-independent version of standard C function strtod() */

double strtod_lc(
  const char *string,   /* A decimal ASCII floating-point number,
                         * optionally preceded by white space.
                         * Must have form "-I.FE-X", where I is the
                         * integer part of the mantissa, F is the
                         * fractional part of the mantissa, and X
                         * is the exponent.  Either of the signs
                         * may be "+", "-", or omitted.  Either I
                         * or F may be omitted, or both.  The decimal
                         * point isn't necessary unless F is present.
                         * The "E" may actually be an "e".  E and X
                         * may both be omitted (but not just one).
                         */
  char **endPtr);       /* If non-NULL, store terminating character's
                         * address here. */

#endif
