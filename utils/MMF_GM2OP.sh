#!/bin/bash
tmpfile=`mktemp /tmp/$1.XXXXXX`

awk -v i=0 '
/~o/ { do getline; while(!/~/) }
/<STATE>/ || /~s/ { 
  print $0 >"'$tmpfile'"
  getline;
  if(/~s/) print $0 >"'$tmpfile'"
  else print "<OBSCOEF> " ++i >"'$tmpfile'"
  next;
}
/<NUMMIXES>/ || /<SWEIGHTS>/ || /<STREAM>/ ||/<MIXTURE>/ || /<GCONST>/ { next }
/<MEAN>/     { getline; next }
/<VARIANCE>/ { getline; next }
{ print $0 >"'$tmpfile'" }
END { print "~o <VECSIZE> " i " <PDFOBSVEC>" }
' $1

cat $tmpfile
rm $tmpfile


