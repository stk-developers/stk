#!/bin/awk -f
/<MIXTURE>/  {$2=counter++}
nmx != 0     {buff=buff$0"\n"}
/<NUMMIXES>/ {nmx=$2; counter=1}
nmx == 0     {print $0}
(/^[ \t]*~s/ || /<STATE>/ || /<TRANSP>/) && nmx {
  if(nmx != counter-1) print NR ": " nmx-counter+1 " mixture(s) missing" > "/dev/stderr"
  print "<NUMMIXES> " counter-1
  printf buff
  buff = ""
  nmx = 0
}
