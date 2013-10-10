#!/bin/bash
./reconf || exit 1
./configure || exit 1
cd src; make || exit 1
echo Done.
