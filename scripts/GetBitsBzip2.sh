#!/bin/bash
N=`ls -la $1 | awk '{print $5;}'`;
echo "Size=$N";
bzip2 $1
C=`ls -la $1.bz2 | awk '{print $5;}'`;
echo "Compressed size=$C";
echo "$C*8/$N" | bc -l
bunzip2 $1.bz2
