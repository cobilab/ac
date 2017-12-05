#!/bin/bash
N=`ls -la $1 | awk '{print $5;}'`;
echo "Size=$N";
lzma $1
C=`ls -la $1.lzma | awk '{print $5;}'`;
echo "Compressed size=$C";
echo "$C*8/$N" | bc -l
lzma -d $1.lzma
