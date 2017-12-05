#!/bin/bash
N=`ls -la $1 | awk '{print $5;}'`;
echo "Size=$N";
7za a $1.7z $1
C=`ls -la $1.7z | awk '{print $5;}'`;
echo "Compressed size=$C";
echo "$C*8/$N" | bc -l
