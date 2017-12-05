#!/bin/bash
N=`ls -la $1 | awk '{print $5;}'`;
echo "Size=$N";
gzip $1
C=`ls -la $1.gz | awk '{print $5;}'`;
echo "Compressed size=$C";
echo "$C*8/$N" | bc -l
gunzip $1.gz
