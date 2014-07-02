#!/bin/bash
SPU=experiments/SPUhull.data
# make sure file is ordered by var-cardinality (column 4) then experiment name (column 1)
cat $1 | sort -k 4,1 -g > .tmp
mv .tmp $1
paste $SPU $1 | awk -f stats.awk
