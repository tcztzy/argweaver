#!/bin/bash

arg=$1
if [[ -d $arg ]]; then
    files=`ls $arg/*.log`
else
    if [[ ! -e $arg ]]; then
	echo "no file $arg"
	exit 1
    fi
    files=$arg
fi

grep -B 1 "sample time" $files |
    awk '$0 ~ /resample_arg_regions/ {subtree=1};
         $0 ~ /resample_arg_leaf/ {subtree=0}
         $NF=="s" {if (subtree==1) {
              totalSubtreeSec += $(NF-1); subtreeCount++
            } else {
              totalLeafSec += $(NF-1); leafCount++;
          }}
         $NF=="m" {if (subtree==1) {
              totalSubtreeSec += $(NF-1)*60; subtreeCount++;
            } else {
              totalLeafSec += $(NF-1)*60; leafCount++;}
         }
         END{avgLeaf=totalLeafSec/60/leafCount;
             avgSubtree=totalSubtreeSec/60/subtreeCount;
             print "leaf: "avgLeaf" m;  subtree: "avgSubtree" m;"}'
