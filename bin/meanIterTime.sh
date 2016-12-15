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

grep "sample time" $files |
    awk '$NF=="s" {totalsec += $(NF-1); count++};
         $NF=="m" {totalsec += $(NF-1)*60; count++};
         END {print totalsec/60/count" m"}'
