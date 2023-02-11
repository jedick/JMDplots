#!/bin/sh
# ARAST/sort_helper.sh
# sort call written in bash to use \t separator
# 20180311

sort -T $1 -t \t -k 1,1 -o $2 $3
