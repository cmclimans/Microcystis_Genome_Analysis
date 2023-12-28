#!/bin/bash


# $1 is positional argument which should be the markers file
expected=$(wc -l $1 | grep -o "[0-9]\+")

# $2 is positional argument which should be the cat files for a marker's extracted sequence
extracted=$(grep -c '>' $2)

if [ "$expected" = "$extracted" ];then
 echo "equal";
else
 echo "not equal"
 exit 1;
fi
