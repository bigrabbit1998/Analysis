#!/bin/bash

# A -strval, increments until B is reached
A=0
B=$1
C=$2

#if [$B<0]; then
#B=0
#else 
#B=$(($B - 1))
#fi

echo "Starting at A=$A and ending at B=$B with $C events!"

make hist
for i in $(seq $A $B);
do
	./hist.exe $A $C > log
	echo "Run $A completed"
	A=$(($A + 1))


done