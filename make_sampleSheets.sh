#!/bin/bash


for i in /home/groupdirs/epibiocore/oreo/keithhenry/20-006/2*
do

	sample=$(basename $i | sed -E 's/re|-re//')
	project_dir=$(echo $i/2*)
	echo -e "$sample\t$project_dir"

done
