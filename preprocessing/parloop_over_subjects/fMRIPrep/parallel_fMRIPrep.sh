#!/bin/bash

for i in `cat subjList.txt`; do
      cd code;
      winpty bash fmriprep_Scripted.sh $i;
done