#!/bin/bash
#Update the subject list name, depending on the dataset collection!
for i in `cat sublist.txt`; do
      bash mriqc_singleSubj.sh $i;
done