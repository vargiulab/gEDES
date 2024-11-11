#!/bin/bash
for i in FRAMESTOSAVE;
do
    catdcd -o frames_hierarchkmeans/$i.pdb -otype pdb -s initial.pdb -stype pdb -first $i -last $i -pdb meta_clusters.pdb
done
