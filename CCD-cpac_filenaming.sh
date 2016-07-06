#!/bin/bash
#
# For each CCD directory move all func data into REST folder

#First start with CCD data in the back up folder

subjects = ls -d /home/jmuraskin/Projects/CCD/data/CCD*/

for subj in $subjects
do

mkdir /home/jmuraskin/Projects/CCD/data/${subj}/rest

#Copy Anatomy
for filetype in 1,2
do
      fbFile=/home/jmuraskin/Projects/CCD/data/$subj/feedback/fb_${filetype}*


      echo "working on ${subj} and file ${filetype}"


      finalDirectory=/home/jmuraskin/Projects/CCD/data/${subj}/rest/feedback_${filetype}
      mkdir ${finalDirectory}
      cp -vn ${fbFile} ${finalDirectory}
done
