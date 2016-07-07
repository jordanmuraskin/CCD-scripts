#!/bin/bash
#
# loop through all of the files that
# end in ANAT.nii.gz that are contained
# in the sub directories of the current
# directories

#First start with CCD data in the back up folder

for directory in "data";
do

  #Copy Anatomy
  for train_file in /data/Projects/CCD/$directory/train/train.nii.gz
  do
      # anat_file contains the full path to
      # an anatomical file, print it out
      # so that we can keep track of progress
      echo "working on ${train_file}"

      # use the dirname function to get the directory
      # name for anat_file
      subj_dir=$( dirname ${train_file} )

      subjBase=$(basename ${subj_dir})

      finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
      # mkdir ${finalDirectory}
      # mkdir ${finalDirectory}/anat
      # rename the anatomical file to anat.nii.gz
      cp -vn ${train_file} ${finalDirectory}/feedback/train.nii.gz
  done

  #Copy Info
  for train_file in /data/Projects/CCD/$directory/train/train_info.txt
  do
      # anat_file contains the full path to
      # an anatomical file, print it out
      # so that we can keep track of progress
      echo "working on ${train_file}"

      # use the dirname function to get the directory
      # name for anat_file
      subj_dir=$( dirname ${train_file} )

      subjBase=$(basename ${subj_dir})

      finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
      # mkdir ${finalDirectory}
      # mkdir ${finalDirectory}/anat
      # rename the anatomical file to anat.nii.gz
      cp -vn ${train_file} ${finalDirectory}/feedback/train_info.txt
  done


done
