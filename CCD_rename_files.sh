#!/bin/bash
#
# loop through all of the files that
# end in ANAT.nii.gz that are contained
# in the sub directories of the current
# directories

#First start with CCD data in the back up folder

for directory in [backup, newCCD];
do

#Copy Anatomy
for anat_file in /data/Projects/CCD/${directory}/*/*ANAT.nii.gz
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/anat
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/anat/anat.nii.gz
done

#Copy Anatomy .txt
for anat_file in /data/Projects/CCD/${directory}/*/*ANAT_info.txt
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/anat
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/anat/anat_info.txt
done

#Copy Tracking TRAIN
for anat_file in /data/Projects/CCD/${directory}/*/*TRAIN.nii.gz
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/train
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/train/train.nii.gz
done

#Copy Tracking TRAIN .txt
for anat_file in /data/Projects/CCD/${directory}/*/*TRAIN_info.txt
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/train
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/train/train_info.txt
done

#Copy Peer1
for anat_file in /data/Projects/CCD/${directory}/*/*PEER1.nii.gz
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/peer
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/peer/peer1.nii.gz
done

#Copy peer1 .txt
for anat_file in /data/Projects/CCD/${directory}/*/*PEER1_info.txt
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/peer
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/peer/peer1_info.txt
done

#Copy Peer2
for anat_file in /data/Projects/CCD/${directory}/*/*PEER2.nii.gz
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/peer
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/peer/peer2.nii.gz
done

#Copy peer2 .txt
for anat_file in /data/Projects/CCD/${directory}/*/*PEER2_info.txt
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subjBase=$(basename ${subj_dir})

    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase}
    mkdir ${finalDirectory}
    mkdir ${finalDirectory}/peer
    # rename the anatomical file to anat.nii.gz
    cp -vn ${anat_file} ${finalDirectory}/peer/peer2_info.txt
done


#Copy Feedback Files but make sure there are 2 copied
#Go through each subject

for subj in /data/Projects/CCD/${directory}/*CCD*;
do
  subjBase=$(basename ${subj});
  count=1;
  for file in ${subj}/*TEST.nii.gz;
  do
    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase};
    mkdir ${finalDirectory};
    mkdir ${finalDirectory}/feedback;

    cp -vn ${file} ${finalDirectory}/feedback/fb_${count}.nii.gz;

    count=2
  done

  count=1;
  for file in ${subj}/*TEST_info.txt;
  do
    finalDirectory=/home/jmuraskin/Projects/CCD/data/${subjBase};
    mkdir ${finalDirectory};
    mkdir ${finalDirectory}/feedback;

    cp -vn ${file} ${finalDirectory}/feedback/fb_${count}_info.txt;

    count=2
  done
done
done
