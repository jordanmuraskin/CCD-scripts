# loop through all of the files that
# end in ANAT.nii.gz that are contained
# in the sub directories of the current
# directories

#First start with CCD data in the back up folder

# backupFolder = /data/Projects/CCD/backup


for anat_file in /data/Projects/CCD/backup/*/*ANAT.nii.gz
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
