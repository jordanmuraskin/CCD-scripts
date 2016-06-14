# loop through all of the files that
# end in ANAT.nii.gz that are contained
# in the sub directories of the current
# directories

#First start with CCD data in the back up folder

backupFolder= "/data/Projects/CCD/backup"

for anat_file in ${backupFolder}/*/*ANAT.nii.gz
do
    # anat_file contains the full path to
    # an anatomical file, print it out
    # so that we can keep track of progress
    echo "working on ${anat_file}"

    # use the dirname function to get the directory
    # name for anat_file
    subj_dir=$( dirname ${anat_file} )

    subj_dir2=$(dirname ${subj_dir})

    finalDirectory="/home/jmuraskin/Projects/CCD/data/${subj_dir2}/anat"
    # rename the anatomical file to anat.nii.gz
    echo "cp -v ${anat_file} ${finalDirectory}/anat.nii.gz"
done
