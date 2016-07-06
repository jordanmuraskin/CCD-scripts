#!/bin/bash
#
# For each CCD directory move all func data into REST folder

#First start with CCD data in the back up folder

#subjects= ls -d /home/jmuraskin/Projects/CCD/data/CCD*/

for subj in /home/jmuraskin/Projects/CCD/data/CCD*/
do

	mkdir ${subj}rest

	#Copy Anatomy
	for filetype in 1 2;
	do
      		fbFile=${subj}feedback/fb_${filetype}*
	

	        echo "working on ${subj} and file ${filetype}"


      		finalDirectory=${subj}rest/feedback_${filetype}
      		mkdir ${finalDirectory}
      		cp -vn ${fbFile} ${finalDirectory}
	done

	#Copy Anatomy
        for filetype in 1 2;
        do
                fbFile=${subj}peer/peer${filetype}*


                echo "working on ${subj} and file ${filetype}"


                finalDirectory=${subj}rest/peer_${filetype}
                mkdir ${finalDirectory}
                cp -vn ${fbFile} ${finalDirectory}
        done

	#Copy Anatomy
        for filetype in 1;
        do
                fbFile=${subj}train/train*


                echo "working on ${subj} and file ${filetype}"


                finalDirectory=${subj}rest/train
                mkdir ${finalDirectory}
                cp -vn ${fbFile} ${finalDirectory}
        done


done
