import os
import shutil
import glob
from nipype.interfaces.fsl import Merge
from nipype.interfaces import fsl

for i in range(1,6):
    for t in ['cope', 'varcope']:
        x=['/home/jmuraskin/Projects/CCD/working/feedback/groupAnalysis/Feedback/cope' + str(i) + '/' + t + str(i) + '_merged.nii.gz',\
        '/home/jmuraskin/Projects/CCD/working/feedback/groupAnalysis/noFeedback/cope' + str(i) + '/' +t + str(i) + '_merged.nii.gz']
        merger = Merge()
        merger.inputs.in_files = x
        merger.inputs.dimension = 't'
        merger.inputs.output_type = 'NIFTI_GZ'
        merger.inputs.merged_file='./' + t + str(i)+'_merged.nii.gz'
        merger.run()
    flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='pairedTTest_15subjects.grp',mask_file='/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='pairedTTest_15subjects.mat',t_con_file='pairedTTest_15subjects.con', run_mode='flame1')
    flameo.run()
    foldername='/home/jmuraskin/Projects/CCD/working/feedback/groupAnalysis/paired-Ttest/cope' + str(i)
    os.mkdir(foldername)
    shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
    shutil.move('varcope' + str(i) + '_merged.nii.gz',foldername)
    shutil.move('stats',foldername)
