import os
import shutil
import glob
from nipype.interfaces.fsl import Merge
from nipype.interfaces import fsl

def subjectinfo(subject_id,getFeedback=True):
    #Get whether scan is a feedback scan or not
    from pandas import read_csv

    SubjInfo = read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/NARSAD_stimulus_JM.csv')
    SubjInfo.set_index('JM_INTERNAL',inplace=True)
    scan1=SubjInfo.loc[subject_id]['SCAN_1_FEEDBACK']
    if scan1:
        feedback=0
        noFeedback=1
    else:
        feedback=1
        noFeedback=0
    if getFeedback:
        return feedback
    if not getFeedback:
        return noFeedback

#Create subject list
CCD_numbers=[15,17,18,21,22,33,34,40,42,52,59,60,61,63,64,66,74,76,83,89,95,98,99]
subject_list=[]
for ccd in CCD_numbers:
    subject_list.append('CCD0%s' % ccd)
#
# secondlevel_folder_names=['noFeedback','Feedback']
#
from nipype.interfaces.fsl import MultipleRegressDesign
# model = MultipleRegressDesign()
# model.inputs.contrasts = [['group mean', 'T',['reg1'],[1]]]
# model.inputs.regressors = dict(reg1=[1]*len(CCD_numbers))
# model.run()
#
# for i in range(1,6):
#     for fb in [0,1]:
#         for t in ['cope', 'varcope']:
#             x=[]
#             for subj in subject_list:
#                 fbLoc=subjectinfo(subj,fb)
#                 fname = '/home/jmuraskin/Projects/CCD/working_v1/feedback_run-%d/feedback/_subject_id_%s/modelestimate/mapflow/_modelestimate0/results/%s%d.nii.gz' % (fbLoc,subj,t,i)
#                 x.append(fname)
#             subjs = len(x)
#             merger = Merge()
#             merger.inputs.in_files = x
#             merger.inputs.dimension = 't'
#             merger.inputs.output_type = 'NIFTI_GZ'
#             merger.run()
#         flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file='/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='design.mat',t_con_file='design.con', run_mode='flame1')
#         flameo.run()
#         foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/' + secondlevel_folder_names[fb] + '/cope' + str(i)
#         os.mkdir(foldername)
#         shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
#         shutil.move('varcope' + str(i) + '_merged.nii.gz',foldername)
#         shutil.move('stats',foldername)



pairedmodel = MultipleRegressDesign()
pairedmodel.inputs.contrasts = [['A>B', 'T',['reg1'],[1 -1]]]
#make paired ttest model
modelX=[0]*2*len(CCD_numbers)
modelXAB=modelX
modelXAB[:len(CCD_numbers)]=1
modelDict=dict(reg1=modelXAB)
for indx,subj in enumerate(CCD_numbers):
    modeltmp=[0]*2*len(CCD_numbers)
    modeltmp[indx]=1
    modeltmp[indx+len(CCD_numbers)]=1
    modelDict['s%d' % indx]= modeltmp
pairedmodel.inputs.regressors = modelDict
pairedmodel.run()



for i in range(1,6):
    for t in ['cope', 'varcope']:
        x=['/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/Feedback/cope' + str(i) + '/' + t + str(i) + '_merged.nii.gz',\
        '/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/noFeedback/cope' + str(i) + '/' +t + str(i) + '_merged.nii.gz']
        merger = Merge()
        merger.inputs.in_files = x
        merger.inputs.dimension = 't'
        merger.inputs.output_type = 'NIFTI_GZ'
        merger.inputs.merged_file='./' + t + str(i)+'_merged.nii.gz'
        merger.run()

    flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file='/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='design.mat',t_con_file='design.con', run_mode='flame1')
    flameo.run()
    foldername='/home/jmuraskin/Projects/CCD/working/feedback/groupAnalysis/paired-Ttest/cope' + str(i)
    os.mkdir(foldername)
    shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
    shutil.move('varcope' + str(i) + '_merged.nii.gz',foldername)
    shutil.move('stats',foldername)

#
#
# from nipype.interfaces import fsl
# import os
# import shutil
#
# for i in range(1,5):
#   os.mkdir('./cope' + str(i))
#   shutil.move('cope' + str(i) + '_merged.nii.gz','./cope' + str(i) + '/cope' + str(i) + '_merged.nii.gz')
#   shutil.move('varcope' + str(i) + '_merged.nii.gz','./cope' + str(i) + '/varcope' + str(i) + '_merged.nii.gz')
#   flameo = fsl.FLAMEO(cope_file='./cope' + str(i) + '/cope'+str(i)+'_merged.nii.gz',var_cope_file='./cope' + str(i) + '/varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file='/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='design.mat',t_con_file='design.con', run_mode='flame1')
#
#   flameo.run()
#   shutil.move('stats','./cope' + str(i))
