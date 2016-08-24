import os
import shutil
import glob
import pandas as pd
import numpy as np
from nipype.interfaces.fsl import Merge
from nipype.interfaces import fsl
from subprocess import call
from nipype.interfaces.fsl import MultipleRegressDesign
from scipy.stats import zscore


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
        return feedback+1
    if not getFeedback:
        return noFeedback+1


#Decide if running all subjects or just good subjects
runAll=False

#load subject list
# motionTest=pd.read_csv('CCD_meanFD.csv',names=['Subject_ID','FB','meanFD'])
motionTest=pd.read_csv('./CCD_meanFD.csv',names=['Subject_ID','FB','scanorder','meanFD'])

fbNames=['NOFEEDBACK','FEEDBACK']

if runAll:
    subject_list=np.unique(motionTest.Subject_ID)
    motionDir='all'
else:
    motionThresh=0.2
    allsubj=np.unique(motionTest['Subject_ID'])
    motionReject=np.unique((motionTest[motionTest.meanFD>motionThresh]['Subject_ID']))
    subject_list=np.setdiff1d(allsubj,motionReject)
    motionDir='motionThresh-%f' % motionThresh


#create second level folders
folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/DMN/'



if not os.path.exists(folderbase):
    os.mkdir(folderbase)

pairedFolder= folderbase + 'scanorder/'

if not os.path.exists(pairedFolder):
    os.mkdir(pairedFolder)

pairedFolder= pairedFolder + motionDir

if not os.path.exists(pairedFolder):
    os.mkdir(pairedFolder)



meanFBFolder=folderbase + 'meanFB'
if not os.path.exists(meanFBFolder):
    os.mkdir(meanFBFolder)

meanNFBFolder=folderbase + 'meanNFB'
if not os.path.exists(meanNFBFolder):
    os.mkdir(meanNFBFolder)

meanFBFolder=meanFBFolder + '/' + motionDir
meanNFBFolder=meanNFBFolder + '/' + motionDir

if not os.path.exists(meanFBFolder):
    os.mkdir(meanFBFolder)

if not os.path.exists(meanNFBFolder):
    os.mkdir(meanNFBFolder)

runWithRandomise = True
nperms=10000
# runPair=True
run1Sample=True

if run1Sample:
    for fb in [0,1]:
        x=[]
        for subj in subject_list:
            fbLoc=subjectinfo(subj,fb)
            fname= '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/dr_tempreg_maps_files_to_standard_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/_fwhm_6/_dr_tempreg_maps_files_smooth_03/temp_reg_map_0003_antswarp_maths.nii.gz' % (subj,fbLoc)
            # fname = '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/dr_tempreg_maps_files_to_standard_smooth/_scan_feedback_%d/%s%d.nii.gz' % (fbLoc,subj,t,i)
            x.append(fname)
        subjs = len(x)
        merger = Merge()
        merger.inputs.in_files = x
        merger.inputs.dimension = 't'
        merger.inputs.output_type = 'NIFTI_GZ'
        merger.inputs.merged_file = './DMN_merged_%s.nii.gz' % fbNames[fb]
        merger.run()
        #get meanFD values for each subject and add as covariate
        meanFD=zscore(motionTest[motionTest.FB==fbNames[fb]][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
        model = MultipleRegressDesign()
        model.inputs.contrasts = [['group mean', 'T',['reg1'],[1]],['group neg mean', 'T',['reg1'],[-1]]]
        model.inputs.regressors = dict(reg1=list(motionTest[motionTest.FB==fbNames[fb]][motionTest.Subject_ID.isin(goodsubj)]['scanorder']-1.5),FD=list(meanFD))
        model.run()

        if runWithRandomise:
            os.mkdir(fbNames[fb])
            randomiseCommand='./randomise_forpython.sh -i %s -o ./%s/fb -d design.mat -t design.con -e design.grp -m %s -T -n %d' % ('DMN_merged_%s.nii.gz' % fbNames[fb],fbNames[fb],'/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',nperms)
            os.system(randomiseCommand)
            shutil.move(fbNames[fb],meanFBFolder + '/' + fbNames[fb] if fb else meanNFBFolder + '/' + fbNames[fb])
            shutil.move('DMN_merged_%s.nii.gz' % fbNames[fb],meanFBFolder + '/DMN_merged_%s.nii.gz' % fbNames[fb] if fb else meanNFBFolder + '/DMN_merged_%s.nii.gz' % fbNames[fb])




# if runPair:
#     pairedmodel = MultipleRegressDesign()
#     pairedmodel.inputs.contrasts = [['S1>S2', 'T',['S1','S2'],[1,-1]],['S2>S1', 'T',['S1','S2'],[-1,1]]]
#     pairedmodel.inputs.groups = [1]*len(subject_list)
#     #make paired ttest model
#     modelX=[0]*2*len(subject_list)
#     modelXAB=modelX
#     modelXAB[0:len(subject_list)]=[1]*len(subject_list)
#     modelDict=dict(reg1=modelXAB)
#     for indx,subj in enumerate(subject_list):
#         modeltmp=[0]*2*len(subject_list)
#         modeltmp[indx]=1
#         modeltmp[indx+len(subject_list)]=1
#         modelDict['s%d' % indx]= modeltmp
#     modelDict['FD'] = list(zscore(list(motionTest[motionTest.FB=='FEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
#     + list(motionTest[motionTest.FB=='NOFEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD'])))
#     pairedmodel.inputs.regressors = modelDict
#     pairedmodel.run()
#
#
#     x=[meanFBFolder + '/DMN_merged_FEEDBACK.nii.gz',\
#     meanNFBFolder + '/DMN_merged_NOFEEDBACK.nii.gz']
#     merger = Merge()
#     merger.inputs.in_files = x
#     merger.inputs.dimension = 't'
#     merger.inputs.output_type = 'NIFTI_GZ'
#     merger.inputs.merged_file='./DMN_pair_merged.nii.gz'
#     merger.run()
#
#     os.mkdir('DMN_pair')
#     randomiseCommand='./randomise_forpython.sh -i %s -o ./DMN_pair/paired -d design.mat -t design.con -e design.grp -m %s -T -n %d' % ('DMN_pair_merged.nii.gz','/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',nperms)
#     os.system(randomiseCommand)
#
#
#     shutil.move('DMN_pair',pairedFolder + '/DMN_pair')
#     shutil.move('DMN_pair_merged.nii.gz',pairedFolder + '/DMN_pair_merged.nii.gz')
