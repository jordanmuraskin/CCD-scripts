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
import argparse


parser = argparse.ArgumentParser(description='Run Second Level Results for DR Resting State Networks')
parser.add_argument('-rwr', help='Option to run with Randomise',required=False,default=True,type=bool)
parser.add_argument('-n',help='Number of Permutations to Run', required=False,default=10000,type=int)
parser.add_argument('-r1samp', help='Option to run 1 sample t-test',required=False,default=True,type=bool)
parser.add_argument('-rall', help='Option to run all subjects or good motion subjects',required=False,default=True,type=bool)
parser.add_argument('-rsn', help='List of Resting-State netorks to run',nargs='+', type=int,required=False,default=[3])
parser.add_argument('-pheno', help='Phenotype Measure to Run', type=str,required=False,default='V1_CCDRSQ_75')
args = parser.parse_args()


#Decide if running all subjects or just good subjects
runWithRandomise =args.rwr
# runFlame= args.rwf
nperms=args.n
run1Sample=args.r1samp
runAll=args.rall
addScanOrder=False
rsn=args.rsn
# runWithRandomise = True
# nperms=10000
# runPair=True
# run1Sample=True
runPairNew=False
# #Decide if running all subjects or just good subjects
# runAll=True

#if run with performance
runWithPerformance=False


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



#load subject list
motionTest=pd.read_csv('CCD_meanFD.csv',names=['Subject_ID','FB','scanorder','meanFD'])
performance=pd.read_csv('CCD_performance.csv',names=['Subject_ID','FB','scanorder','R'])
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

#load phenotypic data
phenoFile='/home/jmuraskin/Projects/CCD/Pheno/narsad+vt_new.csv'
pheno=pd.read_csv(phenoFile)
pheno=pheno.set_index('participant')
pheno_measure_name=args.pheno
pheno_measure = zscore(pheno.loc[subject_list][pheno_measure_name])

secondlevel_folder_names=['noFeedback','Feedback']

# #create second level folders
# folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis'
# for runType in ['randomise','flame']:
#     foldername=folderbase + '/' + runType + '/paired-Ttest/' +  motionDir
#     if not os.path.exists(foldername):
#         os.mkdir(foldername)
#     foldername= foldername + '/' + pheno_measure_name
#     if not os.path.exists(foldername):
#         os.mkdir(foldername)
#
#     for fb in secondlevel_folder_names:
#         foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir
#         if not os.path.exists(foldername):
#             os.mkdir(foldername)
#         foldername= foldername + '/' + pheno_measure_name
#         if not os.path.exists(foldername):
#             os.mkdir(foldername)





for RSN in rsn:

    #create second level folders
    folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/RSN%d/' % RSN



    if not os.path.exists(folderbase):
        os.mkdir(folderbase)

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


    meanFBFolder=meanFBFolder + pheno_measure_name

    if not os.path.exists(meanFBFolder):
        os.mkdir(meanFBFolder)

    meanNFBFolder=meanNFBFolder + pheno_measure_name

    if not os.path.exists(meanNFBFolder):
        os.mkdir(meanNFBFolder)



    if run1Sample:
        for fb in [0,1]:
            x=[]
            for subj in subject_list:
                fbLoc=subjectinfo(subj,fb)
                fname= '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/dr_tempreg_maps_files_to_standard_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/_fwhm_6/_dr_tempreg_maps_files_smooth_0%d/temp_reg_map_000%d_antswarp_maths.nii.gz' % (subj,fbLoc,RSN,RSN)
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
            # if runWithPerformance:
            perf=zscore(performance[performance.FB==fbNames[fb]][performance.Subject_ID.isin(subject_list)]['R'])
            model.inputs.contrasts = [['group mean', 'T',['R'],[1]],['group neg mean', 'T',['R'],[-1]]]
            model.inputs.regressors = dict(reg1=[1]*len(subject_list),FD=list(meanFD),R=list(pheno_measure))

            model.run()

            if runWithRandomise:
                os.mkdir(fbNames[fb])
                if not os.path.exists(fbNames[fb]):
                    os.mkdir(fbNames[fb])
                os.system('mv ./design.* ./%s' % fbNames[fb])
                randomiseCommand='./randomise_forpython.sh -i %s -o ./%s/fb -d ./%s/design.mat -t ./%s/design.con -e ./%s/design.grp -m %s -T -n %d' % ('DMN_merged_%s.nii.gz' % fbNames[fb],fbNames[fb],fbNames[fb],fbNames[fb],fbNames[fb],'/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',nperms)
                os.system(randomiseCommand)
                shutil.move(fbNames[fb],meanFBFolder + '/' + fbNames[fb] if fb else meanNFBFolder + '/' + fbNames[fb])
                shutil.move('DMN_merged_%s.nii.gz' % fbNames[fb],meanFBFolder + '/DMN_merged_%s.nii.gz' % fbNames[fb] if fb else meanNFBFolder + '/DMN_merged_%s.nii.gz' % fbNames[fb])
