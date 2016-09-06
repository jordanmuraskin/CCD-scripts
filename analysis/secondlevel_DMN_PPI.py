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




parser = argparse.ArgumentParser(description='Run Second Level Results for PPI DMN CCD')
parser.add_argument('-rwr', help='Option to run with Randomise',required=False,default=True,type=bool)
parser.add_argument('-rwf', help='Option to run with FLAME',required=False,default=False,type=bool)
parser.add_argument('-n',help='Number of Permutations to Run', required=False,default=10000,type=int)
parser.add_argument('-r1samp', help='Option to run 1 sample t-test',required=False,default=True,type=bool)
parser.add_argument('-rpair', help='Option to run paired t-test',required=False,default=True,type=bool)
parser.add_argument('-rall', help='Option to run all subjects or good motion subjects',required=False,default=True,type=bool)
parser.add_argument('-copes', help='List of copes to run',nargs='+', type=int,required=False,default=[3])
args = parser.parse_args()

#Decide if running all subjects or just good subjects
runWithRandomise =args.rwr
runFlame= args.rwf
nperms=args.n
runPair=args.rpair
run1Sample=args.r1samp
runAll=args.rall
addScanOrder=False
copesToRun=args.copes

# #Decide if running all subjects or just good subjects
# runAll=False
# runWithRandomise = True
# runFlame=False
# nperms=10000
# runPair=True
# run1Sample=True


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




#load subject list
motionTest=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_meanFD.csv',names=['Subject_ID','FB','scanorder','meanFD'])
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


# #Create subject list
# CCD_numbers=[12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,31,32,33,34,40,42,51,
# 53,59,60,61,62,63,64,65,66,67,71,72,73,74,75,76,80,81,82,83,84,85,86,87,88,89,
# 90,91,93,95,97,98,99]

# subject_list=[]
# for ccd in CCD_numbers:
#     subject_list.append('CCD0%s' % ccd)
#
secondlevel_folder_names=['noFeedback','Feedback']

#create second level folders
folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/PPI_DMN'
if not os.path.exists(folderbase):
    os.mkdir(folderbase)



for runType in ['randomise','flame']:
    foldername=folderbase + '/' + runType
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    foldername=folderbase + '/' + runType + '/paired-Ttest/'
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    foldername=folderbase + '/' + runType + '/paired-Ttest/' +  motionDir
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    for fb in secondlevel_folder_names:
        foldername=folderbase + '/' + runType + '/' + fb
        if not os.path.exists(foldername):
            os.mkdir(foldername)

        foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir
        if not os.path.exists(foldername):
            os.mkdir(foldername)





if run1Sample:


    for i in copesToRun:
        for fb in [0,1]:
            for t in ['cope', 'varcope']:
                x=[]
                for subj in subject_list:
                    fbLoc=subjectinfo(subj,fb)
                    fname = '/home/jmuraskin/Projects/CCD/working_v1/PPI_DMN/feedback_run-%d/PPI_model_fit/_subject_id_%s/modelestimate/mapflow/_modelestimate0/results/%s%d.nii.gz' % (fbLoc,subj,t,i)
                    x.append(fname)
                subjs = len(x)
                merger = Merge()
                merger.inputs.in_files = x
                merger.inputs.dimension = 't'
                merger.inputs.output_type = 'NIFTI_GZ'
                merger.run()
            #get meanFD values for each subject and add as covariate
            meanFD=zscore(motionTest[motionTest.FB==fbNames[fb]][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
            model = MultipleRegressDesign()
            model.inputs.contrasts = [['group mean', 'T',['reg1'],[1]],['group neg mean', 'T',['reg1'],[-1]]]
            model.inputs.regressors = dict(reg1=[1]*len(subject_list),FD=list(meanFD))
            model.run()

            if runFlame:
                flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file='/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='design.mat',t_con_file='design.con', run_mode='flame1')
                flameo.run()
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/PPI_DMN/flame/' + secondlevel_folder_names[fb] + '/' + motionDir + '/cope' + str(i)
                if os.path.exists(foldername):
                    shutil.rmtree(foldername)
                    os.mkdir(foldername)
                else:
                    os.mkdir(foldername)
                if not runWithRandomise:
                    shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
                shutil.move('varcope' + str(i) + '_merged.nii.gz',foldername)
                shutil.move('stats',foldername + '/stats')
            if runWithRandomise:
                if not os.path.exists('cope%d' % i):
                    os.mkdir('cope%d' % i)
                os.system('mv ./design.* ./cope%d' % i)
                # shutil.move('./design.*','cope%d' % i)
                randomiseCommand='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/randomise_forpython.sh -i %s -o ./cope%d/cope%d -d ./cope%d/design.mat -t ./cope%d/design.con -e ./cope%d/design.grp -m %s -T -n %d' % ('cope' + str(i) + '_merged.nii.gz',i,i,i,i,i,'/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',nperms)
                os.system(randomiseCommand)

                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/PPI_DMN/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/cope' + str(i)

                if os.path.exists(foldername):
                    shutil.rmtree(foldername)
                    os.mkdir(foldername)
                else:
                    os.mkdir(foldername)
                filename='cope%d' % i
                shutil.move(filename, os.path.join(foldername, filename))
                # shutil.move('cope%d' % i,foldername)
                filename='cope' + str(i) + '_merged.nii.gz'
                # shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
                shutil.move(filename, os.path.join(foldername, filename))




if runPair:


    for i in copesToRun:

        pairedmodel = MultipleRegressDesign()
        pairedmodel.inputs.contrasts = [['A>B', 'T',['reg1'],[1]],['B>A', 'T',['reg1'],[-1]]]
        if runFlame:
            pairedmodel.inputs.groups = [1]*len(subject_list)*2
        else:
            pairedmodel.inputs.groups = range(1,len(subject_list)+1) + range(1,len(subject_list)+1)
        #make paired ttest model
        modelX=[0]*2*len(subject_list)
        modelXAB=modelX
        modelXAB[0:len(subject_list)]=[1]*len(subject_list)
        modelDict=dict(reg1=modelXAB)
        for indx,subj in enumerate(subject_list):
            modeltmp=[0]*2*len(subject_list)
            modeltmp[indx]=1
            modeltmp[indx+len(subject_list)]=1
            modelDict['s%d' % indx]= modeltmp
        modelDict['FD'] = list(zscore(list(motionTest[motionTest.FB=='FEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
        + list(motionTest[motionTest.FB=='NOFEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD'])))
        pairedmodel.inputs.regressors = modelDict
        pairedmodel.run()



        for t in ['cope', 'varcope']:
            try:
                x=['/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/PPI_DMN/randomise/Feedback/' + motionDir +'/cope' + str(i) + '/' + t + str(i) + '_merged.nii.gz',\
                '/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/PPI_DMN/randomise/noFeedback/' + motionDir +'/cope' + str(i) + '/' +t + str(i) + '_merged.nii.gz']
                merger = Merge()
                merger.inputs.in_files = x
                merger.inputs.dimension = 't'
                merger.inputs.output_type = 'NIFTI_GZ'
                merger.inputs.merged_file='./' + t + str(i)+'_merged.nii.gz'
                merger.run()
            except:
                print 'No Varcope'

        if runFlame:
            flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file='/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='design.mat',t_con_file='design.con', run_mode='flame1')
            flameo.run()
            foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/PPI_DMN/flame/paired-Ttest/' + motionDir + '/cope' + str(i)
            if os.path.exists(foldername):
                shutil.rmtree(foldername)
                os.mkdir(foldername)
            else:
                os.mkdir(foldername)
            if not runWithRandomise:
                shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
            shutil.move('varcope' + str(i) + '_merged.nii.gz',foldername)
            shutil.move('stats',foldername)
        if runWithRandomise:
            if not os.path.exists('cope%d' % i):
                os.mkdir('cope%d' % i)
            # shutil.move('./design.*','cope%d' % i)
            os.system('mv ./design.* ./cope%d' % i)
            randomiseCommand='./randomise_forpython.sh -i %s -o ./cope%d/cope%d -d ./cope%d/design.mat -t ./cope%d/design.con -e ./cope%d/design.grp -m %s -T -n %d' % ('cope' + str(i) + '_merged.nii.gz',i,i,i,i,i,'/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain_mask.nii.gz',nperms)
            os.system(randomiseCommand)

            foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/PPI_DMN/randomise/paired-Ttest/' + motionDir + '/cope' + str(i)
            if os.path.exists(foldername):
                shutil.rmtree(foldername)
                os.mkdir(foldername)
            else:
                os.mkdir(foldername)
            filename='cope%d' % i
            shutil.move(filename, os.path.join(foldername, filename))
            # shutil.move('cope%d' % i,foldername)
            filename='cope' + str(i) + '_merged.nii.gz'
            # shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
            shutil.move(filename, os.path.join(foldername, filename))
