import os
import shutil
import glob
import pandas as pd
import numpy as np
from nipype.interfaces.fsl import Merge
from nipype.interfaces import fsl
from subprocess import call
from nipype.interfaces.fsl import MultipleRegressDesign
from nilearn.image import threshold_img
from scipy.stats import zscore
import sys
import argparse



parser = argparse.ArgumentParser(description='Run Second Level Results for CCD')
parser.add_argument('-rwr', help='Option to run with Randomise',required=False,default=1,type=int)
parser.add_argument('-rwf', help='Option to run with FLAME',required=False,default=0,type=int)
parser.add_argument('-n',help='Number of Permutations to Run', required=False,default=10000,type=int)
parser.add_argument('-r1samp', help='Option to run 1 sample t-test',required=False,default=1,type=int)
parser.add_argument('-rpair', help='Option to run paired t-test',required=False,default=1,type=int)
parser.add_argument('-rall', help='Option to run all subjects==1, option run subjects based on Max Relative RMS<1==2, option to run only Mean FD<0.2==0',required=False,default=1,type=int)
parser.add_argument('-copes', help='List of copes to run',nargs='+', type=int,required=False,default=range(5))
parser.add_argument('-a', help='Option to add subject age to model',required=False,default=0,type=int)
parser.add_argument('-g', help='Option to add subject gender to model',required=False,default=0,type=int)
parser.add_argument('-perfSplit', help='Option run by performance split (0-No Split,1-Top Tier,2-Middle Tier,3-Lowest Tier)',required=False,default=0,type=int)
parser.add_argument('-surface', help='Option to make surface plot (need to be on screen of computer running code)',required=False,default=0,type=int)
parser.add_argument('-runFC',help='Optiom to run FC instead of Cope', default=0,required=False,type=int)
parser.add_argument('-fc', help = 'Functional Connectivity ROI to run second level analysis on (overrides cope information)',required=False,default='R_AI',type=str)
parser.add_argument('-gmThresh',help='Grey Matter Threshold Value',default=0.2,type=float)
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
age=args.a
gender=args.g
perfSplit=args.perfSplit
surface=args.surface
fc=args.fc
runFC=args.runFC
gmThresh=args.gmThresh




mask_img=threshold_img('/home/jmuraskin/Projects/CCD/working_v1/seg_probabilities/grey_matter_mask.nii.gz',threshold=gmThresh)
mask_name='/home/jmuraskin/Projects/CCD/working_v1/seg_probabilities/grey_matter_mask-%d-percent.nii.gz' % gmThresh*100
mask_img.to_filename(mask_name)

if runFC:
    copesToRun=[0]
else:
    fc=''

if surface:
    from CCD_packages import make_pysurfer_images

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
motionTest=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_meanFD.csv')
# scanorderInfo=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_scanorder.csv',names=['Subject_ID','FB','meanFD'])
performance=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_performance.csv',names=['Subject_ID','FB','scanorder','R'])

fbNames=['NOFEEDBACK','FEEDBACK']

if runAll==1:
    subject_list=np.unique(motionTest.Subject_ID)
    motionDir='all'
elif runAll==2:
    depressed=np.array(['CCD072','CCD098','CCD083','CCD062','CCD061','CCD051','CCD087'])
    motionThresh=1
    allsubj=np.unique(motionTest['Subject_ID'])
    motionReject=np.unique((motionTest[motionTest.Max_Relative_RMS_Displacement>motionThresh]['Subject_ID']))
    subject_list=np.setdiff1d(np.setdiff1d(allsubj,motionReject),depressed)
    motionDir='motionRMS-%f' % motionThresh
else:
    motionThresh=0.2
    allsubj=np.unique(motionTest['Subject_ID'])
    motionReject=np.unique((motionTest[motionTest.meanFD>motionThresh]['Subject_ID']))
    subject_list=np.setdiff1d(allsubj,motionReject)
    motionDir='motionThresh-%f' % motionThresh


if perfSplit>0:
    # sort by performance
    maxModel=performance[performance.Subject_ID.isin(subject_list)].groupby(['Subject_ID'])['R'].max().sort_values(ascending=False)
    sortedOrder=maxModel.index
    numSubjs=len(sortedOrder)
    numSubjsPerGroup=numSubjs/3
    if perfSplit==1:
        subject_list=np.array(sortedOrder[0:numSubjsPerGroup])
        perf_split_name='_Tier-1'
    if perfSplit==2:
        subject_list=np.array(sortedOrder[numSubjsPerGroup+1:2*numSubjsPerGroup])
        perf_split_name='_Tier-2'
    if perfSplit==3:
        subject_list=np.array(sortedOrder[-numSubjsPerGroup:])
        perf_split_name='_Tier-3'

#load phenotypic data
phenoFile='/home/jmuraskin/Projects/CCD/Pheno/narsad+vt_new.csv'
pheno=pd.read_csv(phenoFile)
pheno=pheno.set_index('participant')
if age:
    ages=zscore(pheno.loc[subject_list]['V1_DEM_001'])
if gender:
    mf=zscore(pheno.loc[subject_list]['V1_DEM_002'])


secondlevel_folder_names=['noFeedback','Feedback']

#create second level folders
folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis'
for runType in ['randomise','flame']:
    foldername=folderbase + '/' + runType + '/paired-Ttest/' +  motionDir + '/' + fc
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    for fb in secondlevel_folder_names:
        foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir + '/' + fc
        if not os.path.exists(foldername):
            os.mkdir(foldername)




if run1Sample:

    for i in copesToRun:
        for fb in [0,1]:
            for t in ['cope']:
                x=[]
                for subj in subject_list:
                    fbLoc=subjectinfo(subj,fb)
                    if len(fc)>0:
                        fname= '/home/jmuraskin/Projects/CCD/working_v1/seed-to-voxel/%s/%s/%s_%s.nii.gz' % (fc,secondlevel_folder_names[fb],fc,subj)
                    else:
                        fname = '/home/jmuraskin/Projects/CCD/working_v1/feedback_run-%d/feedback/_subject_id_%s/modelestimate/mapflow/_modelestimate0/results/%s%d.nii.gz' % (fbLoc,subj,t,i)
                    x.append(fname)
                subjs = len(x)
                merger = Merge()
                merger.inputs.in_files = x
                merger.inputs.dimension = 't'
                merger.inputs.output_type = 'NIFTI_GZ'
                merger.inputs.merged_file = './cope%d_merged.nii.gz' % i
                merger.run()
            #get meanFD values for each subject and add as covariate
            meanFD=zscore(motionTest[motionTest.FB==fbNames[fb]][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
            model = MultipleRegressDesign()
            model.inputs.contrasts = [['group mean', 'T',['reg1'],[1]],['group neg mean', 'T',['reg1'],[-1]]]
            regressors=dict(reg1=[1]*len(subject_list),FD=list(meanFD))
            if age:
                regressors['age']=list(ages)
            if gender:
                regressors['mf']=list(mf)
            model.inputs.regressors = regressors
            model.run()

            if runFlame:
                flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file=mask_name,design_file='design.mat',t_con_file='design.con', run_mode='flame1')
                flameo.run()
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/flame/' + secondlevel_folder_names[fb] + '/' + motionDir + '/cope' + str(i)
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
                filename='cope%d' % i
                if age:
                    filename+='_age'
                if gender:
                    filename+='_gender'

                if perfSplit>0:
                    filename+=perf_split_name

                if not os.path.exists(filename):
                    os.mkdir(filename)
                os.system('mv ./design.* ./%s' % filename)
                os.system('mv cope%d_merged.nii.gz ./%s' % (i,filename))
                # shutil.move('./design.*','cope%d' % i)
                randomiseCommand='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/randomise_forpython.sh -i %s/%s -o ./%s/cope%d -d ./%s/design.mat -t ./%s/design.con -e ./%s/design.grp -m %s -T -n %d' % (filename,'cope' + str(i) + '_merged.nii.gz',filename,i,filename,filename,filename,mask_name,nperms)
                os.system(randomiseCommand)

                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + fc

                # if age:
                #     foldername+='_age'
                # if gender:
                #     foldername+='_gender'
                # if perfSplit>0:
                #     foldername+=perf_split_name

                if not os.path.exists(foldername):
                    os.mkdir(foldername)

                if os.path.exists(os.path.join(foldername,filename)):
                    shutil.rmtree(os.path.join(foldername,filename))

                shutil.move(filename, os.path.join(foldername, filename))
                if surface:
                    make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)




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
        if addScanOrder:
            modelDict['scanorder']= list(zscore(list(motionTest[motionTest.FB=='FEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['scanorder'])
        + list(motionTest[motionTest.FB=='NOFEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['scanorder'])))
        if age:
            modelDict['age']=list(ages)+list(ages)
        if gender:
            modelDict['mf']=list(mf)+list(mf)
        pairedmodel.inputs.regressors = modelDict
        pairedmodel.run()



        for t in ['cope']:
            try:
                feedbackFile='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/Feedback/' + motionDir + '/' + fc + '/cope' + str(i)
                if age:
                    feedbackFile+='_age'
                if gender:
                    feedbackFile+='_gender'
                if perfSplit>0:
                    feedbackFile+=perf_split_name
                feedbackFile+= '/' + t + str(i) + '_merged.nii.gz'
                nofeedbackFile='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/noFeedback/' + motionDir + '/' + fc + '/cope' + str(i)
                if age:
                    nofeedbackFile+='_age'
                if gender:
                    nofeedbackFile+='_gender'
                if perfSplit>0:
                    nofeedbackFile+=perf_split_name
                nofeedbackFile+= '/' + t + str(i) + '_merged.nii.gz'
                x=[feedbackFile,nofeedbackFile]
                merger = Merge()
                merger.inputs.in_files = x
                merger.inputs.dimension = 't'
                merger.inputs.output_type = 'NIFTI_GZ'
                merger.inputs.merged_file='./' + t + str(i)+'_merged.nii.gz'
                merger.run()
            except:
                print 'No Varcope'

        if runFlame:
            flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file=mask_name,design_file='design.mat',t_con_file='design.con', run_mode='flame1')
            flameo.run()
            foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/flame/paired-Ttest/' + motionDir + '/cope' + str(i)
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
            filename='cope%d' % i
            if age:
                filename+='_age'
            if gender:
                filename+='_gender'
            if perfSplit>0:
                filename+=perf_split_name

            if not os.path.exists(filename):
                os.mkdir(filename)
            os.system('mv ./design.* ./%s' % filename)
            os.system('mv cope%d_merged.nii.gz ./%s' % (i,filename))
            # shutil.move('./design.*','cope%d' % i)
            randomiseCommand='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/randomise_forpython.sh -i %s/%s -o ./%s/cope%d -d ./%s/design.mat -t ./%s/design.con -e ./%s/design.grp -m %s -T -n %d' % (filename,'cope' + str(i) + '_merged.nii.gz',filename,i,filename,filename,filename,mask_name,nperms)
            os.system(randomiseCommand)


            if addScanOrder:
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/paired-Ttest/' + motionDir + '/so_cope' + str(i)
            else:
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/paired-Ttest/' + motionDir + '/' + fc

            # if age:
            #     foldername+='_age'
            # if gender:
            #     foldername+='_gender'
            # if perfSplit>0:
            #     foldername+=perf_split_name

            if not os.path.exists(foldername):
                os.mkdir(foldername)

            if os.path.exists(os.path.join(foldername,filename)):
                shutil.rmtree(os.path.join(foldername,filename))

            shutil.move(filename, os.path.join(foldername, filename))
            if surface:
                make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)
