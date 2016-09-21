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
import sys
import argparse
from CCD_packages import make_pysurfer_images





# class Usage(Exception):
#     def __init__(self, msg):
#         self.msg = msg
#
# def main(argv=None):
#     if argv is None:
#         argv = sys.argv
#     try:
#         try:
#             opts, args = getopt.getopt(argv[1:], "h", ["help"])
#         except getopt.error, msg:
#              raise Usage(msg)
#         # more code, unchanged
#     except Usage, err:
#         print >>sys.stderr, err.msg
#         print >>sys.stderr, "for help use --help"
#         return 2
#
# if __name__ == "__main__":
#     sys.exit(main())

parser = argparse.ArgumentParser(description='Run Second Level Results for CCD')
parser.add_argument('-rwr', help='Option to run with Randomise',required=False,default=1,type=int)
parser.add_argument('-rwf', help='Option to run with FLAME',required=False,default=0,type=int)
parser.add_argument('-n',help='Number of Permutations to Run', required=False,default=10000,type=int)
parser.add_argument('-r1samp', help='Option to run 1 sample t-test',required=False,default=1,type=int)
parser.add_argument('-rpair', help='Option to run paired t-test',required=False,default=1,type=int)
parser.add_argument('-rall', help='Option to run all subjects or good motion subjects',required=False,default=1,type=int)
parser.add_argument('-copes', help='List of copes to run',nargs='+', type=int,required=False,default=[1])
parser.add_argument('-a', help='Option to add subject age to model',required=False,default=0,type=int)
parser.add_argument('-g', help='Option to add subject gender to model',required=False,default=0,type=int)
parser.add_argument('-perfSplit', help='Option run by performance split (0-No Split,1-Top Tier,2-Middle Tier,3-Lowest Tier)',required=False,default=0,type=int)
parser.add_argument('-ROIname', help='ROI Foldername',required=True,default='AI_ROI_ts_GSR1',type=str)
parser.add_argument('-pheno', help='Phenotype Measure to Run', type=str,required=False,default='V1_CCDRSQ_75')
parser.add_argument('-perf',help='Run with Performance instead of Phenotype',type=int,required=False,default=0)


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
roiFolder=args.ROIname
runWithPerformance=args.perf


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
# scanorderInfo=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_scanorder.csv',names=['Subject_ID','FB','meanFD'])
performance=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_performance.csv',names=['Subject_ID','FB','scanorder','R'])

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
pheno_measure_name=args.pheno
pheno_measure = zscore(pheno.loc[subject_list][pheno_measure_name])
if runWithPerformance:
    pheno_measure_name='Performance'

if age:
    ages=zscore(pheno.loc[subject_list]['V1_DEM_001'])
if gender:
    mf=zscore(pheno.loc[subject_list]['V1_DEM_002'])
# pheno_measure = zscore(pheno.loc[subject_list][pheno_measure_name])

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
folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis'
for runType in ['randomise','flame']:
    foldername=folderbase + '/' + runType + '/paired-Ttest/' +  motionDir + '/' + roiFolder
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    foldername=folderbase + '/' + runType + '/paired-Ttest/' +  motionDir + '/' + roiFolder + '/' + pheno_measure_name
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    for fb in secondlevel_folder_names:
        foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir + '/' + roiFolder
        if not os.path.exists(foldername):
            os.mkdir(foldername)

        foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir + '/' + roiFolder + '/' + pheno_measure_name
        if not os.path.exists(foldername):
            os.mkdir(foldername)




if run1Sample:

    for i in copesToRun:
        for fb in [0,1]:
            for t in ['cope', 'varcope']:
                x=[]
                for subj in subject_list:
                    fbLoc=subjectinfo(subj,fb)
                    # fname = '/home/jmuraskin/Projects/CCD/working_v1/feedback_run-%d/feedback/_subject_id_%s/modelestimate/mapflow/_modelestimate0/results/%s%d.nii.gz' % (fbLoc,subj,t,i)
                    fname = '/home/jmuraskin/Projects/CCD/working_v1/%s/feedback_run-%d/ROI_model_fit/_subject_id_%s/modelestimate/mapflow/_modelestimate0/results/%s%d.nii.gz' % (roiFolder,fbLoc,subj,t,i)
                    x.append(fname)
                subjs = len(x)
                merger = Merge()
                merger.inputs.in_files = x
                merger.inputs.dimension = 't'
                merger.inputs.output_type = 'NIFTI_GZ'
                merger.run()
            #get meanFD values for each subject and add as covariate
            meanFD=zscore(motionTest[motionTest.FB==fbNames[fb]][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
            if runWithPerformance:
                pheno_measure = zscore(np.arctan(performance[performance.FB==fbNames[fb]][performance.Subject_ID.isin(subject_list)]['R']))
            model = MultipleRegressDesign()
            model.inputs.contrasts = [['pheno pos', 'T',['pheno'],[1]],['pheno neg', 'T',['pheno'],[-1]]]
            regressors=dict(pheno=list(pheno_measure),FD=list(meanFD))
            if age:
                regressors['age']=list(ages)
            if gender:
                regressors['mf']=list(mf)
            model.inputs.regressors = regressors
            model.run()

            if runFlame:
                flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file='/home/jmuraskin/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='design.mat',t_con_file='design.con', run_mode='flame1')
                flameo.run()
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/flame/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + roiFolder + '/cope' + str(i)
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
                randomiseCommand='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/randomise_forpython.sh -i %s/%s -o ./%s/cope%d -d ./%s/design.mat -t ./%s/design.con -e ./%s/design.grp -m %s -T -n %d' % (filename,'cope' + str(i) + '_merged.nii.gz',filename,i,filename,filename,filename,'/home/jmuraskin/standard/MNI152_T1_3mm_brain_mask.nii.gz',nperms)
                os.system(randomiseCommand)

                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + roiFolder + '/' + pheno_measure_name

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
                make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)




if runPair:

    for i in copesToRun:

        meanFD=zscore(np.array(motionTest[motionTest.FB=='FEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD'])+np.array(motionTest[motionTest.FB=='NOFEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD']))
        model = MultipleRegressDesign()
        model.inputs.contrasts = [['pheno pos', 'T',['pheno'],[1]],['pheno neg', 'T',['pheno'],[-1]]]
        if runWithPerformance:
            pheno_measure = zscore(np.array(np.arctan(performance[performance.FB=='FEEDBACK'][performance.Subject_ID.isin(subject_list)]['R']))-np.array(np.arctan(performance[performance.FB=='NOFEEDBACK'][performance.Subject_ID.isin(subject_list)]['R'])))
        regressors=dict(pheno=list(pheno_measure),FD=list(meanFD))
        if age:
            regressors['age']=list(ages)
        if gender:
            regressors['mf']=list(mf)
        model.inputs.regressors = regressors
        model.run()

        if age:
            addage='_age'
        else:
            addage=''
        if gender:
            addgender='_gender'
        else:
            addgender=''

        x=['/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/Feedback/' + motionDir + '/' + roiFolder + '/' + pheno_measure_name +  '/cope' + str(i) + addage + addgender + '/cope' + str(i) + '_merged.nii.gz',\
        '/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/noFeedback/' + motionDir + '/' + roiFolder + '/' + pheno_measure_name +  '/cope' + str(i) + addage + addgender + '/cope' + str(i) + '_merged.nii.gz']


        fslMathsCommand='fslmaths %s -sub %s cope%d_merged' % (x[0],x[1],i)
        os.system(fslMathsCommand)

        if runWithRandomise:
            filename='cope%d' % i
            filename+='_%s' % pheno_measure_name
            if age:
                filename+='_age'
            if gender:
                filename+='_gender'

            if not os.path.exists(filename):
                os.mkdir(filename)
            os.system('mv ./design.* ./%s' % filename)
            os.system('mv cope%d_merged.nii.gz ./%s' % (i,filename))
            # shutil.move('./design.*','cope%d' % i)
            randomiseCommand='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/randomise_forpython.sh -i %s/%s -o ./%s/cope%d -d ./%s/design.mat -t ./%s/design.con -e ./%s/design.grp -m %s -T -n %d -D' % (filename,'cope' + str(i) + '_merged.nii.gz',filename,i,filename,filename,filename,'/home/jmuraskin/standard/MNI152_T1_3mm_brain_mask.nii.gz',nperms)
            os.system(randomiseCommand)

            foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/paired-Ttest/' + motionDir + '/'  + roiFolder + '/' + pheno_measure_name
            print 'Making folder: %s' % foldername
            if not os.path.exists(foldername):
                os.mkdir(foldername)

            if os.path.exists(os.path.join(foldername,filename)):
                shutil.rmtree(os.path.join(foldername,filename))

            shutil.move(filename,os.path.join(foldername,filename))

            make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)
