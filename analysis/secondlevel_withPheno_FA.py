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


parser = argparse.ArgumentParser(description='Run Second Level Results with PhenoTypic Measures')
parser.add_argument('-rwr', help='Option to run with Randomise',required=False,default=1,type=int)
parser.add_argument('-rwf', help='Option to run with FLAME',required=False,default=0,type=int)
parser.add_argument('-n',help='Number of Permutations to Run', required=False,default=10000,type=int)
parser.add_argument('-r1samp', help='Option to run 1 sample t-test',required=False,default=1,type=int)
parser.add_argument('-rpair', help='Option to run paired t-test',required=False,default=0,type=int)
parser.add_argument('-rall', help='Option to run all subjects or good motion subjects',required=False,default=1,type=int)
parser.add_argument('-copes', help='List of copes to run',nargs='+', type=int,required=False,default=[1,3,4,5])
# parser.add_argument('-pheno', help='Phenotype Measure to Run', type=str,required=False,default='V1_CCDRSQ_75')
# parser.add_argument('-perf',help='Run with Performance instead of Phenotype',type=int,required=False,default=0)
parser.add_argument('-a', help='Option to add subject age to model',required=False,default=0,type=int)
parser.add_argument('-g', help='Option to add subject gender to model',required=False,default=0,type=int)
parser.add_argument('-surface', help='Option to make surface plot (need to be on screen of computer running code)',required=False,default=0,type=int)
parser.add_argument('-RSN', help='Option to run with RSN instead of cope, RSN>0)',required=False,default=0,type=int)

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
# runWithPerformance=args.perf
age=args.a
gender=args.g
surface=args.surface
RSN=args.RSN

if RSN>0:
    rsn_name='RSN%d' % (RSN-1)
    rsn=RSN-1
    copesToRun=[0]
else:
    rsn_name=''

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

#
# #Decide if running all subjects or just good subjects
# runAll=True

#load subject list
motionTest=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_meanFD.csv')
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


#load phenotypic data
phenoFile='/home/jmuraskin/Projects/CCD/Pheno/narsad+vt_new.csv'
pheno=pd.read_csv(phenoFile)
pheno=pheno.set_index('participant')

pheno_measure_name='Factor_Analysis'

FA_file='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/FA_subjectloadings.csv'

FA=pd.read_csv(FA_file,index_col=0)
#

# # pheno_measure = zscore(pheno.loc[subject_list][pheno_measure_name])
# if runWithPerformance:
#     pheno_measure_name='Performance'

if age:
    ages=zscore(pheno.loc[subject_list]['V1_DEM_001'])
if gender:
    mf=zscore(pheno.loc[subject_list]['V1_DEM_002'])


secondlevel_folder_names=['noFeedback','Feedback']

#create second level folders
folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis'
for runType in ['randomise','flame']:
    foldername=folderbase + '/' + runType + '/paired-Ttest/' +  motionDir + '/' + rsn_name
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    foldername= foldername + '/' + pheno_measure_name
    if not os.path.exists(foldername):
        os.mkdir(foldername)

    for fb in secondlevel_folder_names:
        if RSN>0:
            foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir + '/' + rsn_name
        else:
            foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir
        if not os.path.exists(foldername):
            os.mkdir(foldername)
        foldername= foldername + '/' + pheno_measure_name
        if not os.path.exists(foldername):
            os.mkdir(foldername)


if run1Sample:


    for i in copesToRun:
        for fb in [0,1]:
            for t in ['cope']:
                x=[]
                for subj in subject_list:
                    fbLoc=subjectinfo(subj,fb)
                    if t=='cope' and RSN>0:
                        fname= '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/dr_tempreg_maps_files_to_standard_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/_fwhm_6/_dr_tempreg_maps_files_smooth_0%d/temp_reg_map_000%d_antswarp_maths.nii.gz' % (subj,fbLoc+1,rsn,rsn)
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
            # if runWithPerformance:
            #     pheno_measure = zscore(np.arctan(performance[performance.FB==fbNames[fb]][performance.Subject_ID.isin(subject_list)]['R']))
            model = MultipleRegressDesign()
            model.inputs.contrasts = [['pheno pos', 'T',['comp1'],[1]],['pheno neg', 'T',['comp1'],[-1]],['pheno 2 pos', 'T',['comp2'],[1]],['pheno 2 neg', 'T',['comp2'],[-1]]]
            regressors=dict(comp1=list(FA.loc[subject_list]['FirstComp']),comp2=list(FA.loc[subject_list]['SecondComp']),FD=list(meanFD))
            if age:
                regressors['age']=list(ages)
            if gender:
                regressors['mf']=list(mf)
            model.inputs.regressors = regressors
            model.run()

            if runFlame:
                flameo = fsl.FLAMEO(cope_file='./cope'+str(i)+'_merged.nii.gz',var_cope_file='./varcope'+str(i)+'_merged.nii.gz',cov_split_file='design.grp',mask_file='/home/jmuraskin/standard/MNI152_T1_3mm_brain_mask.nii.gz',design_file='design.mat',t_con_file='design.con', run_mode='flame1')
                flameo.run()
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/flame/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + pheno_measure_name + '/cope' + str(i)
                if os.path.exists(foldername):
                    shutil.rmtree(foldername)
                    os.mkdir(foldername)
                else:
                    os.mkdir(foldername)
                if not runWithRandomise:
                    shutil.move('cope' + str(i) + '_merged.nii.gz',foldername)
                shutil.move('stats',foldername + '/stats')
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
                if RSN>0:
                    foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + rsn_name + '/' + pheno_measure_name
                else:
                    foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + pheno_measure_name

                if not os.path.exists(foldername):
                    os.mkdir(foldername)

                if os.path.exists(os.path.join(foldername,filename)):
                    shutil.rmtree(os.path.join(foldername,filename))

                shutil.move(filename,foldername+ '/' + filename )
                if surface:
                    make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)





if runPair:


    for i in copesToRun:


        meanFD=zscore(np.array(motionTest[motionTest.FB=='FEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD'])+np.array(motionTest[motionTest.FB=='NOFEEDBACK'][motionTest.Subject_ID.isin(subject_list)]['meanFD']))
        model = MultipleRegressDesign()
        model.inputs.contrasts = [['pheno pos', 'T',['comp1'],[1]],['pheno neg', 'T',['comp1'],[-1]],['pheno 2 pos', 'T',['comp2'],[1]],['pheno 2 neg', 'T',['comp2'],[-1]]]
        regressors=dict(comp1=list(FA.loc[subject_list]['FirstComp']),comp2=list(FA.loc[subject_list]['SecondComp']),FD=list(meanFD))
        if age:
            regressors['age']=list(ages)
        if gender:
            regressors['mf']=list(mf)
        model.inputs.regressors = regressors
        model.run()

        x=['/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/Feedback/' + motionDir +'/cope' + str(i) + '/cope' + str(i) + '_merged.nii.gz',\
        '/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/noFeedback/' + motionDir +'/cope' + str(i) + '/cope' + str(i) + '_merged.nii.gz']


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

            foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/paired-Ttest/' + motionDir + '/' + pheno_measure_name
            print 'Making folder: %s' % foldername
            if not os.path.exists(foldername):
                os.mkdir(foldername)

            if os.path.exists(os.path.join(foldername,filename)):
                shutil.rmtree(os.path.join(foldername,filename))

            shutil.move(filename,os.path.join(foldername,filename))
            if surface:
                make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)
            # shutil.move('cope%d_pair_diff.nii.gz' % i,foldername)
