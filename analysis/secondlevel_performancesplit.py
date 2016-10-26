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
parser.add_argument('-perfSplit', help='Option run by performance split (0-No Split,1-Top Tier,2-Middle Tier,3-Lowest Tier)',required=False,default=1,type=int)
parser.add_argument('-surface', help='Option to make surface plot (need to be on screen of computer running code)',required=False,default=0,type=int)
parser.add_argument('-runFC',help='Optiom to run FC instead of Cope', default=0,required=False,type=int)
parser.add_argument('-fc', help = 'Functional Connectivity ROI to run second level analysis on (overrides cope information)',required=False,default='R_AI',type=str)
parser.add_argument('-gmThresh',help='Grey Matter Threshold Value',default=0.2,type=float)
parser.add_argument('-train', help = 'Run RSN on train data not FB or NoFB',required=False,default=0,type=int)
parser.add_argument('-train_vs',help='Run train performance with FB or No FB',required=False,default=0,type=int)
parser.add_argument('-fbtorun', help = 'Which FB scans to run',required=False,nargs='+',default=[0,1],type=int)
parser.add_argument('-RSN', help='Option to run with RSN instead of cope, RSN>0)',required=False,default=0,type=int)
parser.add_argument('-combine', help='Option to Combine FB and NFB)',required=False,default=0,type=int)


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
fbtorun=args.fbtorun
train=args.train
RSN=args.RSN
combine=args.combine
traindiff=0
train_vs=2

mask_name='/home/jmuraskin/Projects/CCD/working_v1/seg_probabilities/grey_matter_mask-%d-percent.nii.gz' % int(gmThresh*100)

if not os.path.exists(mask_name):
    from nilearn.image import threshold_img
    mask_img=threshold_img('/home/jmuraskin/Projects/CCD/working_v1/seg_probabilities/grey_matter_mask.nii.gz',threshold=gmThresh)
    mask_img_data=mask_img.get_data()
    mask_img_data[mask_img_data>0]=1
    mask_img=new_img_like(mask_img,mask_img_data)
    mask_img.to_filename(mask_name)

if runFC:
    copesToRun=[0]
    rsn_name=''
else:
    fc=''


if surface:
    from CCD_packages import make_pysurfer_images

if RSN>0:
    rsn_name='RSN%d' % (RSN-1)
    rsn=RSN-1
    copesToRun=[0]
    fc=''
else:
    rsn_name=''

if train:
    fbtorun=[2]

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


def getSubjectButtonResponses():
    filelist=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/NARSAD_stimulus_JM.csv')

    for indx,f in enumerate(filelist['JM_INTERNAL']):
        for r in range(1,3):
            if int(f[-2:])<30:
                luminaFlag=0
            else:
                luminaFlag=1
            numberofbuttonPresses=getSubjectButtonPressScore('/home/jmuraskin/Projects/CCD/NARSAD-DMN-clean/%s_run%d.txt' % (f,r),luminaFlag)
            out={'number':numberofbuttonPresses,'filename':f}
            out['filename']=f
            if (indx+r)==1:
                df=pd.DataFrame(out,index=[0])
                df['subject']=f
                df['run']=r
            else:
                tmp=pd.DataFrame(out,index=[0])
                tmp['subject']=f
                tmp['run']=r
                df=pd.concat((df,tmp),ignore_index=0)
    return df


def getSubjectButtonPressScore(filename,luminaFlag):
    config=pd.read_table(filename,delimiter=';',comment='#')
    numButton=0
    for indx in config[config[' Stim Text']==' Push Button'].index[:]:
        numTmp=0
        for n in range(5):
            if luminaFlag:
                if config.iloc[indx+n][' STIM']==' LUMINA' and numTmp==0:
                    numButton+=1
                    numTmp+=1
            else:
                if config.iloc[indx+n][' STIM']!='53' and numTmp==0:
                    numButton+=1
                    numTmp+=1
    return numButton


#load subject list
motionTest=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_meanFD.csv')
# scanorderInfo=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_scanorder.csv',names=['Subject_ID','FB','meanFD'])
performance=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_performance.csv',names=['Subject_ID','FB','scanorder','R'])

fbNames=['NOFEEDBACK','FEEDBACK','TRAIN']

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

elif runAll==3:
    #Reject Depressed subjects
    depressed=np.array(['CCD072','CCD098','CCD083','CCD062','CCD061','CCD051','CCD087'])

    df=getSubjectButtonResponses()
    tmp=df.groupby('subject')['number'].sum()
    poor_performers=np.array(tmp[tmp<14].index[:])

    motionThresh=1
    allsubj=np.unique(motionTest['Subject_ID'])
    motionReject=np.unique((motionTest[motionTest.Max_Relative_RMS_Displacement>motionThresh]['Subject_ID']))
    subject_list=np.setdiff1d(np.setdiff1d(np.setdiff1d(allsubj,motionReject),depressed),poor_performers)
    motionDir='motionRMS-%f-subjperf' % motionThresh


else:
    motionThresh=0.2
    allsubj=np.unique(motionTest['Subject_ID'])
    motionReject=np.unique((motionTest[motionTest.meanFD>motionThresh]['Subject_ID']))
    subject_list=np.setdiff1d(allsubj,motionReject)
    motionDir='motionThresh-%f' % motionThresh


if perfSplit==1:
    #  sort by performance

    modelInfo=np.arctanh(performance[performance.FB=='FEEDBACK'][performance.Subject_ID.isin(subject_list)]['R'])
    NOFB=np.arctanh(performance[performance.FB=='NOFEEDBACK'][performance.Subject_ID.isin(subject_list)]['R'])
    top=performance.iloc[modelInfo[np.all([zscore(modelInfo)>0,zscore(NOFB)>0],axis=0)].index[:]]['Subject_ID']
    bottom=performance.iloc[modelInfo[np.all([zscore(modelInfo)<0,zscore(NOFB)<0],axis=0)].index[:]]['Subject_ID']


    subject_list=np.array(list(top) + list(bottom))
    topRegressor=[0]*len(subject_list)
    topRegressor[0:len(top)]=[1]*len(top)
    botRegressor=[0]*len(subject_list)
    botRegressor[len(top):]=[1]*len(bottom)
    perf_split_name ='performance_grouping_TopvBottom'


#load phenotypic data
phenoFile='/home/jmuraskin/Projects/CCD/Pheno/narsad+vt_new.csv'
pheno=pd.read_csv(phenoFile)
pheno=pheno.set_index('participant')
if age:
    ages=zscore(pheno.loc[subject_list]['V1_DEM_001'])
if gender:
    mf=zscore(pheno.loc[subject_list]['V1_DEM_002'])


secondlevel_folder_names=['noFeedback','Feedback','train','Combine']

#create second level folders
folderbase='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis'
for runType in ['randomise']:
    foldername=folderbase + '/' + runType + '/paired-Ttest/' +  motionDir + '/' + rsn_name + fc
    if not os.path.exists(foldername):
        os.makedirs(foldername)

    for fb in secondlevel_folder_names:
        foldername=folderbase + '/' + runType + '/' + fb + '/' + motionDir + '/' + rsn_name + fc
        if not os.path.exists(foldername):
            os.makedirs(foldername)




if run1Sample:

    for i in copesToRun:
        for fb in fbtorun:
            for t in ['cope']:
                x=[]
                for subj in subject_list:
                    if t=='cope' and RSN>0 and fb==2:
                        fname= '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/dr_tempreg_maps_files_to_standard_smooth/_scan_tra/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/_fwhm_6/_dr_tempreg_maps_files_smooth_0%d/temp_reg_map_000%d_antswarp_maths.nii.gz' % (subj,rsn,rsn)
                    else:
                        fbLoc=subjectinfo(subj,fb)
                        if t=='cope' and RSN>0:
                            fname= '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/dr_tempreg_maps_files_to_standard_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/_fwhm_6/_dr_tempreg_maps_files_smooth_0%d/temp_reg_map_000%d_antswarp_maths.nii.gz' % (subj,fbLoc+1,rsn,rsn)
                        elif len(fc)>0:
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

            if train:
                meanFD=zscore(motionTest[motionTest.FB==fbNames[0]][motionTest.Subject_ID.isin(subject_list)]['train_meanFD'])
            else:
                meanFD=zscore(motionTest[motionTest.FB==fbNames[fb]][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
            model = MultipleRegressDesign()
            model.inputs.contrasts = [['top>bottom', 'T',['top','bot'],[1,-1]],['bottom>top', 'T',['top','bot'],[-1,1]]]
            regressors=dict(top=topRegressor,bot=botRegressor,FD=list(meanFD))
            if age:
                regressors['age']=list(ages)
            if gender:
                regressors['mf']=list(mf)
            model.inputs.regressors = regressors
            model.run()

            if runWithRandomise:
                filename='cope%d' % i
                if age:
                    filename+='_age'
                if gender:
                    filename+='_gender'

                if perfSplit>0:
                    filename+=perf_split_name

                if not os.path.exists(filename):
                    os.makedirs(filename)
                os.system('mv ./design.* ./%s' % filename)
                os.system('mv cope%d_merged.nii.gz ./%s' % (i,filename))
                # shutil.move('./design.*','cope%d' % i)
                randomiseCommand='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/randomise_forpython.sh -i %s/%s -o ./%s/cope%d -d ./%s/design.mat -t ./%s/design.con -e ./%s/design.grp -m %s -T -n %d' % (filename,'cope' + str(i) + '_merged.nii.gz',filename,i,filename,filename,filename,mask_name,nperms)
                os.system(randomiseCommand)

                if RSN>0:
                    foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + rsn_name  + '/' + perf_split_name
                    if train and traindiff:
                        foldername=foldername + '/' + 'Difference'
                    elif train:
                        foldername=foldername + '/' + fbNames[train_vs]
                elif runFC:
                    foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + fc + '/' + perf_split_name
                    if train and traindiff:
                        foldername=foldername + '/' + 'Difference'
                    elif train:
                        foldername=foldername + '/' + fbNames[train_vs]
                else:
                    foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[fb] + '/' + motionDir + '/' + perf_split_name

                if not os.path.exists(foldername):
                    os.makedirs(foldername)

                if os.path.exists(os.path.join(foldername,filename)):
                    shutil.rmtree(os.path.join(foldername,filename))

                shutil.move(filename, os.path.join(foldername, filename))
                if surface:
                    make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)




if runPair:

    for i in copesToRun:

        if train:
            meanFD=zscore(motionTest[motionTest.FB==fbNames[0]][motionTest.Subject_ID.isin(subject_list)]['train_meanFD'])
        else:
            meanFD=zscore(motionTest[motionTest.FB==fbNames[fb]][motionTest.Subject_ID.isin(subject_list)]['meanFD'])
        model = MultipleRegressDesign()
        model.inputs.contrasts = [['top>bottom', 'T',['top','bot'],[1,-1]],['bottom>top', 'T',['top','bot'],[-1,1]]]
        regressors=dict(top=topRegressor,bot=botRegressor,FD=list(meanFD))
        if age:
            regressors['age']=list(ages)
        if gender:
            regressors['mf']=list(mf)
        model.inputs.regressors = regressors
        model.run()



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
                fslMathsCommand='fslmerge -t cope%d_merged %s %s' % (i,x[0],x[1])
                os.system(fslMathsCommand)
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


if combine:

    for i in copesToRun:
        for t in ['cope']:
            x=[]
            for subj in subject_list:
                for fb in [0,1]:
                    fbLoc=subjectinfo(subj,fb)
                    if t=='cope' and RSN>0:
                        fname= '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/dr_tempreg_maps_files_to_standard_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/_fwhm_6/_dr_tempreg_maps_files_smooth_0%d/temp_reg_map_000%d_antswarp_maths.nii.gz' % (subj,fbLoc+1,rsn,rsn)
                    elif len(fc)>0:
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

        meanFD=list(zscore(motionTest[motionTest.FB==fbNames[0]][motionTest.Subject_ID.isin(subject_list)]['meanFD']))+list(zscore(motionTest[motionTest.FB==fbNames[1]][motionTest.Subject_ID.isin(subject_list)]['meanFD']))
        model = MultipleRegressDesign()
        model.inputs.contrasts = [['top>bottom', 'T',['top','bot'],[1,-1]],['bottom>top', 'T',['top','bot'],[-1,1]]]
        regressors=dict(top=topRegressor,bot=botRegressor,FD=list(meanFD))
        if age:
            regressors['age']=list(ages)+list(ages)
        if gender:
            regressors['mf']=list(mf)+list(mf)
        model.inputs.regressors = regressors
        model.run()

        if runWithRandomise:
            filename='cope%d' % i
            if age:
                filename+='_age'
            if gender:
                filename+='_gender'

            if perfSplit>0:
                filename+=perf_split_name

            filename+='_Grouped'

            if not os.path.exists(filename):
                os.makedirs(filename)
            os.system('mv ./design.* ./%s' % filename)
            os.system('mv cope%d_merged.nii.gz ./%s' % (i,filename))
            # shutil.move('./design.*','cope%d' % i)
            randomiseCommand='/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/randomise_forpython.sh -i %s/%s -o ./%s/cope%d -d ./%s/design.mat -t ./%s/design.con -e ./%s/design.grp -m %s -T -n %d' % (filename,'cope' + str(i) + '_merged.nii.gz',filename,i,filename,filename,filename,mask_name,nperms)
            os.system(randomiseCommand)

            if RSN>0:
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[3] + '/' + motionDir + '/' + rsn_name  + '/' + perf_split_name
                if train and traindiff:
                    foldername=foldername + '/' + 'Difference'
                elif train:
                    foldername=foldername + '/' + fbNames[train_vs]
            elif runFC:
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[3] + '/' + motionDir + '/' + fc + '/' + perf_split_name
                if train and traindiff:
                    foldername=foldername + '/' + 'Difference'
                elif train:
                    foldername=foldername + '/' + fbNames[train_vs]
            else:
                foldername='/home/jmuraskin/Projects/CCD/working_v1/groupAnalysis/randomise/' + secondlevel_folder_names[3] + '/' + motionDir + '/' + perf_split_name

            if not os.path.exists(foldername):
                os.makedirs(foldername)

            if os.path.exists(os.path.join(foldername,filename)):
                shutil.rmtree(os.path.join(foldername,filename))

            shutil.move(filename, os.path.join(foldername, filename))
            if surface:
                make_pysurfer_images(folder=os.path.join(foldername, filename),suffix='cope%d' % i)
