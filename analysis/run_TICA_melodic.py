import nipype.interfaces.fsl as fsl
import pandas as pd
from CCD_packages import fb_subjectinfo




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



motionTest=pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_meanFD.csv')
depressed=np.array(['CCD072','CCD098','CCD083','CCD062','CCD061','CCD051','CCD087'])
df=getSubjectButtonResponses()
tmp=df.groupby('subject')['number'].sum()
poor_performers=np.array(tmp[tmp<22].index[:])


motionThresh=1
allsubj=np.unique(motionTest['Subject_ID'])
motionReject=np.unique((motionTest[motionTest.Max_Relative_RMS_Displacement>motionThresh]['Subject_ID']))
subject_list=np.setdiff1d(np.setdiff1d(np.setdiff1d(allsubj,motionReject),depressed),poor_performers)

gmThresh=0.20

fnames=[]
fb=1
for subject in subject_list:
    run=fb_subjectinfo(subject,fb)+1
    fnames.append('/home/jmuraskin/Projects/CCD/working_v1/feedback_run-%d/_subject_id_%s/addMeanImage/mapflow/_addMeanImage0/residual_antswarp_maths_maths.nii.gz' % (run,subject))


melodic_setup = fsl.model.MELODIC()
melodic_setup.inputs.approach = 'symm'
melodic_setup.inputs.in_files = fnames
melodic_setup.inputs.no_bet = True
melodic_setup.inputs.bg_threshold = 10
melodic_setup.inputs.mask = '/home/jmuraskin/Projects/CCD/working_v1/seg_probabilities/grey_matter_mask-%d-percent.nii.gz' % int(gmThresh*100)
melodic_setup.inputs.tr_sec = 2
melodic_setup.inputs.mm_thresh = 0.5
melodic_setup.inputs.out_stats = True
melodic_setup.inputs.t_des = 'TENSOR_GLM.mat'
melodic_setup.inputs.t_con = 'TENSOR_GLM.con'
melodic_setup.inputs.s_des = 'design-Feedback-perf.mat'
melodic_setup.inputs.s_con = 'design-Feedback-perf.con'
melodic_setup.inputs.out_dir = '/home/jmuraskin/Projects/CCD/working_v1/group_FEEDBACK_TICA.out'
melodc_setup.inputs.report = True

melodic_setup.run()
