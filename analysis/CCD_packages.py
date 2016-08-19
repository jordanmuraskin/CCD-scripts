import networkx as nx
import pandas as pd
import numpy as np
from numpy import unique
from scipy.stats import zscore,spearmanr,pearsonr
import seaborn as sns
import matplotlib.pylab as plt
import os.path
from scipy.signal import butter,filtfilt

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y



def getCCDSubjectData(filterOn=False,zscoreOn=True,lowpass=0.1,globalNR=0,saveMotionInfo=False,verbose=False):


    SubjInfo = pd.read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/NARSAD_stimulus_JM.csv')
    # SubjInfo.set_index('JM_INTERNAL',inplace=True)
    dmnIdeal=pd.read_csv('/home/jmuraskin/Projects/NFB/analysis/DMN_ideal_2.csv')

    drFileLocation='/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1'

    GroupDF=[]
    numberOfICs=10
    columnNames=[]
    for rsnNumber in range(numberOfICs):
        columnNames.append('RSN%d' % rsnNumber)

    for indx,row in SubjInfo.iterrows():

        subj = row['JM_INTERNAL']
        if verbose:
            print 'Collecting Subject %s' % subj
        for scan in range(1,3,1):
            drFilePath = '%s/%s_data_/spatial_map_timeseries_for_DR/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global%d.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/spatial_map_timeseries.txt' % (drFileLocation,subj,scan,globalNR)
            df=[]
            subjHasBoth=True
            if os.path.isfile('%s/%s_data_/spatial_map_timeseries_for_DR/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global%d.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/spatial_map_timeseries.txt' % (drFileLocation,subj,1,globalNR)) and os.path.isfile('%s/%s_data_/spatial_map_timeseries_for_DR/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global%d.motion1.quadratic1.gm0.compcor1.csf1/_spatial_map_PNAS_Smith09_rsn10/spatial_map_timeseries.txt' % (drFileLocation,subj,2,globalNR)):
                try:
                    df = pd.read_csv(drFilePath,header=None,names=columnNames,delim_whitespace=True)
                    df['Subject_ID'] = subj
                    df['Subject'] = indx
                    df.index.name = 'TR'
                    df.reset_index(level=0,inplace=True)
                    if row['SCAN_%d_PARADIGM' % scan]==1 or row['SCAN_%d_PARADIGM' % scan]==3:
                        for rsn in columnNames:
                            if filterOn:
                                if zscoreOn:
                                    df[rsn]=pd.Series(-1*zscore(butter_lowpass_filter(df[rsn][:],lowpass,0.5)))
                                else:
                                    df[rsn]=pd.Series(-1*butter_lowpass_filter(df[rsn][:],lowpass,0.5))

                            else:
                                if zscoreOn:
                                    df[rsn]=pd.Series(-1*zscore(df[rsn][:]))
                                else:
                                    df[rsn]=pd.Series(-1*df[rsn][:])


                        df['flip']=-1
                    else:
                        for rsn in columnNames:
                            if filterOn:
                                if zscoreOn:
                                    df[rsn]=pd.Series(zscore(butter_lowpass_filter(df[rsn][:],lowpass,0.5)))
                                else:
                                    df[rsn]=pd.Series(butter_lowpass_filter(df[rsn][:],lowpass,0.5))

                            else:
                                if zscoreOn:
                                    df[rsn]=pd.Series(zscore(df[rsn][:]))
                                else:
                                    df[rsn]=pd.Series(df[rsn][:])


                        df['flip']=1
                    df['FB'] = 'FEEDBACK' if row['SCAN_%d_FEEDBACK' % scan]==1 else 'NOFEEDBACK'
                    df['scanorder']=scan
                    df['modelcorr']=pearsonr(dmnIdeal['Wander']-dmnIdeal['Focus'],df['RSN3'])[0]
    #                 df['DMN']=pd.Series(zscore(nuisanceRegression(df[list(set(columnNames)-set(['RSN3']))],df['RSN3'])))
                    #load meanFD scores
                    fdFilePath='%s/%s_data_/frame_wise_displacement/_scan_feedback_%d/FD.1D' % (drFileLocation,subj,scan)
                    fd=pd.read_csv(fdFilePath,header=None,names=['fd'],delim_whitespace=True)
                    df['meanFD']=fd.mean()[0]
                    df['fd']=fd

                    if len(GroupDF)==0:
                        GroupDF=df
                    else:
                        GroupDF=pd.concat((GroupDF,df),ignore_index=True)
                except:
                    print 'No DR .txt file found or error'

    GroupDF.reset_index(inplace=True)

    motionInfo=GroupDF.groupby(['Subject_ID','FB']).mean()['meanFD']
    if saveMotionInfo:
        motionInfo.to_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_meanFD.csv')

    return GroupDF,motionInfo

def getSubjectList(GroupDF,RejectMotion=True,motionThresh=0.2):

    #reject large motion subjects
    allsubj=unique(GroupDF['Subject_ID'])
    motionReject=unique((GroupDF[GroupDF.meanFD>motionThresh]['Subject_ID']))
    if RejectMotion:
        goodsubj=np.setdiff1d(allsubj,motionReject)
    else:
        goodsubj=allsubj

    return goodsubj

def createTimeSeriesPlots(GroupDF,goodsubj,DMN_name='RSN3',title='DMN_Activity',ylabel=''):

    sns.set_context("paper")
    #plt.subplots(2,1,figsize=(12, 6))
    f, axarr = plt.subplots(1, sharex=True,figsize=(18, 9))
    sns.set(style="white")
    dmnPlot=sns.tsplot(data=GroupDF[GroupDF.Subject_ID.isin(goodsubj)],time='TR',unit='Subject',condition='FB',value=DMN_name,ci=68)
    #get ideal DMN time line
    dmnIdeal=pd.read_csv('/home/jmuraskin/Projects/NFB/analysis/DMN_ideal_2.csv')
    dmnPlot.plot((dmnIdeal['Wander']-dmnIdeal['Focus'])/(3*max(dmnIdeal['Wander'])),'k--')
    #dmnPlot.plot(dmnIdeal['Focus'][4:]/(3*max(dmnIdeal['Focus'])),'r--')
    # dmnPlot.set_ylim([-.8,.8])
    dmnPlot.set_ylabel(ylabel,{'fontsize':18})
    dmnPlot.set_xlabel('')
    dmnPlot.set_title(title,{'fontsize':24})
    f.savefig('%s_timeseries.pdf' % DMN_name, dpi=600)
