import pandas as pd
import numpy as np
from numpy import unique
from scipy.stats import zscore,spearmanr,pearsonr
import seaborn as sns
import matplotlib.pylab as plt
import os.path
from scipy.signal import butter,filtfilt
import matplotlib as mpl
from matplotlib import cm
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import networkx as nx
from sklearn import linear_model
from sklearn import cross_validation
from sklearn import metrics
from scipy.stats import ttest_1samp
from mne.stats.multi_comp import fdr_correction
import plotly.plotly as py
from plotly.graph_objs import *
from nilearn import plotting
from nilearn import image


class MplColorHelper:

  def __init__(self, cmap_name, start_val, stop_val):
    self.cmap_name = cmap_name
    self.cmap = plt.get_cmap(cmap_name)
    self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
    self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

  def get_rgb(self, val):
    return self.scalarMap.to_rgba(val)

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

    motionInfo=GroupDF.groupby(['Subject_ID','FB','scanorder']).mean()['meanFD']
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

def createTimeSeriesPlots(GroupDF,goodsubj,DMN_name='RSN3',title='DMN_Activity',ylabel='',figsize=(18,9),savefig=True):

    sns.set_context("paper")
    #plt.subplots(2,1,figsize=(12, 6))
    f, axarr = plt.subplots(1, sharex=True,figsize=figsize)
    sns.set(style="white")
    dmnPlot=sns.tsplot(data=GroupDF[GroupDF.Subject_ID.isin(goodsubj)],time='TR',unit='Subject',condition='FB',value=DMN_name,ci=68)
    #get ideal DMN time line
    dmnIdeal=pd.read_csv('/home/jmuraskin/Projects/NFB/analysis/DMN_ideal_2.csv')
    dmnPlot.plot((dmnIdeal['Wander']-dmnIdeal['Focus'])/(3*max(dmnIdeal['Wander'])),'k--')
    #dmnPlot.plot(dmnIdeal['Focus'][4:]/(3*max(dmnIdeal['Focus'])),'r--')
    # dmnPlot.set_ylim([-.8,.8])
    dmnPlot.set_ylabel(ylabel,{'fontsize':18})
    dmnPlot.set_xlabel('TR')
    dmnPlot.set_title(title,{'fontsize':24})
    if savefig:
        f.savefig('%s_timeseries.pdf' % DMN_name, dpi=600)


def createSubjectModelBarPlot(GroupDF,goodsubj,figsize=(18,9),savefig=True):

    f, axarr = plt.subplots(1, sharex=True,figsize=figsize)
    sns.set(style="white")

    maxModel=GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['Subject'])['modelcorr'].max().sort_values(ascending=False)
    sortedOrder=maxModel.index

    sns.barplot(data=GroupDF[GroupDF.Subject_ID.isin(goodsubj)],x='Subject',y='modelcorr',hue='FB',order=sortedOrder)

    if savefig:
        f.savefig('Subject_ModelCorrelations.pdf', dpi=600)

def createScanOrderBarPlot(GroupDF,goodsubj,BV=False,savefig=True):
    plt.figure()
    if BV:
        sns.factorplot(data=GroupDF[GroupDF.Subject_ID.isin(goodsubj)],x='FB',y='modelcorr',hue='scanorder',kind='bar',units='Subject',ci=68)
    else:
        sns.violinplot(data=GroupDF[GroupDF.Subject_ID.isin(goodsubj)],x='FB',y='modelcorr',hue='scanorder',split='True',bw=.4,inner='quartile')
    if savefig:
        plt.savefig('ScanOrder_ModelCorrelations.pdf',dpi=600)

def printModelCorrelations(GroupDF,goodsubj,DMN_name='RSN3'):
    dmnIdeal=pd.read_csv('/home/jmuraskin/Projects/NFB/analysis/DMN_ideal_2.csv')
    print 'No Feedback Focus Correlation= %0.2f' % GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['FB','TR']).mean()[DMN_name].loc['NOFEEDBACK'].corr(dmnIdeal['Focus'])
    print 'Feedback Focus Correlation= %0.2f' % GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['FB','TR']).mean()[DMN_name].loc['FEEDBACK'].corr(dmnIdeal['Focus'])
    print 'No Feedback Wander Correlation= %0.2f' % GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['FB','TR']).mean()[DMN_name].loc['NOFEEDBACK'].corr(dmnIdeal['Wander'])
    print 'Feedback Wander Correlation= %0.2f' % GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['FB','TR']).mean()[DMN_name].loc['FEEDBACK'].corr(dmnIdeal['Wander'])
    print 'No Feedback Overall Correlation= %0.2f' % GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['FB','TR']).mean()[DMN_name].loc['NOFEEDBACK'].corr(dmnIdeal['Wander']-dmnIdeal['Focus'])
    print 'Feedback Overall Correlation= %0.2f' % GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['FB','TR']).mean()[DMN_name].loc['FEEDBACK'].corr(dmnIdeal['Wander']-dmnIdeal['Focus'])


def generateHeatMaps(GroupDF,goodsubj):

    numberOfICs=10
    columnNames=[]
    for rsnNumber in range(numberOfICs):
            columnNames.append('RSN%d' % rsnNumber)

    heatmapDF=GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['Subject_ID','FB','TR']).mean()
    hmDiff=np.zeros((10,10,len(unique(GroupDF[GroupDF.Subject_ID.isin(goodsubj)]['Subject_ID']))))
    hmFB=hmDiff.copy()
    hmNFB=hmDiff.copy()

    for indx,subj in enumerate(unique(GroupDF[GroupDF.Subject_ID.isin(goodsubj)]['Subject_ID'])):
        hmFB[:,:,indx]=heatmapDF.loc[subj,'FEEDBACK'][columnNames].corr()
        hmNFB[:,:,indx]=heatmapDF.loc[subj,'NOFEEDBACK'][columnNames].corr()
        hmDiff[:,:,indx]=np.arctan(heatmapDF.loc[subj,'FEEDBACK'][columnNames].corr())*np.sqrt(405)-np.arctan(heatmapDF.loc[subj,'NOFEEDBACK'][columnNames].corr())*np.sqrt(405)

    return hmFB,hmNFB,hmDiff



def dist (A,B):
        return np.linalg.norm(np.array(A)-np.array(B))

def get_idx_interv(d, D):
    k=0
    while(d>D[k]):
        k+=1
    return  k-1

class InvalidInputError(Exception):
    pass

def deCasteljau(b,t):
    N=len(b)
    if(N<2):
        raise InvalidInputError("The  control polygon must have at least two points")
    a=np.copy(b) #shallow copy of the list of control points
    for r in range(1,N):
        a[:N-r,:]=(1-t)*a[:N-r,:]+t*a[1:N-r+1,:]
    return a[0,:]

def BezierCv(b, nr=5):
    t=np.linspace(0, 1, nr)
    return np.array([deCasteljau(b, t[k]) for k in range(nr)])

def makeChordDiagram(G,cmap='coolwarm',plotName='ChordDiagram',scale=[-1.0,1.0],title='',savefig=True):

    widthScale=50.0/max(scale)
    COL = MplColorHelper(cmap, scale[0], scale[1])
    #GET LABEL NAMES IF THERE ARE ANY
    labels=G.node.keys()
    Edges=G.edge
    #Get all edge weights
    E=[]
    Weights=[]
    for indx1,j in enumerate(Edges):
        for indx2,k in enumerate(Edges[j]):
            E.append([j,k])
            if j==k:
                Weights.append(0)
            else:
                Weights.append(Edges[j][k]['weight'])
    layt=nx.drawing.layout.circular_layout(G)
    L=len(layt)

    dist(layt[0], layt[5])

    Dist=[0, dist([1,0], 2*[np.sqrt(2)/2]), np.sqrt(2),
    dist([1,0],  [-np.sqrt(2)/2, np.sqrt(2)/2]), 2.0]
    params=[1.2, 1.5, 1.8, 2.1]

    #set node color
    minColor=np.array(COL.get_rgb(scale[0]))*255.0
    maxColor=np.array(COL.get_rgb(scale[1]))*255.0
    node_color=['rgba(%f,%f,%f,1)' % (minColor[0],minColor[1],minColor[2]),
                'rgba(%f,%f,%f,1)' % (maxColor[0],maxColor[1],maxColor[2])]*5

    #Get Node Positions
    Xn=[layt[k][0] for k in range(L)]
    Yn=[layt[k][1] for k in range(L)]


    lines=[]# the list of dicts defining   edge  Plotly attributes
    edge_info=[]# the list of points on edges where  the information is placed

    for j, e in enumerate(E):
        A=np.array(layt[e[0]])
        B=np.array(layt[e[1]])
        d=dist(A, B)
        K=get_idx_interv(d, Dist)
        b=[A, A/params[K], B/params[K], B]
        # color=edge_colors[0]
        pts=BezierCv(b, nr=5)
        mark=list(deCasteljau(b,0.9))
        rgb=np.array(COL.get_rgb(Weights[j]))*255.0

        lines.append(Scatter(x=list(pts[:,0]),
                             y=list(pts[:,1]),
                             mode='lines',
                             line=Line(
                                      shape='spline',
                                      color='rgba(%d,%d,%d,.9)' % (rgb[0],rgb[1],rgb[2]),
                                      width=abs(Weights[j])*widthScale#The  width is proportional to the edge weight
                                     ),
                            hoverinfo='none'
                           )
                    )




    trace2=Scatter(x=Xn,
               y=Yn,
               mode='markers',
               name='',
               marker=Marker(symbol='dot',
                             size=.5,
                             color=node_color,
                             cmin=scale[0],cmax=scale[1],
                             colorscale='coolwarm',
                             colorbar = dict(tickangle=20,thickness=15,x=1.2)
                             ),
               text=labels,
               hoverinfo='text',

               )

    axis=dict(showline=False, # hide axis line, grid, ticklabels and  title
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title=''
              )
    width=600
    height=575

    layout=Layout(title=title,
                  paper_bgcolor='rgba(0,0,0,0)',
                  plot_bgcolor='rgba(0,0,0,0)',
                  font= Font(size=12),
                  showlegend=False,
                  autosize=False,
                  width=width,
                  height=height,
                  xaxis=XAxis(axis),
                  yaxis=YAxis(axis),
                  margin=Margin(l=120,
                                r=120,
                                b=92,
                                t=135,
                              ),
                  hovermode='closest',
                  images=[dict(
                        source="./RSN/RSN0.png",
                        xref="paper", yref="paper",
                        x=.93, y=.4,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN1.png",
                        xref="paper", yref="paper",
                        x=.85, y=.7,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN2.png",
                        xref="paper", yref="paper",
                        x=.55, y=.94,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN3.png",
                        xref="paper", yref="paper",
                        x=.25, y=.94,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN4.png",
                        xref="paper", yref="paper",
                        x=-0.02, y=.7,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN5.png",
                        xref="paper", yref="paper",
                        x=-.1, y=.4,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN6.png",
                        xref="paper", yref="paper",
                        x=-0.03, y=.1,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN7.png",
                        xref="paper", yref="paper",
                        x=.25, y=-0.14,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN8.png",
                        xref="paper", yref="paper",
                        x=.55, y=-0.14,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="./RSN/RSN9.png",
                        xref="paper", yref="paper",
                        x=.85, y=.1,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      )]
                  )

    data=Data(lines+edge_info+[trace2])
    fig=Figure(data=data, layout=layout)
    if savefig:
        py.image.save_as(fig, filename='%s.png' % plotName)
    return fig


def heatmap2Chord(matrix,plotName='ChordDiagram',title='',savefig=True,scale=[-1,1]):

    # check size of matrix
    matSize=np.shape(matrix)
    if len(matSize)==3:
        #get mean of matrix
        matrix=np.mean(matrix,axis=2)

    G = nx.from_numpy_matrix(matrix)

    fig=makeChordDiagram(G,cmap='coolwarm',plotName='ChordDiagram',scale=scale,title=title,savefig=savefig)

    return fig



def LinRegression(X,y):
    regr = linear_model.LinearRegression()
    regr.fit(X,y)

    score=regr.score(X,y)
    resids=y-regr.predict(X)

    return score,resids,regr.coef_,regr.predict(X)


def leaveOneOutCV(clf,X,y,LOO=False,numFolds=10):
    from sklearn.cross_validation import LeaveOneOut,KFold
    coefs=np.zeros((X.shape[1],))
    intercept=0.0
    predicted=np.zeros((len(y),))
    if LOO:
        loo = LeaveOneOut(n=len(y))
        numFolds=len(y)
    else:
        loo = KFold(n=len(y),n_folds=numFolds)
        # numFolds=10
    for train_index, test_index in loo:
        clf.fit(X[train_index,:],y[train_index])
        predicted[test_index]=clf.predict(X[test_index,:])
        intercept+=clf.intercept_
        coefs+=clf.coef_

    intercept=intercept/numFolds
    coefs=coefs/numFolds
    return predicted,intercept,coefs

def bayesianRidge(X,y):
    clf = linear_model.BayesianRidge(compute_score=True)
    predicted = cross_validation.cross_val_predict(clf, X,y,cv=408)
    return clf,predicted

def GroupRegression(GroupDF,goodsubj,feedback,numFolds=10):

    numberOfICs=10
    columnNames=[]
    for rsnNumber in range(numberOfICs):
            columnNames.append('RSN%d' % rsnNumber)
    dmnIdeal=pd.read_csv('/home/jmuraskin/Projects/NFB/analysis/DMN_ideal_2.csv')

    SubjectDF = GroupDF[GroupDF.Subject_ID.isin(goodsubj)].groupby(['Subject_ID','FB','TR']).mean()
    clf = linear_model.LinearRegression()

    for indx,subj in enumerate(unique(GroupDF['Subject_ID'])):
        predicted,intercepts,coef = leaveOneOutCV(clf,np.array(SubjectDF.loc[subj,feedback][columnNames]),dmnIdeal['Wander']-dmnIdeal['Focus'],numFolds=numFolds)
        if indx==0:
            groupGLM=pd.DataFrame({'TR':range(408),'predicted':predicted,'subj':[subj]*408})
            coefs=pd.DataFrame({'Coef':coef,'pe':range(10),'subj':[subj]*10})
            performance=pd.DataFrame({'R':[pearsonr(dmnIdeal['Wander']-dmnIdeal['Focus'],predicted)[0]],'subj':[subj]})
        else:
            df=pd.DataFrame({'TR':range(408),'predicted':predicted,'subj':[subj]*408})
            groupGLM=pd.concat((groupGLM,df),ignore_index=True)
            coefs=pd.concat((coefs,pd.DataFrame({'Coef':coef,'pe':range(10),'subj':[subj]*10})),ignore_index=True)
            performance=pd.concat((performance,pd.DataFrame({'R':[pearsonr(dmnIdeal['Wander']-dmnIdeal['Focus'],predicted)[0]],'subj':[subj]})),ignore_index=True)

    return groupGLM,coefs,performance

def linearRegressionData(GroupDF,goodsubj,numFolds=10):
    print 'Running Feedback on Regressions'
    fb_pred,fb_coefs,fb_performance=GroupRegression(GroupDF[GroupDF.Subject_ID.isin(goodsubj)],goodsubj,'FEEDBACK',numFolds=numFolds)
    print 'Finished...'
    print 'Running Feedback off Regressions'
    nfb_pred,nfb_coefs,nfb_performance=GroupRegression(GroupDF[GroupDF.Subject_ID.isin(goodsubj)],goodsubj,'NOFEEDBACK',numFolds=numFolds)
    print 'Finished...'

    fb_pred['fb']='FEEDBACK'
    nfb_pred['fb']='NOFEEDBACK'
    predictions=pd.concat((fb_pred,nfb_pred),ignore_index=True)

    fb_coefs['fb']='FEEDBACK'
    nfb_coefs['fb']='NOFEEDBACK'
    coefs=pd.concat((fb_coefs,nfb_coefs),ignore_index=True)

    fb_performance['fb']='FEEDBACK'
    nfb_performance['fb']='NOFEEDBACK'
    performance=pd.concat((fb_performance,nfb_performance),ignore_index=True)

    return predictions,coefs,performance,fb_coefs,nfb_coefs


def createRegressionPlots(predictions,performance,coefs,fb_coefs,nfb_coefs,GroupDF,goodsubj,savefig=True):
    f=plt.figure(figsize=(22,12))
    ax1=plt.subplot2grid((2,4),(0,0), colspan=3)
    ax2=plt.subplot2grid((2,4),(0,3))
    ax3=plt.subplot2grid((2,4),(1,0), colspan=2)
    ax4=plt.subplot2grid((2,4),(1,2), colspan=2)

    dmnIdeal=pd.read_csv('/home/jmuraskin/Projects/NFB/analysis/DMN_ideal_2.csv')

    sns.tsplot(data=predictions,time='TR',value='predicted',unit='subj',condition='fb',ax=ax1)
    ax1.plot((dmnIdeal['Wander']-dmnIdeal['Focus'])/3,'k--')
    ax1.set_title('Average Predicted Time Series')

    g=sns.violinplot(data=performance,x='fb',y='R',split=True,bw=.3,inner='quartile',ax=ax2)
    # plt.close(g.fig)
    ax2.set_title('Mean Subject Time Series Correlations')
    g=sns.violinplot(data=coefs,x='pe',y='Coef',hue='fb',split=True,bw=.3,inner='quartile',ax=ax3)
    g.plot([-1,10],[0,0],'k--')
    g.set_xlim([-.5,9.5])


    t,p=ttest_1samp(np.array(fb_coefs['Coef']-nfb_coefs['Coef']).reshape(len(unique(GroupDF[GroupDF.Subject_ID.isin(goodsubj)]['Subject_ID'])),10),0)
    p05,padj=fdr_correction(p,0.05)

    sns.barplot(x=range(len(t)),y=t,ax=ax4,color='Red')
    for idx,pFDR in enumerate(p05):
        if pFDR:
            ax4.scatter(idx,t[idx]+ np.sign(t[idx])*0.2,marker='*',s=75)
    ax4.set_xlim([-0.5,10])
    ax4.set_xlabel('pe')
    ax4.set_ylabel('t-value')
    ax4.set_title('FB vs. nFB PE')

    for ax in [ax1,ax2,ax3,ax4]:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(18)
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(12)

    f.tight_layout()
    if savefig:
        f.savefig('RSN_LinearRegPrediction.pdf',dpi=300)


def createTFCEfMRIOverlayImages(TFCEposImg,posImg,TFCEnegImg,negImg,title='',vmax=8,display_mode='z',slices=range(-20,50,10),threshold=0.94999,plotToAxis=False,f=[],axes=[],colorbar=True,tight_layout=False):

    bg_img='./Templates/MNI152_.5mm_masked_edged.nii.gz'
    # threshold=0.949
    pos=image.math_img("np.multiply(img1,img2)",
                         img1=image.threshold_img(TFCEposImg,threshold=threshold),img2=posImg)
    neg=image.math_img("np.multiply(img1,img2)",
                         img1=image.threshold_img(TFCEnegImg,threshold=threshold),img2=negImg)
    fw=image.math_img("img1-img2",img1=pos,img2=neg)

    if plotToAxis:
        display=plotting.plot_stat_map(fw,display_mode=display_mode,threshold=0,
                                       cut_coords=slices,vmax=vmax,colorbar=colorbar,bg_img=bg_img,black_bg=False,title=title,dim=0,figure=f,axes=axes)
    else:
        display=plotting.plot_stat_map(fw,display_mode=display_mode,threshold=0,
        cut_coords=slices,vmax=vmax,colorbar=colorbar,bg_img=bg_img,
        black_bg=False,title=title,dim=0)
    if tight_layout:
        display.tight_layout()
    return display


def runRLMR(y,X,modelNames=[],RLM=True,addconstant=True,plotFigure=True):
    import statsmodels.api as sm
    if addconstant:
        X=sm.add_constant(X)
    if RLM:
        model = sm.RLM(y, X)
    else:
        model = sm.OLS(y,X)
    results = model.fit()
    print results.summary()

    if plotFigure:
        #first figure out how many plots
        numX=X.shape[1]-1
        if numX>2:
            fig, axarr = plt.subplots(int(np.ceil(numX/3.0)),3,figsize=(20,20))
        else:
            fig, axarr = plt.subplots(1,2,figsize=(20,20))
        row=0
        column=0
        for n in range(1,numX+1):
#             fig, axarr = plt.subplots(int(np.ceil(numX/3.0)),3,figsize=(10,10))
            if numX>2:
                sm.graphics.plot_ccpr(results,results.model.exog_names[n], ax = axarr[row][column])
            else:
                sm.graphics.plot_ccpr(results,results.model.exog_names[n], ax = axarr[column])

            axarr[row][column].set_title(modelNames[n-1])
            column+=1
            if column==3:
                column=0
                row+=1
    return results
