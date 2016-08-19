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

def createScanOrderBarPlot(GroupDF,goodsubj,savefig=True):
    plt.figure()
    sns.factorplot(data=GroupDF[GroupDF.Subject_ID.isin(goodsubj)],x='FB',y='modelcorr',hue='scanorder',kind='bar',units='Subject',ci=68)
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
    import plotly.plotly as py
    from plotly.graph_objs import Data,Layout,Figure
    py.sign_in('jordan.muraskin', 'q8p1xy3l62')


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
    minColor=array(COL.get_rgb(scale[0]))*255.0
    maxColor=array(COL.get_rgb(scale[1]))*255.0
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
        color=edge_colors[0]
        pts=BezierCv(b, nr=5)
        mark=list(deCasteljau(b,0.9))
        rgb=array(COL.get_rgb(Weights[j]))*255.0

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
    width=1200
    height=1150

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
                  margin=Margin(l=240,
                                r=240,
                                b=185,
                                t=270,
                              ),
                  hovermode='closest',
                  images=[dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN0.png",
                        xref="paper", yref="paper",
                        x=.93, y=.4,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN1.png",
                        xref="paper", yref="paper",
                        x=.85, y=.7,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN2.png",
                        xref="paper", yref="paper",
                        x=.55, y=.94,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN3.png",
                        xref="paper", yref="paper",
                        x=.25, y=.94,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN4.png",
                        xref="paper", yref="paper",
                        x=-0.02, y=.7,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN5.png",
                        xref="paper", yref="paper",
                        x=-.1, y=.4,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN6.png",
                        xref="paper", yref="paper",
                        x=-0.03, y=.1,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN7.png",
                        xref="paper", yref="paper",
                        x=.25, y=-0.14,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN8.png",
                        xref="paper", yref="paper",
                        x=.55, y=-0.14,
                        sizex=.2, sizey=.2,
                        xanchor="left", yanchor="bottom",
                        layer='below'

                      ),
                    dict(
                        source="https://dl.dropboxusercontent.com/u/13758/RSN/RSN9.png",
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
    return iplot(fig)


def heatmap2Chord(matrix,plotName='ChordDiagram',title='',savefig=True,scale=[-1,1]):
    import networkx as nx
    # check size of matrix
    matSize=np.shape(matrix)
    if len(matSize)==3:
        #get mean of matrix
        matrix=np.mean(matrix,axis=2)

    G = nx.from_numpy_matrix(matrix)

    plot_out=makeChordDiagram(G,cmap='coolwarm',plotName='ChordDiagram',scale=scale,title=title,savefig=True)

    return plot_out
