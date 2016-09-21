import os                                    # system functions
import yaml
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.fsl as fsl          # fsl
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.algorithms.modelgen as model   # model generation
from CPAC.registration import create_wf_apply_ants_warp
from nipype.workflows.fmri.fsl import (create_featreg_preproc,
                                       create_modelfit_workflow,
                                       create_reg_workflow)





template = '/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain.nii.gz'



# CCD_numbers=[15,17,18,21,23,33,40,52,59,64,66,74,76,83,89,95]
CCD_numbers=[12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,31,32,33,34,40,41,42,51,52,
53,59,60,61,62,63,64,65,66,67,71,72,73,74,75,76,80,81,82,83,84,85,86,87,88,89,
90,91,92,93,94,95,96,97,98,99]
# CCD_numbers=[16]
# Specify the subject directories

# subject_list = ['CCD060','CCD066','CCD089']
# subject_list = ['CCD015','CCD015','CCD017','CCD066','CCD089','CCD052','CCD076','CCD059','CCD064','CCD083']
subject_list=[]
for ccd in CCD_numbers:
    subject_list.append('CCD0%s' % ccd)

# scan_order1=list(SubjInfo.loc[subject_list]['V1_NSI_001'])
# scan_order2=list(SubjInfo.loc[subject_list]['V1_NSI_005'])


globalSR=0
if globalSR:
    GSR='GSR1'
else:
    GSR='GSR0'



# Specify the location of the data.
data_dir = os.path.abspath('/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1')


ROI_file='/home/jmuraskin/Projects/CCD/working_v1/ROIs/R_AI.nii.gz'

for feedbackRun in range(2):

    workflow = pe.Workflow(name= "ROI_ts")
    workflow.base_dir = os.path.abspath('/home/jmuraskin/Projects/CCD/')
    workflow.config = {"execution": {"crashdump_dir":os.path.abspath('%s/crashdumps' % workflow.base_dir)}}
    working_dir = os.path.abspath('%s/working_v1/SMG_ROI_ts_%s' % (workflow.base_dir,GSR))



    # Workflow base directory
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    workflow = pe.Workflow(name='feedback_run-%d' % feedbackRun, base_dir=working_dir)


    # Map field names to individual subject runs.
    info = dict(func=[['subject_id', ['functional_mni_other_resolutions_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_apply_isoxfm_3.0/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global%d.motion1.quadratic1.gm0.compcor1.csf1/_fwhm_6/residual_antswarp_maths' % (feedbackRun+1,globalSR)]]],
    funcMean=[['subject_id',['mean_functional_in_mni/_scan_feedback_%d/fb_%d_calc_tshift_resample_volreg_calc_tstat_antswarp' % (feedbackRun+1,feedbackRun+1)]]])

    infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name="infosource")
    infosource.iterables = ('subject_id', subject_list)

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['func','funcMean']),
                         name = 'datasource')
    datasource.inputs.base_directory = data_dir
    datasource.inputs.template = '%s_data_/%s.nii.gz'
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True
    workflow.connect(infosource, 'subject_id', datasource, 'subject_id')


    meanTs=pe.Node(interface=fsl.utils.ImageMeants(),name='ExtractTimeSeries',iterfield=['in_file'])
    meanTs.inputs.mask=ROI_file

    workflow.connect([(datasource,meanTs,[('func','in_file')])])

    workflow.run(plugin='MultiProc',plugin_args={'n_procs':15})

    addMeanImage =  pe.MapNode(interface=fsl.maths.MultiImageMaths(),name='addMeanImage',iterfield=['in_file'])
    addMeanImage.inputs.op_string = "-add %s"
    workflow.connect([(datasource,addMeanImage,[('func','in_file')]),
        (datasource,addMeanImage,[('funcMean','operand_files')])])
    #modelspec
    TR = 2
    ##  moral dilemma
    modelspec = pe.Node(interface=model.SpecifyModel(),name="modelspec")
    modelspec.inputs.input_units = 'secs'
    modelspec.inputs.time_repetition = TR
    modelspec.inputs.high_pass_filter_cutoff = 100

    workflow.connect(addMeanImage, 'out_file', modelspec,'functional_runs')

    def subjectinfo(subject_id,r,working_dir):
        from pandas import read_csv
        from nipype.interfaces.base import Bunch
        from scipy.stats import zscore

        # drFileLocation='/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1'
        # # numberOfICs=10
        # # columnNames=[]
        # # for rsnNumber in range(numberOfICs):
        # #     columnNames.append('RSN%d' % rsnNumber)
        filterOn=False
        zscoreOn=True
        lowpass=0.1
        globalNR=0
        #load DMN_network
        ROIFilePath = '%s/feedback_run-%d/_subject_id_%s/ExtractTimeSeries/residual_antswarp_maths_ts.txt' % (working_dir,r,subject_id)
        df = read_csv(ROIFilePath,header=None,names=['ROI'],delim_whitespace=True)
        df['Subject_ID'] = subject_id
        # df['Subject'] = indx
        df.index.name = 'TR'
        df.reset_index(level=0,inplace=True)


        # regressors=read_csv('./PPI.csv',sep=',',header=0)
        # regressor_names=list(regressors.keys().values)
        # regressor_values=list(regressors.values.transpose().tolist())

        #Make subject specific EVs given feedback ordering
        output=[]
        SubjInfo = read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/NARSAD_stimulus_JM.csv')
        SubjInfo.set_index('JM_INTERNAL',inplace=True)
        if r==0:
            paradigmType=SubjInfo.loc[subject_id]['SCAN_1_PARADIGM']
        else:
            paradigmType=SubjInfo.loc[subject_id]['SCAN_2_PARADIGM']
        if paradigmType==0 or paradigmType == 2:
            signFlip=1
        elif paradigmType==1 or paradigmType == 3:
            signFlip=-1
        regressors=read_csv('./PPI.csv',sep=',',header=0)
        regressors['Cont']=regressors['Cont']*signFlip
        regressors['Cont_Deriv']=regressors['Cont_Deriv']*signFlip
        regressor_names=list(regressors.keys().values)
        regressor_values=list(regressors.values.transpose().tolist())
        PPI=list(signFlip*regressors['Cont']*df['ROI'])
        regressor_names+=['PHYS','PPI']
        regressor_values.append(list(df['ROI']))
        regressor_values.append(PPI)
        # regressor_values.append(list(signFlip*df['RSN3']))
        output.insert(0,Bunch(regressor_names=regressor_names,regressors=regressor_values))

        return output

    modelfit = create_modelfit_workflow(name='ROI_model_fit')
    modelfit.inputs.inputspec.interscan_interval = TR
    modelfit.inputs.inputspec.model_serial_correlations = True
    modelfit.inputs.inputspec.bases = {'dgamma': {'derivs': False}}
    cont1 = ['ROI Correlation','T', ['PPI'],[1]]

    modelfit.inputs.inputspec.contrasts = [cont1]

    workflow.connect([(infosource,modelspec,[(('subject_id',subjectinfo,feedbackRun,working_dir),'subject_info')])])

    workflow.connect(modelspec, 'session_info', modelfit, 'inputspec.session_info')

    #workflow.connect(datasource, 'func', modelfit, 'inputspec.functional_data')
    workflow.connect(addMeanImage,'out_file', modelfit, 'inputspec.functional_data')


    workflow.run(plugin='MultiProc',plugin_args={'n_procs':20})


    # def createOperandFileName(infoDict):
    #     print infoDict[0]
    #     print infoDict[1]
    #     filename= '%s_data_/%s.nii.gz' % (infoDict[0],infoDict[1])
    #     return filename

    # # add mean image to fmri
    # addMeanImage =  pe.MapNode(interface=fsl.maths.MultiImageMaths(),name='addMeanImage',iterfield=['in_file'])
    # addMeanImage.inputs.op_string = "-add %s"
    # workflow.connect([(datasource,addMeanImage,[('func','in_file')]),
    #     (datasource,addMeanImage,[('funcMean','operand_files')])])
