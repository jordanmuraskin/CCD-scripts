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
CCD_numbers=[12,14,15,16,17,18,20,21,22,27,32,33,34,40,42,51,52,59,60,61,63,64,65,66,67,73,74,75,76,81,83,85,88,89,91,93,95,97,98,99]

# Specify the subject directories

# subject_list = ['CCD060','CCD066','CCD089']
# subject_list = ['CCD015','CCD015','CCD017','CCD066','CCD089','CCD052','CCD076','CCD059','CCD064','CCD083']
subject_list=[]
for ccd in CCD_numbers:
    subject_list.append('CCD0%s' % ccd)

# scan_order1=list(SubjInfo.loc[subject_list]['V1_NSI_001'])
# scan_order2=list(SubjInfo.loc[subject_list]['V1_NSI_005'])






# Specify the location of the data.
data_dir = os.path.abspath('/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1')




for feedbackRun in range(2):

    workflow = pe.Workflow(name= "level1")
    workflow.base_dir = os.path.abspath('/home/jmuraskin/Projects/CCD/')
    workflow.config = {"execution": {"crashdump_dir":os.path.abspath('%s/crashdumps' % workflow.base_dir)}}
    working_dir = os.path.abspath('%s/working_v1' % workflow.base_dir)



    # Workflow base directory
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    workflow = pe.Workflow(name='feedback_run-%d' % feedbackRun, base_dir=working_dir)


    # Map field names to individual subject runs.
    info = dict(func=[['subject_id', ['functional_mni_other_resolutions_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_apply_isoxfm_3.0/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_fwhm_6/residual_antswarp_maths' % (feedbackRun+1)]]],
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
    #
    # def createOperandFileName(infoDict):
    #     print infoDict[0]
    #     print infoDict[1]
    #     filename= '%s_data_/%s.nii.gz' % (infoDict[0],infoDict[1])
    #     return filename

    # add mean image to fmri
    addMeanImage =  pe.MapNode(interface=fsl.maths.MultiImageMaths(),name='addMeanImage',iterfield=['in_file'])
    addMeanImage.inputs.op_string = "-add %s"
    # addMeanImage.inputs.operand_files = ['%s_data_/%s.nii.gz']
    # addMeanImage.inputs.out_file =
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


    def subjectinfo(subject_id,r):
        from pandas import read_csv
        from nipype.interfaces.base import Bunch

        Order1_onsets = [22,150,248,406,564,662]
        Order1_durations = [30,60,90,60,30,90]
        Order2_onsets = [56,214,342,470,598,756]
        Order2_durations = [90,30,60,90,60,30]
        # button_onsets = [52,144,206,238,330,392,454,546,578,640,732]
        # button_duration =[2]
        #Make subject specific EVs given feedback ordering
        output=[]
        names=['Focus','Wander']
        SubjInfo = read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/NARSAD_stimulus_JM.csv')
        SubjInfo.set_index('JM_INTERNAL',inplace=True)
        if r==0:
            paradigmType=SubjInfo.loc[subject_id]['SCAN_1_PARADIGM']
        else:
            paradigmType=SubjInfo.loc[subject_id]['SCAN_2_PARADIGM']
        if paradigmType==0 or paradigmType == 2:
            focus_onset = Order2_onsets
            focus_durations = Order2_durations
            wander_onset = Order1_onsets
            wander_durations = Order1_durations
        elif paradigmType==1 or paradigmType == 3:
            focus_onset = Order1_onsets
            focus_durations = Order1_durations
            wander_onset = Order2_onsets
            wander_durations = Order2_durations
        output.insert(r,Bunch(conditions=names,onsets=[focus_onset, wander_onset],
                              durations=[focus_durations, wander_durations], regressors=None))
        return output
    ## end moral dilemma


    # workflow.connect(datasource, 'func', modelspec,'functional_runs')
    # workflow.connect(smooth, 'out_file', modelspec,'functional_runs')


    #modelfit
    modelfit = create_modelfit_workflow(name='feedback')
    modelfit.inputs.inputspec.interscan_interval = TR
    modelfit.inputs.inputspec.model_serial_correlations = True
    modelfit.inputs.inputspec.bases = {'dgamma': {'derivs': True}}
    cont1 = ['Focus>Wander','T', ['Focus','Wander'],[1,-1]]
    cont2 = ['Wander>Focus','T', ['Focus', 'Wander'],[-1,1]]
    cont3 = ['Mean Focus','T',['Focus'],[1]]
    cont4 = ['Mean Wander','T',['Wander'],[1]]
    cont5 = ['Average Activation', 'T', ['Focus', 'Wander'],[.5,.5]]
    # cont3 = ['Task','F', [cont1,cont2]]

    modelfit.inputs.inputspec.contrasts = [cont1, cont2, cont3,cont4,cont5]

    workflow.connect([(infosource,modelspec,[(('subject_id',subjectinfo,feedbackRun),'subject_info')])])

    workflow.connect(modelspec, 'session_info', modelfit, 'inputspec.session_info')

    #workflow.connect(datasource, 'func', modelfit, 'inputspec.functional_data')
    workflow.connect(addMeanImage,'out_file', modelfit, 'inputspec.functional_data')


    workflow.run(plugin='MultiProc')
