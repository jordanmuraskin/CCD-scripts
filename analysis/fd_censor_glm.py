import os                                    # system functions
import yaml
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.fsl as fsl          # fsl
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util




template = '/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain.nii.gz'



# CCD_numbers=[15,17,18,21,23,33,40,52,59,64,66,74,76,83,89,95]
CCD_numbers=[12]
# ,14,15,16,17,18,19,20,21,22,23,24,25,26,27,31,32,33,34,40,41,42,51,52,
# 53,59,60,61,62,63,64,65,66,67,71,72,73,74,75,76,80,81,82,83,84,85,86,87,88,89,
# 90,91,92,93,94,95,96,97,98,99]
# # CCD_numbers=[16]
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
csv_dir=os.path.abspath('/home/jmuraskin/Projects/CCD/working_v1/censor_directory')




for feedbackRun in range(2):

    workflow = pe.Workflow(name= "censorAndRerunDR")
    workflow.base_dir = os.path.abspath('/home/jmuraskin/Projects/CCD/')
    workflow.config = {"execution": {"crashdump_dir":os.path.abspath('%s/crashdumps' % workflow.base_dir)}}
    working_dir = os.path.abspath('%s/working_v1/censorAndRerunDR' % workflow.base_dir)



    # Workflow base directory
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    workflow = pe.Workflow(name='feedback_run-%d' % feedbackRun, base_dir=working_dir)

    # Map field names to individual subject runs.
    info = dict(func=[['subject_id', ['functional_mni_other_resolutions_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_apply_isoxfm_3.0/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf1/_fwhm_6/residual_antswarp_maths' % (feedbackRun+1)]]],
    mask=[['subject_id',['functional_brain_mask_to_standard_other_resolutions/_scan_feedback_%d/_apply_isoxfm_3.0/fb_%d_calc_tshift_resample_volreg_mask_antswarp' % (feedbackRun+1,feedbackRun+1)]]])


    infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name="infosource")
    infosource.iterables = ('subject_id', subject_list)

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['func','mask']),
                         name = 'datasource')
    datasource.inputs.base_directory = data_dir
    datasource.inputs.template = '%s_data_/%s.nii.gz'
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True
    workflow.connect(infosource, 'subject_id', datasource, 'subject_id')

    info_2 = dict(csv=['subject_id', ['feedback_run_%d.csv' % (feedbackRun+1)]])

    infosource2 = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name="infosource_csv")
    infosource2.iterables = ('subject_id', subject_list)

    datasource2 = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['csv']),
                         name = 'datasource_csv')
    datasource2.inputs.base_directory = csv_dir
    datasource2.inputs.template = '%s/%s'
    datasource2.inputs.template_args = info
    datasource2.inputs.sort_filelist = True
    workflow.connect(infosource2, 'subject_id', datasource2, 'subject_id')




    glmNode = pe.Node(interface=fsl.GLM(),name='GLM')
    glmNode.inputs.out_res='residuals.nii.gz'

    workflow.connect([(datasource,glmNode,[('func','in_file')]),
    (datasource,glmNode,[('mask','mask')]),
    (datasource2,glmNode,[('csv','design')])])

    workflow.run()
