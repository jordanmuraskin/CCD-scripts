"""
Producing single subject maps of seed-to-voxel correlation
==========================================================

This example shows how to produce seed-to-voxel correlation maps for a single
subject based on resting-state fMRI scans. These maps depict the temporal
correlation of a **seed region** with the **rest of the brain**.

This example is an advanced one that requires manipulating the data with numpy.
Note the difference between images, that lie in brain space, and the
numpy array, corresponding to the data inside the mask.
"""

# author: Franz Liem and Jordan Muraskin


import argparse
import numpy as np
from nilearn import input_data
import os
from pandas import DataFrame,read_csv


parser = argparse.ArgumentParser(description='Run First Level Functional Connectivity for Neurofeedback Data')
parser.add_argument('-globalSR', help='Option to run with global signal regression',required=False,default=0,type=int)
parser.add_argument('-name', help='ROI name for foldernaming',required=True,default='ROI',type=str)
parser.add_argument('-sphere', help='Sphere size',required=False,default=6,type=int)

args = parser.parse_args()

globalSR=args.globalSR


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


def subjectinfo(subject_id,getFeedback=True):


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
        return feedback+1
    if not getFeedback:
        return noFeedback+1

def signFlipInfo(subject_id,r):
    SubjInfo = read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/NARSAD_stimulus_JM.csv')
    SubjInfo.set_index('JM_INTERNAL',inplace=True)

    if r==1:
        paradigmType=SubjInfo.loc[subject_id]['SCAN_1_PARADIGM']
    elif r==2:
        paradigmType=SubjInfo.loc[subject_id]['SCAN_2_PARADIGM']
    if paradigmType==0 or paradigmType == 2:
        signFlip=1
    elif paradigmType==1 or paradigmType == 3:
        signFlip=-1

    return signFlip

topDir='/home/jmuraskin/Projects/CCD/working_v1/ROI-timeseries/'
if not os.path.exists(topDir):
    os.mkdir(topDir)

topDir=topDir + '/' + args.name
if not os.path.exists(topDir):
    os.mkdir(topDir)


for indx,fb in enumerate(['noFeedback','Feedback']):

    baseDir=topDir + '/' + fb
    if not os.path.exists(baseDir):
        os.mkdir(baseDir)

    for subject_id in subject_list:

        ##########################################################################
        # Getting the data
        # ----------------

        # We will work with the first subject of the adhd data set.
        # adhd_dataset.func is a list of filenames. We select the 1st (0-based)
        # subject by indexing with [0]).
        # from nilearn import datasets
        scan=subjectinfo(subject_id,getFeedback=indx)
        signFlip=signFlipInfo(subject_id,scan)

        func_filename = '/home/jmuraskin/Projects/CCD/CPAC-out/pipeline_CCD_v1/%s_data_/functional_mni_other_resolutions_smooth/_scan_feedback_%d/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_apply_isoxfm_3.0/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global%d.motion1.quadratic1.gm0.compcor1.csf1/_fwhm_6/residual_antswarp_maths.nii.gz' % (subject_id,scan,globalSR)


        ##########################################################################
        # Time series extraction
        # ----------------------
        #
        # We are going to extract signals from the functional time series in two
        # steps. First we will extract the mean signal within the **seed region of
        # interest**. Second, we will extract the **brain-wide voxel-wise time series**.
        #
        # We will be working with one seed sphere in the Posterior Cingulate Cortex,
        # considered part of the Default Mode Network.

        labels=['Posterior Cingulate','Angular Gyrus','Paracingulate Gyrus','Middle Cingulate','Superior Parietal Lobule','Thalamus','Anterior Insula','Anteior Cingulate']

        coords = [(0,-52,28),(54,-51,15),(0,50,-6),(3,-15,42),(36,-48,45),(15,-9,9),(42,21,-3),(0,36,24)]

        ##########################################################################
        # We use :class:`nilearn.input_data.NiftiSpheresMasker` to extract the
        # **time series from the functional imaging within the sphere**. The
        # sphere is centered at pcc_coords and will have the radius we pass the
        # NiftiSpheresMasker function (here 8 mm).
        #
        # The extraction will also detrend, standardize, and bandpass filter the data.
        # This will create a NiftiSpheresMasker object.


        seed_masker = input_data.NiftiSpheresMasker(
            coords, radius=args.sphere, standardize=True,t_r=2., verbose=1)

        ##########################################################################
        # Then we extract the mean time series within the seed region while
        # regressing out the confounds that
        # can be found in the dataset's csv file
        seed_time_series = seed_masker.fit_transform(func_filename)

        ##########################################################################
        # Next, we can proceed similarly for the **brain-wide voxel-wise time
        # series**, using :class:`nilearn.input_data.NiftiMasker` with the same input
        # arguments as in the seed_masker in addition to smoothing with a 6 mm kernel
        brain_masker = input_data.NiftiMasker(standardize=True, t_r=2.,verbose=1)

        ##########################################################################
        # Then we extract the brain-wide voxel-wise time series while regressing
        # out the confounds as before
        brain_time_series = brain_masker.fit_transform(func_filename)

        df=DataFrame(data=brain_time_series,columns=labels)
        df.to_csv('%s/%s.csv' % (baseDir,subject_id))
