import CCD_packages
import numpy as np
from nilearn import image
from pandas import read_csv

GroupDF,motionInfo=CCD_packages.getCCDSubjectData(saveMotionInfo=False)
goodsubj,badsubj = CCD_packages.getSubjectList(GroupDF=GroupDF,motionThresh=1,poor_performer=15)

zmap_filenames=[]
fc='PC'
secondlevel_folder_names=['noFeedback','Feedback','train']
fb=2
behavioral_target=np.arctanh(np.array(GroupDF[np.all([GroupDF.Subject_ID.isin(goodsubj),GroupDF.FB=='FEEDBACK'],axis=0)].groupby('Subject_ID')['modelcorr'].mean()))


for subj in goodsubj:
    zmap_filenames.append('/home/jmuraskin/Projects/CCD/working_v1/seed-to-voxel/%s/%s/%s_%s.nii.gz' % (fc,secondlevel_folder_names[fb],fc,subj))

mask_filename='/home/jmuraskin/Projects/CCD/working_v1/seg_probabilities/grey_matter_mask-20-percent.nii.gz'

from scipy.stats import zscore
#load phenotypic data
phenoFile='/home/jmuraskin/Projects/CCD/Pheno/narsad+vt_new.csv'
pheno=read_csv(phenoFile)
pheno=pheno.set_index('participant')

ages=zscore(pheno.loc[goodsubj]['V1_DEM_001'])

mf=zscore(pheno.loc[goodsubj]['V1_DEM_002'])

motionTest=read_csv('/home/jmuraskin/Projects/CCD/CCD-scripts/analysis/CCD_meanFD.csv')
meanFD=zscore(motionTest[motionTest.FB=='FEEDBACK'][motionTest.Subject_ID.isin(goodsubj)]['train_meanFD'])


imgs=image.concat_imgs(zmap_filenames)

clean_imgs=image.clean_img(imgs,confounds=[ages,mf,meanFD],detrend=False,standardize=True)


from nilearn.decoding import SpaceNetRegressor

decoder = SpaceNetRegressor(mask=mask_filename, penalty="tv-l1",
                            eps=1e-1,  # prefer large alphas
                            memory="nilearn_cache")

decoder.fit(clean_imgs, behavioral_target)
