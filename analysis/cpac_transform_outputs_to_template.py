#!/usr/bin/env python

def gather_filepaths(output_folder_path, outputs_list=None):

    import os
    filepaths = []

    for root, dirs, files in os.walk(output_folder_path):
        # loops through every file in the directory
        for filename in files:
            # checks if the file is a nifti (.nii.gz)
            if outputs_list:
                for output in outputs_list:
                    if output in os.path.join(root, filename):
                        filepaths.append(os.path.join(root, filename))
            else:
                filepaths.append(os.path.join(root, filename))

    if len(filepaths) == 0:
        err = "\n\n[!] No filepaths were found given the output folder!\n\n"
        raise Exception(err)

    return filepaths


def filepath_list_to_dict(filepaths_list, output_folder):

    filepath_dict = {}

    for filepath in filepaths_list:
        half = filepath.split(output_folder)[1]
        dir_levels = [x for x in half.split("/") if x != ""]
        sub_id = dir_levels[0]
        output_name = dir_levels[1]

        try:
            filepath_dict[sub_id].update({output_name: filepath})
        except KeyError:
            filepath_dict[sub_id] = {output_name: filepath}

    return filepath_dict


def antswarp_output_to_template(args_tuple): #filepath, warps_list, ref_image, outfile, interp="Linear"):

    import subprocess

    filepath = args_tuple[0]
    warps_list = args_tuple[1]
    ref_image = args_tuple[2]
    outfile = args_tuple[3]
    interp = args_tuple[4]

    print "Warping %s.." % filepath

    if not interp:
        interp = "Linear"

    ants_cmd = ["antsApplyTransforms", "--default-value 0", \
                "--dimensionality 3", "--input", filepath, \
                "--input-image-type 0", "--interpolation", interp, \
                "--output", outfile, "--reference-image", ref_image]

    for warp in warps_list:
        ants_cmd.append("--transform")
        ants_cmd.append(warp)

    try:
        retcode = subprocess.check_output(ants_cmd)
    except Exception as e:
        err = "\n\nLocals %s\n\n[!] antsApplyTransforms did not complete " \
              "successfully.\n\nCommand:\n%s\n\nError details: %s\n\n" \
              % (locals(), ants_cmd, e)
        raise Exception(err)


def main():

    import os
    import argparse
    from multiprocessing import Pool

    parser = argparse.ArgumentParser()
 
    parser.add_argument("cpac_output_directory", type=str, \
                            help="the path to the CPAC run's output " \
                                 "directory")
 
    parser.add_argument("reference_image", type=str, \
                            help="the filepath to the template image you " \
                                 "want to use for applying the transforms")
 
    parser.add_argument('--output_folder', type=str, \
                            help="the directory to write the warped files to")

    parser.add_argument("--interpolation", type=str, \
                            help="which interpolation to use when applying " \
                                 "the ANTS warps - available: 'Linear', " \
                                 "'NearestNeighbor', 'MultiLabel" \
                                 "[<sigma=imageSpacing>,<alpha=4.0>]', " \
                                 "'Gaussian[<sigma=imageSpacing>," \
                                 "<alpha=1.0>]', 'BSpline[<order=3>]', " \
                                 "'CosineWindowedSinc','WelchWindowedSinc', "\
                                 "'HammingWindowedSinc', " \
                                 "'LanczosWindowedSinc' - default: Linear")

    parser.add_argument("--num_cores", type=int, \
                            help="number of cores to use (number of warps " \
                                 "to calculate in parallel - default: 1")
 
    args = parser.parse_args()

    # defaults
    if not args.output_folder:
        output_folder = os.getcwd()
    else:
        output_folder = args.output_folder

    if not args.num_cores:
        num_cores = 1
    else:
        num_cores = args.num_cores

    # run it!
    outputs_list = ["anatomical_gm_mask", "anatomical_csf_mask", \
                    "anatomical_wm_mask", "anatomical_to_mni_nonlinear_xfm", \
                    "ants_affine_xfm", "ants_initial_xfm", "ants_rigid_xfm"]

    to_warp_list = ["anatomical_gm_mask", "anatomical_csf_mask", \
                    "anatomical_wm_mask"]

    filepaths = gather_filepaths(args.cpac_output_directory, outputs_list)

    file_dict = filepath_list_to_dict(filepaths, args.cpac_output_directory)

    pool = Pool(processes=num_cores)

    # per file per subject!
    args_list = []
    for sub_id in file_dict.keys():

        # reverse order for ANTS
        warps_list = [file_dict[sub_id]["anatomical_to_mni_nonlinear_xfm"],
                      file_dict[sub_id]["ants_affine_xfm"],
                      file_dict[sub_id]["ants_rigid_xfm"],
                      file_dict[sub_id]["ants_initial_xfm"]]

        for resource in to_warp_list:

            filepath = file_dict[sub_id][resource]

            outname = "_".join([resource, "to_template"])
            outdir = os.path.join(output_folder, "warp_outputs", sub_id, 
                outname)
            outfile = os.path.join(outdir, "".join([outname, ".nii.gz"]))

            if not os.path.exists(outdir):
                os.makedirs(outdir)

            args_list.append((filepath, warps_list, args.reference_image, outfile, args.interpolation))
    
    #for args_tuple in args_list:
    #    antswarp_output_to_template(args_list)
    pool.map(antswarp_output_to_template, args_list)


if __name__ == "__main__":
    main()
