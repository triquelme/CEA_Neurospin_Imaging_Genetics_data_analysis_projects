import os
import glob


if __name__ == '__main__':
    #set inputdir, rootdir, outdir and list file names in inputdir
    INDIR = ("/neurospin/ukb/PLINK_OUTPUT/25k_diffusion_TBSS/imp_variants/ch1/"
    "PLINK_parsed/")
    ROOTDIR = "/neurospin/ukb/PLINK_OUTPUT/25k_diffusion_TBSS/imp_variants/"
    fnames = glob.glob(INDIR+"*.pval")
    outdir = os.path.join(ROOTDIR, "concatenated")

    #create output dir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    #for file in filenames, get basename and name output file with basename
    for fname in fnames:
        bname = os.path.basename(fname)
        output = os.path.join(outdir, bname)
        # copy file content from ch1 to concatenated output file
        os.system("cp "+fname+" "+output)
        # loop over files in other ch directories with same names:
        for ch in range(2,23):
            fname = os.path.join(ROOTDIR, "ch"+str(ch), "PLINK_parsed",
                                 bname)
            # tail -n +K output lines starting with the Kth
            # "tail -n+2" output lines starting at line 2 
            # in order to skip file header
            # copy this content to concatenated outputfile
            os.system("tail -n+2 "+fname+" >> "+output)
