
import os
import glob
import pandas as pd
import numpy as np


if __name__ == '__main__':
    DIR = '/neurospin/ukb/PLINK_OUTPUT/25k_diffusion_TBSS/imp_variants/concatenated'
    fnames = glob.glob(DIR+'/*.sel')
    OUTDIR = '/neurospin/tmp/yleguen/formatted'
    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR)
    OUTDIR2 = '/neurospin/ukb/tmp/THOMAS'
    fnames2 = []
    for fname in fnames:
        df = pd.read_csv(fname, sep='\t')
        df2 = df[["CHR", "SNP", "BP", "P"]]
        df2 = df2.rename(columns={"BP": "POS"})
        basename = os.path.basename(fname)
        basename, extension = os.path.splitext(basename)
        basename = basename.replace("covar_SexAgeMRI_10PCSArray_25k.all_", "")
 	phenotype = basename.split(".",1)[1]
 	df2["PHENOTYPE"] = phenotype
        # add column group
 	# strip basename to get only measure type (FA, ICVF, MO...)
        measure_type = basename.split("_",1)[0]
 	df2["GROUP"] = measure_type
        output = os.path.join(OUTDIR, os.path.basename(fname))
 	df2.to_csv(output, sep='\t', header=True, index=False)
        fnames2.append(output)


    output = os.path.join(OUTDIR, "concatenated_all_sel.sel")
    df = pd.concat(pd.read_csv(fname, sep='\t') for fname in fnames2)
    df.to_csv(output, sep='\t', header=True, index=False)

    for ch in df['CHR'].unique():
        df0 = df.loc[df.CHR == ch]
        df1  = pd.DataFrame()
        df1['SNP'] = list(set(df0['SNP']))
        output = os.path.join(OUTDIR, 'chr'+str(ch)+'_list_snps.txt')
        df1.to_csv(output, sep='\t', header=False, index=False)


        
    DIR_PLINK = '/neurospin/ukb/GENETIC_DATA_500k/IMPUTED_DATA/plink_extracted_25k'
    for ch in df['CHR'].unique():

        list_snp = os.path.join(OUTDIR, 'chr'+str(ch)+'_list_snps.txt')
        geno = os.path.join(DIR_PLINK, 'chr'+str(ch)+'_imp')
        geno_filt = os.path.join(OUTDIR2, 'chr'+str(ch)+'_imp_filtered')

        cmd = " ".join(["plink",
                        "--bfile %s"%geno,
                        "--extract %s"% list_snp,
                        "--make-bed",
                        "--out %s"% geno_filt])
        #print cmd
        #os.system(cmd)
        output = os.path.join(OUTDIR, "ld_matrix_ch"+str(ch))
        cmd = " ".join(["plink",
                        "--bfile %s"% geno_filt,
                        "--matrix",
                        "--r2",
                        "--out %s"% output])
        #print cmd
        #os.system(cmd)
        
        #m_ld = np.loadtxt(output+'.ld')
        
    BIM_DIR = "/neurospin/ukb/tmp/THOMAS/"
    LD_DIR = "/neurospin/tmp/yleguen/formatted/"

    # reduce to top snps and removes LD snps
    input_bim = (BIM_DIR + "chr10_imp_filtered.bim")
    bim = pd.read_csv(input_bim, sep='\t', header=None)

    input_ld = (LD_DIR + "ld_matrix_ch10.ld")
    ld = pd.read_csv(input_ld, sep=' ', header=None)
    ld = ld.drop(len(ld.columns) - 1, 1)
    ld.index = bim[1]
    ld.columns = bim[1]

    # Too many variables undeclared in the lines below
    # from Thomas script please review ..
    """
    for lab, row in ld.iterrows():
        if (row[lab] >= 5):
            ld_snp.append(ld[lab])
        ld_snp_list.append(ld_snp)
        ld_snp = []

    for lab, row in ld.iterrows():
        for i in row:
            if i > 0.5:
                ld_dict[row.index] = lab
    """
    df.index = df.SNP
    LD_THRESHOLD = 0.001
    # DataFrame that will contain all the data kept
    df_result = pd.DataFrame()
    # So the goal here is to only keep top snps not in ld with others for each group and each tract
    # NB for now we assume the top SNP in a LD region will be the same for various phenotypes
    # if that's not the case we will adjust this maybe to merge the top SNPs as a single one
    
    for ch in df['CHR'].unique():
        # Work on one chromosome at the time
        df0 = df.loc[df.CHR == ch]
        input_bim = (BIM_DIR + "chr"+str(ch)+"_imp_filtered.bim")
        bim = pd.read_csv(input_bim, sep='\t', header=None)

        # Load LD matrix
        input_ld = (LD_DIR + "ld_matrix_ch"+str(ch)+".ld")
        ld = pd.read_csv(input_ld, sep=' ', header=None)
        ld = ld.drop(len(ld.columns) - 1, 1)
        # Identified rows and columns in this LD matrix
        ld.index = bim[1]
        ld.columns = bim[1]

        # Then let's work feature by feature
        for gr in df0['GROUP'].unique():
            df1 = df0.loc[df0.GROUP == gr]
            # Then work tract by tract
            for ph in df1['PHENOTYPE'].unique():
                df2 = df1.loc[df1.PHENOTYPE == ph]
                # Sort pvalues
                df2 = df2.sort_values('P', ascending=True)
                df2.index = df2.SNP
                count = 0
                # Then loop through snps and keep the one not in ld
                while len(list(df2.SNP)) > count:
                    # Select current top snp
                    snp = list(df2.SNP)[count]
                    for snp2 in list(df2.SNP)[count+1:]:
                        # Weirdly there are duplicates SNP in bim
                        # This is probably due to the conversion
                        # from .bgen to .bed/.bim/.fam format
                        # for now we ignore this problem that should be dealt
                        # with by looking at the reference that should be
                        # different for each line
                        """
                        # if unique snp this lines are sufficient
                        if ld.loc[snp][snp2] > LD_THRESHOLD:
                            # Drop SNP in LD
                            df2 = df2.drop(snp2)
                        """
                        # Normal case SNP is not duplicated
                        if ld.loc[snp][snp2].size == 1:
                            if ld.loc[snp][snp2] > LD_THRESHOLD:
                                # Drop SNP in LD
                                df2 = df2.drop(snp2)
                        # Unusual case SNP is duplicated
                        elif max(ld.loc[snp][snp2]) > LD_THRESHOLD:
                            print ld.loc[snp][snp2].size
                            print ld.loc[snp][snp2]
                            # Check that we did not already removed it
                            if snp2 in df2.index:
                                # Drop SNP in LD
                                df2 = df2.drop(snp2)
                        

                        
                    count += 1
                
                df_result = df_result.append(df2)
    
    output = os.path.join(OUTDIR, 'final_result_LD0.001.sel')
    df_result = df_result.sort_values(['CHR', 'POS'])
    df_result.to_csv(output, sep='\t', index=False, header=True)
    df_result.index = range(df_result.shape[0])
    df_result = df_result.loc[df_result.P < 5e-8/(48*6)]
    df_result = df_result.sort_values(['CHR', 'POS'])
    output = os.path.join(OUTDIR, 'final_result_bonf_corrected_LD0.001.sel')
    df_result.to_csv(output, sep='\t', index=False, header=True)

    # DataFrame that will contain all the data kept
    df_result2 = pd.DataFrame()
    
    # Here we gonna replace the snps that are in LD by a unique top snp
    # see the NB remark above
    for ch in df_result['CHR'].unique():
        # Work on one chromosome at the time
        df0 = df_result.loc[df_result.CHR == ch]
        input_bim = (BIM_DIR + "chr"+str(ch)+"_imp_filtered.bim")
        bim = pd.read_csv(input_bim, sep='\t', header=None)

        # Load LD matrix
        input_ld = (LD_DIR + "ld_matrix_ch"+str(ch)+".ld")
        ld = pd.read_csv(input_ld, sep=' ', header=None)
        ld = ld.drop(len(ld.columns) - 1, 1)
        # Identified rows and columns in this LD matrix
        ld.index = bim[1]
        ld.columns = bim[1]
        
        # Sort pvalues
        df0 = df0.sort_values('P', ascending=True)
        df0.index = range(df0.shape[0])
        list_done = []
        list_SNP = ['0']*df0.shape[0]
        
        for k in range(df0.shape[0]):

            if not k in list_done:
                snp_k = df0.loc[k].SNP
                list_SNP[k] = snp_k
                for j in range(k+1, df0.shape[0]):
                    snp_j = df0.loc[j].SNP
                    # Normal case SNP is not duplicated
                    if ld.loc[snp_k][snp_j].size == 1:
                        if ld.loc[snp_k][snp_j] > LD_THRESHOLD:
                            #df0.loc[j].SNP = snp_k
                            list_SNP[j] = snp_k
                            list_done.append(j)
                    # Unusual case SNP is duplicated
                    elif max(ld.loc[snp_k][snp_j]) > LD_THRESHOLD:
                        #df0.loc[j].SNP = snp_k
                        list_SNP[j] = snp_k
                        list_done.append(j)
        df0['SNP'] = list_SNP
        df_result2 = df_result2.append(df0)

    df_result2 = df_result2.sort_values(['CHR', 'POS'])
    output = os.path.join(OUTDIR, 'final_result_bonf_corrected_LD0.001_unique.sel')
    df_result2.to_csv(output, sep='\t', index=False, header=True)
    df_result3 = df_result2.drop_duplicates('SNP')
    df_result3 = df_result3.sort_values(['CHR', 'POS'])
    output = os.path.join(OUTDIR, 'final_result_bonf_corrected_LD0.001_UNIQUE.sel')
    df_result3.to_csv(output, sep='\t', index=False, header=True)

    
    # DataFrame that will contain all the data kept
    df_result3 = pd.DataFrame()
    MAX_DISTANCE = 200000
    for ch in df_result2['CHR'].unique():
        # Work on one chromosome at the time
        df0 = df_result2.loc[df_result2.CHR == ch]
        # Sort pvalues
        df0 = df0.sort_values('P', ascending=True)
        df0.index = range(df0.shape[0])
        list_done = []
        list_SNP = ['0']*df0.shape[0]
        for k in range(df0.shape[0]):
            if not (k in list_done):
                snp_k = df0.loc[k].SNP
                list_SNP[k] = snp_k
                for j in range(k+1, df0.shape[0]):
                    if  abs(df0.loc[k].POS-df0.loc[j].POS) < MAX_DISTANCE:
                        #df0.loc[j].SNP = snp_k
                        list_SNP[j] = snp_k
                        list_done.append(j)
        df0['SNP'] = list_SNP           
        df_result3 = df_result3.append(df0)
    df_result3 = df_result3.sort_values(['CHR', 'POS'])
    output = os.path.join(OUTDIR, 'final_result_bonf_corrected_LD0.001_unique_200kBP.sel')
    df_result3.to_csv(output, sep='\t', index=False, header=True)
    df_result3.drop_duplicates('SNP')
    df_result3 = df_result3.drop_duplicates('SNP')
    df_result3 = df_result3.sort_values(['CHR', 'POS'])
    output = os.path.join(OUTDIR, 'final_result_bonf_corrected_LD0.001_200kBP_UNIQUE.sel')
    df_result3.to_csv(output, sep='\t', index=False, header=True)
