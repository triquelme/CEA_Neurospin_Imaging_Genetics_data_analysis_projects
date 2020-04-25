import os

INDIR = "/neurospin/ukb/PLINK_OUTPUT/25k_diffusion_TBSS/imp_variants/"

fnames = []
nb_variants_list = []

for ch in range(1, 23):
    fname = os.path.join(INDIR, "ch" + str(ch), "PLINK_brut",
                         "covar_SexAgeMRI_10PCSArray_25k.all_FA_skeletonised.log")
    fnames.append(fname)

for fname in fnames:
    with open(fname) as f:
        all_lines = f.read().splitlines()
        line = all_lines[37]
        nb_variants = int(line.split(' ')[0])
        nb_variants_list.append(nb_variants)

total_nb_variant = sum(nb_variants_list)

print(total_nb_variant)