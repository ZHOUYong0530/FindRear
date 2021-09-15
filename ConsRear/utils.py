import numpy as np
import pandas as pd


def get_numeric_chr(chr_name):
    chr_num = chr_name.split("chr")[1]
    if chr_num == "X":
        return 23
    elif chr_num == "Y":
        return 24
    else:
        return np.int(chr_num)


def get_alpha(chromosome):
    if chromosome == 23:
        chromosome = "chrX"
    elif chromosome == 24:
        chromosome = "chrY"
    else:
        chromosome = "chr" + np.str(chromosome)
    return chromosome


def nor_coverage_bp(data_fn):
    pd.options.mode.use_inf_as_na = True
    nc200 = pd.read_csv(data_fn, sep=" ", usecols=["chr", "pos", "norm", "normal"])
    nc200["nc"] = nc200["norm"].div(nc200["normal"])
    nc200.fillna(method="pad", limit=2, inplace=True)
    nc200.dropna(inplace=True)
    nc200.index = np.arange(nc200.shape[0])
    nc200 = nc200.round({"nc": 3})
    return nc200


def nor_coverage_kb(smooth_fn):
    pd.options.mode.use_inf_as_na = True
    nc = pd.read_csv(smooth_fn, sep=" ")
    nc.fillna(method="pad", limit=2, inplace=True)
    nc.dropna(inplace=True)
    nc.index = np.arange(nc.shape[0])
    nc = nc.round({"ratio": 3})
    return nc


def write_coverage_file(V, nc200, coverage_outfile, status):
    for v in V:
        chromosome, start, flag = v.get_info().split(":")
        chromosome = np.int(chromosome)
        start = np.int(start)
        chromosome = get_alpha(chromosome)
        nc_chrom = nc200[nc200["chr"] == chromosome]
        n = nc_chrom.shape[0]
##print(v.get_info())
        ix = nc_chrom["pos"].searchsorted(start)
        if ix >= n - 100:
            nc_chrom_target = (np.zeros(101 + (ix - n)) + nc_chrom.iloc[n - 1, 4]).tolist()
            nc_chrom_target = nc_chrom.iloc[ix - 100:n, 4].tolist() + nc_chrom_target
        elif ix >= 100:
            nc_chrom_target = nc_chrom.iloc[ix - 100:ix + 101, 4].tolist()
        elif ix < 100:
            nc_chrom_target = (np.zeros(100 - ix) + nc_chrom.iloc[ix, 4]).tolist()
            nc_chrom_target += nc_chrom.iloc[0:ix + 101, 4].tolist()

        if len(nc_chrom_target) != 201:
            print(v.get_info(), ix, n, len(nc_chrom_target))
## print(v.get_info(),nc_chrom.iloc[ix,1],nc_chrom.iloc[ix,4])
##print(v.get_info(),nc_chrom.iloc[ix,1],nc_chrom.iloc[ix,4],nc_chrom_target,sep="\t")
        if status == "bp":
            v.bp_cov = nc_chrom_target
        elif status == "kb":
            v.kb_cov = nc_chrom_target
    coverage_to_disk = {}
    if status == "bp":
        for v in V:
            coverage_to_disk[v.get_info()] = v.bp_cov
    elif status == "kb":
        for v in V:
            coverage_to_disk[v.get_info()] = v.kb_cov

    pd.DataFrame(coverage_to_disk).to_csv(coverage_outfile, index=True, sep="\t")

def cov_median(bp_cov):
    ix = np.int(len(bp_cov) / 2)
    left = np.median(bp_cov[0:ix])
    right = np.median(bp_cov[ix + 1:])
    return np.around(left, 3), np.around(right, 3)

if __name__ == "__main__":
    chrom = 14
    print(get_alpha(chrom))
