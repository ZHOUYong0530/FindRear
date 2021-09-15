import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def binary_search(chromosome, query, data, shift_l=0, shift_r=1):
    index = np.where(data["chr"] == chromosome)
    data_chr = data.iloc[index]
    db = data_chr["pos"]

    db_len = len(db)
    db_start = 0
    db_end = db_len - 1
    db_middle = np.int((db_start + db_end) / 2)

    iterations = np.int(np.log2(db_len)) + 1

    for i in range(iterations):
        if query <= db.iloc[db_middle]:
            db_end = db_middle
            db_middle = np.int((db_start + db_end) / 2)

        if query > db.iloc[db_middle]:
            db_start = db_middle
            db_middle = np.int((db_start + db_end) / 2)

        if db_middle == db_start:
            break
    return data_chr.iloc[db_middle - shift_l:db_middle + shift_r]


def plot_templated_insertion(bpss, sn, outdir):
    svaba_fn = "/lustre/home/yzhou/ESCC/svaba_somatic_intact/" + sn + ".svaba.somatic.sv.vcf"
    delly_fn = "/lustre/home/yzhou/ESCC/delly_somatic_filter/" + sn + ".somatic.filter.vcf"
    data_fn = "/lustre/home/yzhou/ESCC/CNA_200bp_data/" + sn + ".normalized.txt.gz"
    smooth_fn = "/lustre/home/yzhou/ESCC/CNV_logratio/" + sn + ".txt.gz"

    if os.path.exists(data_fn) and os.path.exists(smooth_fn) and os.path.exists(svaba_fn) and os.path.exists(delly_fn):
        data = pd.read_csv(data_fn, sep=" ", compression='gzip')
        smooth_data = pd.read_csv(smooth_fn, sep=" ", compression="gzip")
    else:
        print(sn, "one of files (data,smooth_data,svaba_fn,delly_fn) not exist")
    for bps in bpss:
        bp_len = len(bps)

        fig = plt.figure(figsize=(20, 4 * bp_len))
        plt.suptitle(sn + "\n")

        bp1 = bps[0]
        bp2 = bps[-1]

        for i in np.arange(bp_len):
            j = i + 1

            chr, start, flag = bps[i].split(":")

            if chr == "23":
                chr = "chrX"
            elif chr == "24":
                chr = "chrY"
            else:
                chr = "chr" + np.str(chr)
            start = np.int(start)

            ##        print(chr)
            smoo_tmp = binary_search(chr, start, smooth_data, 100, 100)

            raw_tmp = binary_search(chr, start, data, 100, 100)

            plt.subplot(bp_len, 2, 2 * i + 1)
            plt.plot(smoo_tmp["pos"], smoo_tmp["ratio"])
            plt.vlines(start, ymin=0, ymax=2, linestyle="--")
            plt.grid()
            plt.title(chr + ":" + np.str(start) + ":" + flag)

            plt.subplot(bp_len, 2, 2 * i + 2)
            plt.plot(raw_tmp["pos"], raw_tmp["norm"] / raw_tmp["normal"])
            plt.vlines(start, ymin=0, ymax=2, linestyle="--")
            plt.title(chr + ":" + np.str(start) + ":" + flag)
            plt.grid()
            fig_name = outdir + "/" + sn + bp1 + "-" + bp2 + ".png"
            plt.savefig(fig_name)

        plt.close()


def plot_templated_insertion_graph(bpss, sn, outdir, V_dict):
    for bps in bpss:
        bp_len = len(bps)

        fig = plt.figure(figsize=(20, 4 * bp_len))
        plt.suptitle(sn + "\n")

        bp1 = bps[0]
        bp2 = bps[-1]

        for i in np.arange(bp_len):
            j = i + 1

            chr, start, flag = bps[i].split(":")

            if chr == "23":
                chr = "chrX"
            elif chr == "24":
                chr = "chrY"
            else:
                chr = "chr" + np.str(chr)
            start = np.int(start)

            plt.subplot(bp_len, 2, 2 * i + 1)
            plt.plot(np.arange(201), V_dict[bps[i]].kb_cov)
            plt.vlines(101, ymin=0, ymax=2, linestyle="--")
            plt.grid()
            plt.title(chr + ":" + np.str(start) + ":" + flag)

            plt.subplot(bp_len, 2, 2 * i + 2)
            plt.plot(np.arange(201), V_dict[bps[i]].bp_cov)
            plt.vlines(101, ymin=0, ymax=2, linestyle="--")
            plt.title(chr + ":" + np.str(start) + ":" + flag)
            plt.grid()
            fig_name = outdir + "/" + sn + bp1 + "-" + bp2 + ".png"
            plt.savefig(fig_name)

        plt.close()


def plot_templated_insertion(bpss, sn, outdir):
    svaba_fn = "/lustre/home/yzhou/ESCC/svaba_somatic_intact/" + sn + ".svaba.somatic.sv.vcf"
    delly_fn = "/lustre/home/yzhou/ESCC/delly_somatic_filter/" + sn + ".somatic.filter.vcf"
    data_fn = "/lustre/home/yzhou/ESCC/CNA_200bp_data/" + sn + ".normalized.txt.gz"
    smooth_fn = "/lustre/home/yzhou/ESCC/CNV_logratio/" + sn + ".txt.gz"

    if os.path.exists(data_fn) and os.path.exists(smooth_fn) and os.path.exists(svaba_fn) and os.path.exists(delly_fn):
        data = pd.read_csv(data_fn, sep=" ", compression='gzip')
        smooth_data = pd.read_csv(smooth_fn, sep=" ", compression="gzip")
    else:
        print(sn, "one of files (data,smooth_data,svaba_fn,delly_fn) not exist")
    for bps in bpss:
        bp_len = len(bps)

        fig = plt.figure(figsize=(20, 4 * bp_len))
        plt.suptitle(sn + "\n")

        bp1 = bps[0]
        bp2 = bps[-1]

        for i in np.arange(bp_len):
            j = i + 1

            chr, start, flag = bps[i].split(":")

            if chr == "23":
                chr = "chrX"
            elif chr == "24":
                chr = "chrY"
            else:
                chr = "chr" + np.str(chr)
            start = np.int(start)

            ##        print(chr)
            smoo_tmp = binary_search(chr, start, smooth_data, 100, 100)

            raw_tmp = binary_search(chr, start, data, 100, 100)

            plt.subplot(bp_len, 2, 2 * i + 1)
            plt.plot(smoo_tmp["pos"], smoo_tmp["ratio"])
            plt.vlines(start, ymin=0, ymax=2, linestyle="--")
            plt.grid()
            plt.title(chr + ":" + np.str(start) + ":" + flag)

            plt.subplot(bp_len, 2, 2 * i + 2)
            plt.plot(raw_tmp["pos"], raw_tmp["norm"] / raw_tmp["normal"])
            plt.vlines(start, ymin=0, ymax=2, linestyle="--")
            plt.title(chr + ":" + np.str(start) + ":" + flag)
            plt.grid()
            fig_name = outdir + "/" + sn + bp1 + "-" + bp2 + ".png"
            plt.savefig(fig_name)

        plt.close()
