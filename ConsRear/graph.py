import numpy as np
import pandas as pd


class BreakPoint:
    color = ""
    pi = ""

    def __init__(self, bp, bp_cov=None, kb_cov=None):
        self.bp = bp
        self.chrom, self.start, self.flag = bp.split(":")
        if bp_cov != "None":
            self.bp_cov = bp_cov
        else:
            self.bp_cov = None

        if kb_cov != "None":
            self.kb_cov = kb_cov
        else:
            self.kb_cov = None

    def get_cov(self):  # cov: coverage of breakpoint: left_coverage, right coverge
        return self.bp_cov

    def get_cpf(self):  # cps : chromsome,position,flag (orientation)
        return self.chrom, self.start, self.flag

    def get_info(self):
        return ":".join([np.str(self.chrom), np.str(self.start), self.flag])


class BreakEdge:
    def __init__(self, bp1, bp2, status):
        self.bp1 = bp1
        self.bp2 = bp2
        self.status = status

    def get_bps(self):
        return self.bp1, self.bp2

    def get_status(self):
        return self.status


def combo(line):
    bp1 = ":".join([np.str(line.chr_1), np.str(line.pos_1), line.flag_1])
    bp2 = ":".join([np.str(line.chr_2), np.str(line.pos_2), line.flag_2])
    return bp1, bp2


def cohort_sort(sv_cohort):
    return sv_cohort.sort_values(['SN', 'chr_1', 'pos_1'], ascending=[True, True, True])


def cohort_sort_dd(sv_cohort):
    sv_cohort_sort = sv_cohort.sort_values(['SN', 'chr_1', 'pos_1'], ascending=[True, True, True])
    new_order = ["SN", "chr_2", "pos_2", "flag_2", "chr_1", "pos_1", "flag_1"]
    sv_cohort_sort = sv_cohort_sort[["SN", "chr_1", "pos_1", "flag_1", "chr_2", "pos_2", "flag_2"]]
    sv_cohort_sort_rev = sv_cohort_sort[new_order]
    sv_cohort_sort_rev.columns = sv_cohort_sort.columns

    frames = [sv_cohort_sort, sv_cohort_sort_rev]
    sv_cohort_dd = pd.concat(frames)
    sv_cohort_dd_sorted = sv_cohort_dd.sort_values(['SN', 'chr_1', 'pos_1'], ascending=[True, True, True])

    return sv_cohort_dd_sorted


def seek_overlaps(chromosome, start, end, svaba_d, invers=0):
    chr1 = svaba_d.loc[:, "chr_1"]
    pos1 = np.int0(svaba_d.loc[:, "pos_1"])
    flag1 = svaba_d.loc[:, "flag_1"]
    chr2 = svaba_d.loc[:, "chr_2"]
    pos2 = np.int0(svaba_d.loc[:, "pos_2"])
    flag2 = svaba_d.loc[:, "flag_2"]

    if invers == 0:
        ix1 = (chr1 == chr2) & (chr1 == chromosome) & (pos2 >= start) & (pos2 <= end)
        ix2 = (chr1 == chr2) & (chr1 == chromosome) & (pos1 >= start) & (pos1 <= end)
        ix = np.ravel(np.where(ix1 | ix2))
        return svaba_d.loc[ix]
    elif invers == 1:
        ix1 = (chr1 == chr2) & (chr1 == chromosome) & (pos2 >= start) & (pos2 <= end)
        ix2 = (chr1 == chr2) & (chr1 == chromosome) & (pos1 >= start) & (pos1 <= end)
        ix3 = (chr1 == chr2) & (chr1 == chromosome) & (pos1 < start) & (pos2 > end)
        ix = np.ravel(np.where(ix1 | ix2 | ix3))
        return svaba_d.iloc[np.ravel(ix)]
    elif invers == 2:
        ix1 = (chr1 == chromosome) & (pos1 >= start) & (pos1 <= end)
        ix2 = (chr2 == chromosome) & (pos2 >= start) & (pos2 <= end)
        ix = np.ravel(np.where(ix1 | ix2))
        return svaba_d.iloc[np.ravel(ix)]
    elif invers == 3:
        ix1 = (chr1 == chromosome) & (pos1 > start) & (pos1 < end)
        ix2 = (chr2 == chromosome) & (pos2 > start) & (pos2 < end)
        ix = np.ravel(np.where(ix1 | ix2))
        return svaba_d.iloc[np.ravel(ix)]
    elif invers == 4:
        ix3 = (chr1 == chr2) & (chr1 == chromosome) & (pos1 < start) & (pos2 > end)
        ix = np.ravel(np.where(ix3))
        return svaba_d.iloc[np.ravel(ix)]
    else:
        ix1 = (chr1 == chr2) & (chr1 == chromosome) & (pos2 >= start) & (pos2 <= end) & (flag1 == flag2) & (
                (pos2 - pos1) > 20000)
        ix2 = (chr1 == chr2) & (chr1 == chromosome) & (pos1 >= start) & (pos1 <= end) & (flag1 == flag2) & (
                (pos2 - pos1) > 20000)
        ix = np.ravel(np.where(ix1 | ix2))
        return svaba_d.loc[ix]


def cohort_rm_simple(sv_cohort):
    sv_cohort_sorted = cohort_sort(sv_cohort)
    # n = sv_cohort_sorted.shape[0]
    simple_flag = []
    for i, line in sv_cohort_sorted.iterrows():
        if line["chr_1"] != line["chr_2"]:
            simple_flag.append("non-simple")
            continue
        elif line["flag_1"] == line["flag_2"]:
            simple_flag.append("non-simple")
            continue
        else:
            chrom = line["chr_1"]
            start = line["pos_1"]
            end = line["pos_2"]
            interval_sv = seek_overlaps(chrom, start, end, sv_cohort_sorted, invers=3)
            if interval_sv.shape[0] == 0:
                simple_flag.append("simple")
            else:
                simple_flag.append("non-simple")

    sv_cohort_sorted["simple"] = simple_flag
    return sv_cohort_sorted


def find_genomicEdges(sv_cohort_sort, V_dict, V, E, gs_allow=2, distance=5000000):
    new_order = ["SN", "chr_2", "pos_2", "flag_2", "chr_1", "pos_1", "flag_1"]
    sv_cohort_sort = sv_cohort_sort[["SN", "chr_1", "pos_1", "flag_1", "chr_2", "pos_2", "flag_2"]]
    sv_cohort_sort_rev = sv_cohort_sort[new_order]
    sv_cohort_sort_rev.columns = sv_cohort_sort.columns

    frames = [sv_cohort_sort, sv_cohort_sort_rev]
    sv_cohort_dd = pd.concat(frames)
    sv_cohort_dd_sorted = sv_cohort_dd.sort_values(['SN', 'chr_1', 'pos_1'], ascending=[True, True, True])

    # find nearby geomic edges
    n_rows = sv_cohort_dd_sorted.shape[0]
    sv_cohort_dd_sorted.index = np.arange(n_rows)

    for i, line in sv_cohort_dd_sorted.iterrows():
        gs_allowed = gs_allow
        bp1, bp2 = combo(line)

        bp1_v = BreakPoint(bp1) if bp1 not in V_dict else V_dict[bp1]
        bp2_v = BreakPoint(bp2) if bp2 not in V_dict else V_dict[bp2]

        V_dict[bp1] = bp1_v
        V_dict[bp2] = bp2_v

        V.append(bp1_v)
        V.append(bp2_v)

        edge1 = BreakEdge(bp1_v, bp2_v, "sv")  # SV edge
        E[(bp1_v, bp2_v)] = edge1

        # search for genomic edge for bp1
        chrom, start, flag = bp1.split(":")
        start = np.int(start)
        chrom = np.int(chrom)
        if flag == "-":
            while gs_allowed:
                ix = i + gs_allow - gs_allowed + 1
                if (ix >= n_rows) or (sv_cohort_dd_sorted.iloc[ix, 1] != chrom):
                    break
                if sv_cohort_dd_sorted.iloc[ix, 3] == "+":
                    bp_next = ":".join(np.str(j) for j in sv_cohort_dd_sorted.iloc[ix, 1:4])

                    if np.abs(np.int0(sv_cohort_dd_sorted.iloc[ix, 2]) - np.int0(
                            sv_cohort_dd_sorted.iloc[i, 2])) > distance:
                        break
                    else:
                        bp_next_v = BreakPoint(bp_next) if bp_next not in V_dict else V_dict[bp_next]
                        V_dict[bp_next] = bp_next_v
                        if (bp1_v, bp_next_v) not in E and (bp_next_v, bp1_v) not in E:
                            E[(bp1_v, bp_next_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                            E[(bp_next_v, bp1_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                gs_allowed = gs_allowed - 1

        else:
            while gs_allowed:
                ix = i - (gs_allow - gs_allowed + 1)

                if ix < 0 or (sv_cohort_dd_sorted.iloc[ix, 1] != chrom):
                    break
                if sv_cohort_dd_sorted.iloc[ix, 3] == "-":
                    bp_next = ":".join(np.str(j) for j in sv_cohort_dd_sorted.iloc[ix, 1:4])

                    if np.abs(np.int0(sv_cohort_dd_sorted.iloc[ix, 2]) - np.int0(
                            sv_cohort_dd_sorted.iloc[i, 2])) > distance:
                        break
                    else:
                        bp_next_v = BreakPoint(bp_next) if bp_next not in V_dict else V_dict[bp_next]
                        V_dict[bp_next] = bp_next_v

                        if (bp1_v, bp_next_v) not in E and (bp_next_v, bp1_v) not in E:
                            E[(bp1_v, bp_next_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                            E[(bp_next_v, bp1_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                gs_allowed = gs_allowed - 1


def find_genomicEdges_relax(sv_cohort_sort, V_dict, V, E, gs_allow=2, distance=5000000):
    new_order = ["SN", "chr_2", "pos_2", "flag_2", "chr_1", "pos_1", "flag_1"]
    sv_cohort_sort = sv_cohort_sort[["SN", "chr_1", "pos_1", "flag_1", "chr_2", "pos_2", "flag_2"]]
    sv_cohort_sort_rev = sv_cohort_sort[new_order]
    sv_cohort_sort_rev.columns = sv_cohort_sort.columns

    frames = [sv_cohort_sort, sv_cohort_sort_rev]
    sv_cohort_dd = pd.concat(frames)
    sv_cohort_dd_sorted = sv_cohort_dd.sort_values(['SN', 'chr_1', 'pos_1'], ascending=[True, True, True])

    ## find nearby geomic edges
    n_rows = sv_cohort_dd_sorted.shape[0]
    sv_cohort_dd_sorted.index = np.arange(n_rows)

    for i, line in sv_cohort_dd_sorted.iterrows():
        gs_allowed = gs_allow
        bp1, bp2 = combo(line)

        bp1_v = BreakPoint(bp1) if bp1 not in V_dict else V_dict[bp1]
        bp2_v = BreakPoint(bp2) if bp2 not in V_dict else V_dict[bp2]

        V_dict[bp1] = bp1_v
        V_dict[bp2] = bp2_v

        V.append(bp1_v)
        V.append(bp2_v)

        edge1 = BreakEdge(bp1_v, bp2_v, "sv")  # SV edge
        E[(bp1_v, bp2_v)] = edge1

        # search for genomic edge for bp1
        chrom, start, flag = bp1.split(":")
        start = np.int(start)
        chrom = np.int(chrom)

        ix = i
        if flag == "-":
            while ix:
                ix = ix + 1
                if (ix >= n_rows) or (sv_cohort_dd_sorted.iloc[ix, 1] != chrom):  ##
                    break
                if sv_cohort_dd_sorted.iloc[ix, 3] == "+":
                    bp_next = ":".join(np.str(j) for j in sv_cohort_dd_sorted.iloc[ix, 1:4])

                    if np.abs(np.int0(sv_cohort_dd_sorted.iloc[ix, 2]) - np.int0(
                            sv_cohort_dd_sorted.iloc[i, 2])) > distance:
                        break
                    else:
                        bp_next_v = BreakPoint(bp_next) if bp_next not in V_dict else V_dict[bp_next]
                        V_dict[bp_next] = bp_next_v
                        if (bp1_v, bp_next_v) not in E and (bp_next_v, bp1_v) not in E:
                            E[(bp1_v, bp_next_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                            E[(bp_next_v, bp1_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                        break

        else:
            while ix:
                ix = ix - 1

                if ix < 0 or (sv_cohort_dd_sorted.iloc[ix, 1] != chrom):
                    break
                if sv_cohort_dd_sorted.iloc[ix, 3] == "-":
                    bp_next = ":".join(np.str(j) for j in sv_cohort_dd_sorted.iloc[ix, 1:4])

                    if np.abs(np.int0(sv_cohort_dd_sorted.iloc[ix, 2]) - np.int0(
                            sv_cohort_dd_sorted.iloc[i, 2])) > distance:
                        break
                    else:
                        bp_next_v = BreakPoint(bp_next) if bp_next not in V_dict else V_dict[bp_next]
                        V_dict[bp_next] = bp_next_v

                        if (bp1_v, bp_next_v) not in E and (bp_next_v, bp1_v) not in E:
                            E[(bp1_v, bp_next_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                            E[(bp_next_v, bp1_v)] = BreakEdge(bp1_v, bp_next_v, "gs")
                        break





if __name__ == "__main__":
    print("zhouyong")
