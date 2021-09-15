import graph
import utils
import DFS
import pandas as pd
import os
import numpy as np
import argparse


def con_graph(test_sn_mark_rm):
    # fill V and Ｅ
    V_dict = {}
    V = []
    E = {}
    graph.find_genomicEdges_relax(test_sn_mark_rm, V_dict, V, E)
    # check V E success or not
    # def check_V_E(V,E,V_dict):
    for u, v in E:
        # print(u.get_info(),v.get_info(),E[(u,v)].get_status())
        if u not in V or v not in V:
            print("erro")

    E_vector = []
    E_vector1 = []
    E_vector2 = []
    for u, v in E:
        E_vector.append(u)
        E_vector.append(v)
        E_vector1.append(u)
        E_vector2.append(v)

    for v in V:
        if v not in E_vector:
            print("error: v not in E")
    print(set(E_vector).difference(set(V)))
    print(set(V).difference(set(E_vector)))

    for info in V_dict:
        if V_dict[info] not in V:
            print("error v not in V_dict()")

    # construct sv
    sv_graph = {}
    for v1, v2 in E:
        if v1 not in sv_graph:
            sv_graph[v1] = []
        sv_graph[v1].append(v2)

    return V, E, V_dict, sv_graph


def fd_search(input_file, output_file, cbp_dict, ckp_dict):
    sv_cohort = pd.read_csv(input_file, sep="\t")
    fo = open(output_file, "w")
    for sn in set(sv_cohort["SN"]):
        test_sn = sv_cohort[(sv_cohort["SN"] == sn) & (sv_cohort["calling_type"] != "TSI_L")][
            ["SN", "chr_1", "pos_1", "flag_1", "chr_2", "pos_2", "flag_2"]]
        # data cleaning
        test_sn_mark = graph.cohort_rm_simple(test_sn)
        test_sn_mark_rm = test_sn_mark.query('simple == "non-simple"')

        # construct graph
        V, E, V_dict, sv_graph = con_graph(test_sn_mark_rm)

        # initialize coverage of breakpoints
        # check
        if sn not in cbp_dict or sn not in ckp_dict:
            print("the coverge of " + sn + " not in covergee file list")
            break
        coverage_outfile_kb = ckp_dict[sn]
        coverage_outfile_bp = cbp_dict[sn]
        if os.path.exists(coverage_outfile_kb) and os.path.exists(coverage_outfile_bp):
            print(coverage_outfile_kb, coverage_outfile_bp)
            cov_200 = pd.read_csv(coverage_outfile_bp, sep="\t", index_col=0)
            cov_20kb = pd.read_csv(coverage_outfile_kb, sep="\t", index_col=0)
            for v in V:
                if v.get_info() not in cov_200.columns or v.get_info() not in cov_20kb.columns:
                    print("error happens: vertices not in coverage file")
                    break
                v.bp_cov = cov_200[v.get_info()]
                v.kb_cov = cov_20kb[v.get_info()]
        else:
            print("{} coverage file not exist".format(sn))


        # identify fold-back inversion
        test_sn_mark_rm_sort_dd = graph.cohort_sort_dd(test_sn_mark_rm)
        bpss_list = {}
        for i in np.arange(test_sn_mark_rm_sort_dd.shape[0] - 1):
            line1 = test_sn_mark_rm_sort_dd.iloc[i]
            line2 = test_sn_mark_rm_sort_dd.iloc[i + 1,]
            bp11, bp12 = graph.combo(line1)
            bp21, bp22 = graph.combo(line2)

            chr11, start11, flag11 = bp11.split(":")
            start11 = np.int(start11)
            chr21, start21, flag21 = bp21.split(":")
            start21 = np.int(start21)

            if bp12 == bp21:
                continue
            elif flag11 == flag21 and chr11 == chr21:
                for u in V:
                    u.color = "white"
                    u.pi = ""
                results = []
                start = V_dict[bp11]
                des = V_dict[bp21]
                DFS.DFS_path(sv_graph, start, des, [], results, E)
                for p in results:
                    if len(p) % 2 == 0 and len(p) >= 4:
                        s_left, s_right = utils.cov_median(p[0].bp_cov)
                        e_left, e_right = utils.cov_median(p[-1].bp_cov)

                        if flag11 == "+" and s_left > s_right + 0.5 and e_left > e_right + 0.5:
                            # and (start21 - start11) < 100000:
                            print(sn, "fb_confident", bp11, bp21, sep="\t", end="\t", file=fo)
                        elif flag11 == "-" and e_right > e_left + 0.5 and s_right > s_left + 0.5:
                            # and (start21 - start11) < 100000:
                            print(sn, "fb_confident", bp11, bp21, sep="\t", end="\t", file=fo)
                        else:
                            # continue
                            print(sn, "fb_unknown", bp11, bp21, sep="\t", end="\t", file=fo)
                        for ii in p:
                            i_left = np.around(np.median(ii.bp_cov[0:100]), 3)
                            i_right = np.around(np.median(ii.bp_cov[101:201]), 3)
                            print(ii.get_info(), i_left, i_right, sep=";", end=" => ", file=fo)

                        if bp11 not in bpss_list:
                            bpss_list[bp11] = []
                        bpss_list[bp11].append([np.str(j.get_info()) for j in p])
                        print("\n", end="", file=fo)

    fo.close()


# identify complex TD/DEL events
def tandem_del_search(input_file, output_file, cbp_dict, ckp_dict):
    sv_cohort = pd.read_csv(input_file, sep="\t")
    fo = open(output_file, "w")
    for sn in sv_cohort["SN"].sort_values().unique():
        test_sn = sv_cohort[(sv_cohort["SN"] == sn) & (sv_cohort["calling_type"] != "TSI_L")][
            ["SN", "chr_1", "pos_1", "flag_1", "chr_2", "pos_2", "flag_2"]]
        # data cleaning
        test_sn_mark = graph.cohort_rm_simple(test_sn)
        test_sn_mark_rm = test_sn_mark.query('simple == "non-simple"')
        # construct graph
        V, E, V_dict, sv_graph = con_graph(test_sn_mark_rm)

        # initialize coverage of breakpoints
        if sn not in cbp_dict or sn not in ckp_dict:
            print("the coverge of " + sn + " not in covergee file list")
            break
        coverage_outfile_kb = ckp_dict[sn]
        coverage_outfile_bp = cbp_dict[sn]
        if os.path.exists(coverage_outfile_kb) and os.path.exists(coverage_outfile_bp):
            print(coverage_outfile_kb, coverage_outfile_bp)
            cov_200 = pd.read_csv(coverage_outfile_bp, sep="\t", index_col=0)
            cov_20kb = pd.read_csv(coverage_outfile_kb, sep="\t", index_col=0)
            for v in V:
                if v.get_info() not in cov_200.columns or v.get_info() not in cov_20kb.columns:
                    print("error happens: vertices not in coverage file")
                    break
                v.bp_cov = cov_200[v.get_info()]
                v.kb_cov = cov_20kb[v.get_info()]
        else:
            print("{} coverage file not exist".format(sn))

        # identify inversions (unblanced inversions)
        test_sn_mark_rm_sort_dd = graph.cohort_sort_dd(test_sn_mark_rm)
        bpss_list = {}
        for i in np.arange(test_sn_mark_rm_sort_dd.shape[0] - 1):
            line1 = test_sn_mark_rm_sort_dd.iloc[i]
            line2 = test_sn_mark_rm_sort_dd.iloc[i + 1,]
            bp11, bp12 = graph.combo(line1)
            bp21, bp22 = graph.combo(line2)

            chr11, start11, flag11 = bp11.split(":")
            start11 = np.int(start11)
            chr21, start21, flag21 = bp21.split(":")
            start21 = np.int(start21)

            dist = abs(np.int(start11) - np.int(start21))

            if dist < 20000:
                shift = min(np.int(dist / 200), 100)
                TI_cov = np.median(V_dict[bp11].bp_cov[100:101 + shift])

            elif dist < 100000:
                TI_cov = (np.median(V_dict[bp11].bp_cov[101:]) + np.median(V_dict[bp21].bp_cov[0:100])) / 2

            else:
                shift = min(np.int(dist / 10000), 100)
                TI_cov = (np.median(V_dict[bp11].kb_cov[101:101 + shift]) + np.median(
                    V_dict[bp21].kb_cov[100 - shift:101])) / 2

            if bp12 == bp21:
                continue
            elif bp11.split(":")[2] == "+" and bp21.split(":")[2] == "-" and bp11.split(":")[0] == bp21.split(":")[0]:
                for u in V:
                    u.color = "white"
                    u.pi = ""
                results = []
                start = V_dict[bp11]
                des = V_dict[bp21]
                DFS.DFS_path(sv_graph, start, des, [], results, E)
                for p in results:
                    if len(p) % 2 == 0:
                        s_left, s_right = utils.cov_median(p[0].bp_cov)
                        e_left, e_right = utils.cov_median(p[-1].bp_cov)

                        if (s_left > TI_cov + 0.06 and e_right > TI_cov + 0.06) and TI_cov < 1:

                            print(sn, "del_confident", bp11, bp21, TI_cov, sep="\t", end="\t", file=fo)
                        else:
                            print(sn, "del_unknown", bp11, bp21, TI_cov, sep="\t", end="\t", file=fo)
                        for ii in p:
                            i_left = np.around(np.median(ii.bp_cov[0:100]), 3)
                            i_right = np.around(np.median(ii.bp_cov[101:201]), 3)
                            print(ii.get_info(), i_left, i_right, sep=";", end=" => ", file=fo)

                        # print("\n", file=fo)


            elif bp11.split(":")[2] == "-" and bp21.split(":")[2] == "+" and bp11.split(":")[0] == bp21.split(":")[0]:
                for u in V:
                    u.color = "white"
                    u.pi = ""
                results = []
                start = V_dict[bp11]
                des = V_dict[bp21]
                DFS.DFS_path(sv_graph, start, des, [], results, E)
                for p in results:
                    if len(p) % 2 == 0:
                        s_left, s_right = utils.cov_median(p[0].bp_cov)
                        e_left, e_right = utils.cov_median(p[-1].bp_cov)

                        if s_left + 0.06 < TI_cov and e_right + 0.06 < TI_cov:
                            print(sn, "TD_confident", bp11, bp21, TI_cov, sep="\t", end="\t", file=fo)
                        else:
                            print(sn, "TD_unknown", bp11, bp21, TI_cov, sep="\t", end="\t", file=fo)
                        for i in p:
                            print(i.get_info(), np.around(np.median(i.bp_cov[0:100]), 3),
                                  np.around(np.median(i.bp_cov[101:201]), 3), sep=";", end=" => ", file=fo)
                        print("\n", end="", file=fo)
    fo.close()


def genomic_chain_search(input_file, output_file, cbp_dict, ckp_dict):
    sv_cohort = pd.read_csv(input_file, sep="\t")
    fo = open(output_file, "w")

    print("sn", "status", "se_cov", "se_kb_cov", "interval_cov", "connected_bps", "coverage", sep="\t", file=fo)

    for sn in sv_cohort["SN"].unique():
        test_sn = sv_cohort[(sv_cohort["SN"] == sn) & (sv_cohort['calling_type'] != 'TSI_L')][
            ["SN", "chr_1", "pos_1", "flag_1", "chr_2", "pos_2", "flag_2"]]
        # data cleaning
        test_sn_mark = graph.cohort_rm_simple(test_sn)
        test_sn_mark_rm = test_sn_mark.query('simple == "non-simple"')

        # construct graph
        V, E, V_dict, sv_graph = con_graph(test_sn_mark_rm)

        # inilinize coverage of breakpoints
        if sn not in cbp_dict or sn not in ckp_dict:
            print("the coverge of " + sn + " not in covergee file list")
            break
        coverage_outfile_kb = ckp_dict[sn]
        coverage_outfile_bp = cbp_dict[sn]
        if os.path.exists(coverage_outfile_kb) and os.path.exists(coverage_outfile_bp):
            print(coverage_outfile_kb, coverage_outfile_bp)
            cov_200 = pd.read_csv(coverage_outfile_bp, sep="\t", index_col=0)
            cov_20kb = pd.read_csv(coverage_outfile_kb, sep="\t", index_col=0)
            for v in V:
                if v.get_info() not in cov_200.columns or v.get_info() not in cov_20kb.columns:
                    print("error happens: vertices not in coverage file")
                    break
                v.bp_cov = cov_200[v.get_info()]
                v.kb_cov = cov_20kb[v.get_info()]
        else:
            print("{} coverage file not exist".format(sn))

        # identify genomic chains
        test_sn_mark_rm_sort_dd = graph.cohort_sort_dd(test_sn_mark_rm)
        bpss_list = {}

        for i in np.arange(test_sn_mark_rm_sort_dd.shape[0] - 1):
            line1 = test_sn_mark_rm_sort_dd.iloc[i]
            line2 = test_sn_mark_rm_sort_dd.iloc[i + 1,]
            bp11, bp12 = graph.combo(line1)
            bp21, bp22 = graph.combo(line2)

            chr11, start11, flag11 = bp11.split(":")
            start11 = np.int(start11)
            chr21, start21, flag21 = bp21.split(":")
            start21 = np.int(start21)

            for u in V:
                u.color = "white"
                u.pi = ""
            results = []
            start = V_dict[bp11]
            des = V_dict[bp21]
            DFS.DFS_visit(sv_graph, start, [], results, E)

            for newpath in results:
                bpss = [i.get_info() for i in newpath]

                end = newpath[-1]

                if len(bpss) >= 4:
                    TI_info, cov_space = DFS.get_TI_info(bpss, V_dict)

                    s_left, s_right = utils.cov_median(start.bp_cov)
                    e_left, e_end = utils.cov_median(end.bp_cov)

                    s_kb_left, s_kb_right = utils.cov_median(start.kb_cov)
                    e_kb_left, e_kb_end = utils.cov_median(end.kb_cov)

                    se_cov = ";".join(np.str(j) for j in [s_left, s_right, e_left, e_end])
                    se_kb_cov = ";".join(np.str(j) for j in [s_kb_left, s_kb_right, e_kb_left, e_kb_end])

                    g = lambda x: ":".join([np.str(j) for j in x])
                    TI_info_output = ";".join([g(j) for j in TI_info])

                    good_TI_num = 0
                    for info in TI_info:
                        if info[2] > (info[1] + 0.06) and info[2] > (info[3] + 0.06):
                            # plot_templated_insertion_inline([bpss],sn,outdir)
                            good_TI_num += 1

                    if good_TI_num == len(TI_info):
                        # print(sn,end="\t")
                        print(sn, "confident", se_cov, se_kb_cov, cov_space,
                              ";".join([np.str(i.get_info()) for i in newpath]), TI_info_output, sep="\t", file=fo)
                    if good_TI_num != len(TI_info):
                        # print(sn,end="\t")
                        print(sn, "unknown", se_cov, se_kb_cov, cov_space,
                              ";".join([np.str(i.get_info()) for i in newpath]), TI_info_output, sep="\t", file=fo)
                        # plot_templated_insertion_inline([bpss],sn,outdir)
    fo.close()


def unbal_inver(gc_file, output_file):
    genomic_chains = pd.read_csv(gc_file, sep="\t")
    fo = open(output_file, "w")
    f_counts = 0
    bpss_set = []
    for i, line in genomic_chains.iterrows():
        sn = line.sn
        bpss = line.connected_bps.split(";")

        typ = ""
        if True:  # re.search("confident",line.status):
            orientation_1 = bpss[0].split(":")[-1]
            chrom_1 = bpss[0].split(":")[0]
            start_1 = np.int(bpss[0].split(":")[1])
            orientation_2 = bpss[-1].split(":")[-1]
            chrom_2 = bpss[-1].split(":")[0]
            end_1 = np.int(bpss[-1].split(":")[1])
            s_left, s_right, e_left, e_right = np.array(line["se_cov"].split(";")).astype(np.float64)
            s_kb_left, s_kb_right, e_kb_left, e_kb_right = np.array(line["se_kb_cov"].split(";")).astype(np.float64)

            bp_start = bpss[0]
            bp_end = bpss[-1]

            if set(bpss) not in bpss_set:  # remove duplicates
                bpss_set.append(set(bpss))
                # print("unique","\t".join(map(str,line)),sep="\t")
            else:
                continue

            if chrom_1 == chrom_2:
                if orientation_1 == "+" and orientation_2 == "-":
                    typ = "del" if start_1 < end_1 else "td"
                elif orientation_1 == "-" and orientation_2 == "+":
                    typ = "td" if start_1 < end_1 else "del"
                else:
                    typ = "inv"
            else:
                typ = "trans"

            TI_len = []
            TI_is_in_DEL = False
            for ii in np.arange(1, len(line["connected_bps"].split(";")) - 2, 2):
                pos1 = line["connected_bps"].split(";")[ii].split(":")[1]
                pos2 = line["connected_bps"].split(";")[ii + 1].split(":")[1]

                TI_len.append(abs(np.int(pos2) - np.int(pos1)))

                if np.int(pos1) > min(start_1, end_1) and np.int(pos1) < max(start_1, end_1) \
                        and np.int(pos2) > min(start_1, end_1) and np.int(pos2) < max(start_1, end_1):
                    TI_is_in_DEL = True

            pos1 = line["connected_bps"].split(";")[0].split(":")[1]
            pos2 = line["connected_bps"].split(";")[-1].split(":")[1]

            pos_min = min(np.int(pos1), np.int(pos2))
            pos_max = max(np.int(pos1), np.int(pos2))

            dist = pos_max - pos_min

            if typ == "del" and orientation_1 == "+" and orientation_2 == "-" and TI_is_in_DEL:
                ix_cov1 = s_right < 1 or e_left < 1
                ix_cov2 = line["interval_cov"] < 1
                ix_cov = ix_cov1 or ix_cov2

                if len(bpss) <= 8 and ix_cov and dist < 10000000:
                    print(sn, typ, line.status, line.se_cov, line.interval_cov, line.connected_bps,
                          sep="\t", file=fo)
                    # plot_templated_insertion_inline([line["connected_bps"].split(";")], line.sn)
                    f_counts += 1

            elif typ == "del" and orientation_1 == "-" and orientation_2 == "+" and TI_is_in_DEL:
                ix_cov1 = s_right < 1 and e_left < 1
                ix_cov2 = line["interval_cov"] < 1
                ix_cov = ix_cov1 or ix_cov2
                if len(bpss) <= 6 and ix_cov and dist < 10000000:
                    print(sn, typ, line.status, line.se_cov, line.interval_cov, line.connected_bps,
                          sep="\t", file=fo)
                    # plot_templated_insertion_inline([line["connected_bps"].split(";")],line.sn,outdir)
                    f_counts += 1
    fo.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process sv and coverage files")
    parser.add_argument('sv', help="SV_cort (7 columns)")
    parser.add_argument('out_rear', help="output rearrangements")
    parser.add_argument('cbp', help="200-rows matrix (200bp-unit) list sn link")
    parser.add_argument('ckp', help="200-rows matrix (10kb-unit) list # sn link")

    parser.add_argument("--cov", help="patchwork_converted coverage")
    parser.add_argument("--type", help="fold-back,del_td,gc,default: b-del")
    args = parser.parse_args()

    print(args.sv, args.out_rear, args.cbp, args.ckp)

    cbp_dict = {}
    ckp_dict = {}
    for line in open(args.cbp):
        sn, link = line.strip().split("\t")
        cbp_dict[sn] = link

    for line in open(args.ckp):
        sn, link = line.strip().split("\t")
        ckp_dict[sn] = link

    if args.type == "fold-back":
        fd_search(args.sv, args.out_rear, cbp_dict, ckp_dict)
    elif args.type == "del_td":
        tandem_del_search(args.sv, args.out_rear, cbp_dict, ckp_dict)
    elif args.type == "gc":
        genomic_chain_search(args.sv, args.out_rear, cbp_dict, ckp_dict)
    elif args.type == "inv":
        genomic_chain_search(args.sv, "gc_tmp.txt", cbp_dict, ckp_dict)
        unbal_inver("gc_tmp.txt", args.out_rear)
    else:
        tandem_del_search(args.sv, args.out_rear, cbp_dict, ckp_dict)
