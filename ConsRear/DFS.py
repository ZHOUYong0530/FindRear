import numpy as np
import utils


def DFS_visit(sv_graph, u, path, results, E):
    path = path + [u]
    if len(path) >= 2:
        u_prior = path[-2]
    if len(path) == 2 and E[(u_prior, u)].get_status() == "gs":
        return None

    u.color = "gray"
    if u not in sv_graph:
        print(u.get_info())
    for v in sv_graph[u]:
        if v.color == "white":
            v.pi = u
            if len(path) >= 3 and E[(u_prior, u)].get_status() == E[(u, v)].get_status():
                return None
            newpath = DFS_visit(sv_graph, v, path, results, E)
            if newpath and len(newpath) >= 4 and len(newpath) % 2 == 0:
                if not results:
                    results.append(newpath)
                else:
                    flag = "non_duplicated"
                    for res in results:
                        if res[0:len(newpath)] == newpath:
                            flag = "duplicated"
                    if flag == "non_duplicated":
                        results.append(newpath)

    u.color = "black"
    return path


def get_TI_info(bpss, V_dict):
    TI = []
    ##print(len(line["path"].split(" => "))-1)
    for i in np.arange(1, len(bpss) - 1, 2):
        bp1 = bpss[i]
        bp2 = bpss[i + 1]

        bp1_Lcov, bp1_Rcov = utils.cov_median(V_dict[bp1].bp_cov)
        bp2_Lcov, bp2_Rcov = utils.cov_median(V_dict[bp2].bp_cov)

        chrom1, pos1, flag1 = bp1.split(":")
        chrom2, pos2, flag2 = bp2.split(":")

        dist = abs(np.int(pos2) - np.int(pos1))

        if dist < 20000:
            shift = min(np.int(dist / 200), 100)
            TI_cov = np.median(V_dict[bp1].bp_cov[100:101 + shift]) if flag1 == "-" else np.median(
                V_dict[bp2].bp_cov[100:101 + shift])

        elif dist < 100000:
            if flag1 == "-":
                TI_cov = (np.median(V_dict[bp1].bp_cov[101:]) + np.median(V_dict[bp2].bp_cov[0:100])) / 2
            elif flag1 == "+":
                TI_cov = (np.median(V_dict[bp1].bp_cov[0:100]) + np.median(V_dict[bp2].bp_cov[101:])) / 2
        else:
            shift = min(np.int(dist / 10000), 100)
            if flag1 == "-":
                TI_cov = (np.median(V_dict[bp1].kb_cov[101:101 + shift]) + np.median(
                    V_dict[bp2].kb_cov[100 - shift:101])) / 2
            elif flag1 == "+":
                TI_cov = (np.median(V_dict[bp1].kb_cov[100 - shift:101]) + np.median(
                    V_dict[bp2].kb_cov[101:101 + shift])) / 2

        if flag1 == "-":
            TI.append([dist, bp1_Lcov, TI_cov, bp2_Rcov])
        else:
            TI.append([dist, bp2_Lcov, TI_cov, bp1_Rcov])

    cov_space = 0
    if True:  ## output coverage of start and end breakpoint coverage for use
        bp1 = bpss[0]
        bp2 = bpss[-1]

        chrom1, pos1, flag1 = bp1.split(":");
        pos1 = np.int(pos1)
        chrom2, pos2, flag2 = bp2.split(":");
        pos2 = np.int(pos2)

        dist = abs(pos1 - pos2)

        if chrom1 == chrom2 and dist < 10000000:

            if dist < 20000:
                shift = min(np.int(dist / 200), 100)
                cov_space = np.median(V_dict[bp1].bp_cov[100:101 + shift]) if pos1 < pos2 else np.median(
                    V_dict[bp2].bp_cov[100:101 + shift])

            elif dist < 100000:
                if pos1 < pos2:
                    cov_space = (np.median(V_dict[bp1].bp_cov[101:]) + np.median(V_dict[bp2].bp_cov[0:100])) / 2
                else:
                    cov_space = (np.median(V_dict[bp1].bp_cov[0:100]) + np.median(V_dict[bp2].bp_cov[101:])) / 2
            else:
                shift = min(np.int(dist / 10000), 100)
                if pos1 < pos2:
                    cov_space = (np.median(V_dict[bp1].kb_cov[101:101 + shift]) + np.median(
                        V_dict[bp2].kb_cov[100 - shift:101])) / 2
                else:
                    cov_space = (np.median(V_dict[bp1].kb_cov[100 - shift:101]) + np.median(
                        V_dict[bp2].kb_cov[101:101 + shift])) / 2

    return TI, np.around(cov_space, 3)


def DFS_path(sv_graph, u, des, path, results, E):
    path = path + [u]
    if len(path) >= 2:
        u_prior = path[-2]
    if len(path) == 2 and E[(u_prior, u)].get_status() == "gs":
        return None

    if u == des:
        return results.append(path)
    u.color = "gray"
    if u not in sv_graph:
        print(u.get_info())
    for v in sv_graph[u]:
        ##print(u.get_info(),v.get_info(),E[(u,v)].status)
        if v.color == "white":
            v.pi = u
            if len(path) >= 3 and E[(u_prior, u)].get_status() == E[(u, v)].get_status():
                continue
            DFS_path(sv_graph, v, des, path, results,E)
    u.color = "black"


if __name__ == "__main__":
    print("hello world, DFS.py")