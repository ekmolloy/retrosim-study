"""
This file is used to prepare a set of trees (same topology) for plotting.

Copyright (c) 2021 Erin K. Molloy
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import numpy
import os
import sys
import treeswift


def parse_branch_label(label, prefix=None):
    """
    Parse ASTRAL's branch label with quantities:
    q1,q2,q3 = normalized quartet frequencies
               for AB|CD, AD|BC, AC|BD
    f1,f2,f3 = quartet frequencies for AB|CD, AD|BC, AC|BD
    pp1,pp2,pp3 - local pp for AB|CD, AD|BC, AC|BD
    QC = total number of quartets around the branch
         i.e. if A, B, C, D denote the sets of taxa
         in the quadpartition defined by branch,
         then QC = |A| * |B| * |C| * |D|.
    EN = effective number of gene trees with information
         about the resolution of the branch
    """
    data = {}
    for x in label[1:-1].split(';'):
        key, value = x.split('=')
        if prefix is not None:
            key = prefix + '_' + key
        try:
            data[key] = float(value)
        except ValueError:
            data[key] = value
    return data


def read_name_map(input):
    nmap = {}
    with open(input, 'r') as f:
        for line in f:
            code, name = line.strip().split(',')
            nmap[code] = name
    return nmap


def initalize_clade_data(tree):
    clade_data = {}
    i = 0
    for node in tree.traverse_preorder():
        if (not node.is_leaf()):
            labels = sorted([l.label for l in node.traverse_leaves()])
            clade = ",".join(labels)

            clade_data[clade] = {}
            clade_data[clade]["TRUE_BRLN"] = node.edge_length
            clade_data[clade]["BRLN"] = []  # AA_RI_MLEBL
            clade_data[clade]["SUPP"] = []  # AA_pp1
            clade_data[clade]["Q1"] = []    # AA_q1
            clade_data[clade]["Q2"] = []    # AA_q2
            clade_data[clade]["Q3"] = []    # AA_q3
            clade_data[clade]["QC"] = []    # AA_QC
            clade_data[clade]["EN"] = []    # AA_EN

            i = i + 1

    return clade_data


def compute_clade_data(clade_data, trees):
    nt = len(trees)
    for tree in trees:
        for node in tree.traverse_preorder():
            if node.is_leaf():
                continue

            labels = sorted([l.label for l in node.traverse_leaves()])
            clade = ",".join(labels)

            in_true_tree = True
            try:
                data = clade_data[clade]
            except KeyError:
                in_true_tree = False

            if not in_true_tree:
                continue

            brln = None

            if node.label is not None:
                if node.label.find("=") == -1:
                    # Treat label as branch support
                    clade_data[clade]["SUPP"].append(float(node.label))
                else:
                    # Parse ASTRAL_BP label
                    data = parse_branch_label(node.label)

                    brln = data["AA_RI_MLEBL"]
                    clade_data[clade]["SUPP"].append(data["AA_pp1"])
                    clade_data[clade]["Q1"].append(data["AA_q1"])
                    clade_data[clade]["Q2"].append(data["AA_q2"])
                    clade_data[clade]["Q3"].append(data["AA_q3"])
                    clade_data[clade]["QC"].append(data["AA_QC"])
                    clade_data[clade]["EN"].append(data["AA_EN"])

                    is_basic = False

            if brln is None:
                brln = node.edge_length
            
            if brln is not None:
                if numpy.isinf(brln):
                    brln = 9.0
                clade_data[clade]["BRLN"].append(brln)

    for clade in clade_data.keys():
        data = clade_data[clade]

        info = {}
        info["TRUE_BRLN"] = clade_data[clade]["TRUE_BRLN"]
        for key, val in data.items():
            if key == "TRUE_BRLN":
                continue

            n = len(val)

            info[key] = {}
            info[key]["N"] = n

            if n == 0:
                info[key]["MEAN"] = None
            else:
                info[key]["MEAN"] = numpy.mean(val)
                info[key]["STD"] = numpy.std(val)
                info[key]["MEDIAN"] = numpy.median(val)
                info[key]["MIN"] = numpy.min(val)
                info[key]["MAX"] = numpy.max(val)

        clade_data[clade] = info

    return clade_data


if __name__ == "__main__":
    # Process inputs...
    tdir = "../../../data/model-trees/"
    edir = "../../../data/astral-bp-trees/Palaeognathae/"

    fnam = edir + "name-map-to-plot.csv"
    nmap = read_name_map(fnam)

    fnam = tdir + "Palaeognathae_model.tre"
    true_tree = treeswift.read_tree_newick(fnam)
    for node in true_tree.traverse_leaves():
        node.label = nmap[node.label]

    for nret in [1000, 5000]:
        esti_trees = []
        for repl in range(1, 26):
            fnam = edir + "astral-bp-Rep" + str(repl) + "-" + str(nret) + "-rooted-t2.tre"
            esti_tree = treeswift.read_tree_newick(fnam)
            for node in esti_tree.traverse_leaves():
                node.label = nmap[node.label]
            esti_trees.append(esti_tree)

        # Get clade data...
        clade_data = initalize_clade_data(true_tree)
        compute_clade_data(clade_data, esti_trees)

        with open("sim-" + str(nret) + "-table.txt", 'w') as f:
            f.write("\\begin{table}\n")
            f.write("\\centering\n")
            f.write("\\caption{Add caption..." + str(nret) + "RIs}\n")
            f.write("\\begin{tabular}{cccccc}\n")
            f.write("\\toprule\n")
            f.write("Clade & $\\sigma^*$ & $\\hat{\\sigma}$ & local PP & $EN$ & \\#STwC\\\\\n")
            f.write("\\midrule\n")
            for clade, data in clade_data.items():
                f.write("%s" % clade)

                if data["SUPP"]["N"] == 0:
                    continue

                val = round(round(data["TRUE_BRLN"], 5), 4)
                f.write(" & $%1.4f$" % val)

                avg = round(round(data["BRLN"]["MEAN"], 5), 4)
                std = round(round(data["BRLN"]["STD"], 5), 4)
                f.write(" & $%1.4f \\pm %1.4f$" % (avg, std))

                avg = round(round(data["SUPP"]["MEAN"], 3), 2)
                std = round(round(data["SUPP"]["STD"], 3), 2)
                f.write(" & $%1.2f \\pm %1.2f$" % (avg, std))

                avg = round(round(data["EN"]["MEAN"], 3), 2)
                std = round(round(data["EN"]["STD"], 3), 2)
                f.write(" & $%1.2f \\pm %1.2f$" % (avg, std))

                f.write(" & $%d$" % data["SUPP"]["N"])

                f.write("\\\\\n")

            f.write("\\bottomrule\n")
            f.write("\\end{tabular}\n")
            f.write("\\end{table}\n")

