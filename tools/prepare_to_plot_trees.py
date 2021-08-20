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
        labels = sorted([l.label for l in node.traverse_leaves()])
        clade = ",".join(labels)
        if len(labels) > 1:
            clade_data[clade] = {}
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
        for key, val in data.items():
            n = len(val)

            info[key] = {}
            info[key]["N"] = n

            if n == 0:
                info[key]["MEAN"] = None
            elif n == nt:
                info[key]["MEAN"] = numpy.mean(val)
                info[key]["STD"] = numpy.std(val)
                info[key]["MEDIAN"] = numpy.median(val)
                info[key]["MIN"] = numpy.min(val)
                info[key]["MAX"] = numpy.max(val)
            else:
                sys.exit("Warning: Not all trees have information on branch!\n")

        clade_data[clade] = info

    return clade_data


def annotate_with_clade_data(tree, clade_data):
    for node in tree.traverse_preorder():
        labels = sorted([l.label for l in node.traverse_leaves()])
        if len(labels) > 1:
            clade = ",".join(labels)
            data = clade_data[clade]
            node.edge_length = data["BRLN"]["MEAN"]
            node.label = data["SUPP"]["MEAN"]


def annotate_with_treelist(main_tree, trees):
    clade_data = initalize_clade_data(main_tree)
    compute_clade_data(clade_data, trees)
    annotate_with_clade_data(main_tree, clade_data)
    return clade_data


def clean_branch_info(tree):
    for node in tree.traverse_preorder():
        node.label = None
        node.edge_length = None


def main(args):
    if args.namemap is not None:
        nmap = read_name_map(args.namemap)
    else:
        nmap = None

    trees =[]
    for tf in args.trees:
        tree = treeswift.read_tree_newick(tf)
        if nmap is not None:
            for node in tree.traverse_leaves():
                node.label = nmap[node.label]
        trees.append(tree)
    main_tree = trees[0]

    if args.topo_only:
        clean_branch_info(main_tree)
    else:
        # Annotate with median branch lengths and 
        # median support values from tree list
        clade_data = annotate_with_treelist(main_tree, trees)

        # Write clade information
        with open(args.output + ".clade_data", 'w') as f:
            for clade, data in clade_data.items():
                f.write("clade = %s\n" % clade)
                for key, val in data.items():
                    xn = val["N"]
                    if xn > 0:
                        xavg = round(round(val["MEAN"], 5), 4)
                        xstd = round(round(val["STD"], 5), 4)
                        xmin = round(round(val["MIN"], 5), 4)
                        xmax = round(round(val["MAX"], 5), 4)
                        f.write("    %s : %1.4f +/- %1.4f, %d elements in range [%1.4f, %1.4f]\n" % (key, xavg, xstd, xn, xmin, xmax))
            f.write("\n")

        # Clean-up labels
        for node in main_tree.traverse_preorder():
            if node.is_leaf():

                continue
            if node.label is not None:
                supp = node.label * 100
                supp = round(round(supp, 1), 0)
                node.label = str("%d" % supp)

            if node.edge_length is not None:
                brln = round(round(node.edge_length, 5), 4)
                node.edge_length = brln

        # Set terminal edges to 0.1
        tlen = 0.1
        for node in main_tree.traverse_leaves():
            node.edge_length = tlen

        is_outg_rooted = False
        children = main_tree.root.child_nodes()
        for child in children:
            leaves = [l for l in child.traverse_leaves()]
            if len(leaves) == 1:
                is_outg_rooted = True

        if is_outg_rooted:
            for child in children:
                child.edge_length = tlen / 2.0

    # Write tree
    with open(args.output + ".tre", 'w') as f:
        f.write(main_tree.newick() + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--trees", type=str, nargs='+',
                        help="Input list of trees to compute median branch lengths / support values",
                        required=True)

    parser.add_argument("-m", "--namemap", type=str,
                        help="Name map CSV : old label,new label",
                        required=False)

    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        required=True)

    parser.add_argument("--topo_only", action="store_true")

    main(parser.parse_args())
