"""
This file is used to prepare a set of trees (same topology) for plotting.

Copyright (c) 2021 Erin K. Molloy
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import copy
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives
import numpy
import os
import sys


def get_clade_set(tree):
    clades = []
    for node in tree.nodes():
        leaves = sorted([l.taxon.label for l in node.leaf_nodes()])
        clades.append(",".join(leaves))
    return set(clades)


def is_different(tree1, tree2):
    clades1 = get_clade_set(tree1)
    clades2 = get_clade_set(tree2)
    if clades1 == clades2:
        return False
    return True


def read_name_map(input):
    nmap = {}
    with open(input, 'r') as f:
        for line in f:
            code, name = line.strip().split(',')
            nmap[code] = name
    return nmap


def clean_branch_info(tree):
    for node in tree.nodes():
        node.label = None
        node.edge.length = None


def is_binary(tree):
    """
    Checks if a tree is binary

    Parameters
    ----------
    tree : dendropy tree object

    Returns
    -------
    True if the tree is binary and False otherwise
    """
    nodes = [n for n in tree.preorder_node_iter()]
    for n in nodes[1:]:
        if not n.is_leaf():
            children = n.child_nodes()
            if len(children) != 2:
                return False
    return True


def is_rooted(tree):
    """
    Checks if a tree is rooted

    Parameters
    ----------
    tree : dendropy tree object

    Returns
    -------
    True if the tree is rooted and False otherwise
    """
    n = tree.seed_node
    children = n.child_nodes()
    if len(children) != 2:
        return False
    return True


def initialize_terminal_branch_lengths(tree):
    # Set terminal edges to 0.1
    for node in tree.leaf_nodes():
        node.edge.length = 0.1

    is_outg_rooted = False
    children = tree.seed_node.child_nodes()
    for child in children:
        leaves = [l for l in child.leaf_nodes()]
        if len(leaves) == 1:
            is_outg_rooted = True

    if is_outg_rooted:
        for child in children:
            child.edge.length = 0.05


def force_ultrametric(tree):
    
    if not is_binary(tree):
        sys.exit("Tree is not binary!")
    if not is_rooted(tree):
        sys.exit("Tree is not rooted!")

    dists = []
    tree.calc_node_root_distances()
    for leaf in tree.leaf_nodes():
        dists.append(leaf.root_distance)
    maxd = numpy.ceil(numpy.max(dists))

    for leaf in tree.leaf_nodes():
        dist = leaf.parent_node.root_distance
        leaf.edge.length = maxd - dist 


def relabel_leaves(tree, nmap):
    for node in tree.leaf_nodes():
        node.taxon.label = nmap[node.taxon.label]


def main(args):
    # Process species tree
    taxa = dendropy.TaxonNamespace()
    stree = dendropy.Tree.get(path=args.stree,
                             schema='newick',
                             rooting='force-rooted',
                             preserve_underscores=True)

    if args.namemap is not None:
        nmap = read_name_map(args.namemap)
        relabel_leaves(stree, nmap)

    for node in stree.nodes():
        node.label = None
        if node.edge.length is not None:
            if numpy.isinf(node.edge.length):
                node.edge.length = 9.0

    initialize_terminal_branch_lengths(stree)
    force_ultrametric(stree)

    st = stree.as_string(schema="newick", unquoted_underscores=True)[5:]

    # Proces gene trees
    clean_branch_info(stree)
    gtrees = [stree]
    gts = [stree.as_string(schema="newick", unquoted_underscores=True)[5:]]
    if args.gtrees is not None:
        for gf in args.gtrees:
            gtree = dendropy.Tree.get(path=gf,
                                      schema='newick',
                                      rooting='force-rooted',
                                      preserve_underscores=True)

            clean_branch_info(gtree)
            gtrees.append(gtree)
            gts.append(gtree.as_string(schema="newick", unquoted_underscores=True)[5:])

    # Write nexus file
    with open(args.output, 'w') as f:
        f.write("#NEXUS\n")
        f.write("BEGIN TREES;\n")
        for i, gt in enumerate(gts):
            f.write("  tree gt%d = %s\n" % (i, gt))
        f.write("END;\n\n")
        f.write("BEGIN NETWORKS;\n")
        f.write("  Network net1 = %s\n" % st)
        f.write("END;\n\n")
        f.write("BEGIN PHYLONET;\n")
        for i, gt in enumerate(gts):
            if i == 0:
                f.write("CalGTProb net1 (gt%d);\n" % i)
            else:
                if is_different(stree, gtrees[i]):
                    f.write("CalGTProb net1 (gt%d);\n" % i)
        f.write("END;\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--stree", type=str,
                        help="Model species tree",
                        required=True)

    parser.add_argument("-g", "--gtrees", type=str, nargs='+',
                        help="Gene tree topologies",
                        required=False)

    parser.add_argument("-m", "--namemap", type=str,
                        help="Name map CSV : stree label,gtree label",
                        required=False)

    parser.add_argument("-o", "--output", type=str,
                        help="Output file",
                        required=True)

    main(parser.parse_args())
