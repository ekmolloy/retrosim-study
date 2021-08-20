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


def extract_clade_set(tree):
    clades = set([])
    i = 0
    for node in tree.traverse_preorder():
        leaf_labels = sorted([l.label for l in node.traverse_leaves()])
        clades.add(",".join(leaf_labels))
    return clades


def is_outgroup_rooted(tree):
    children = tree.root.child_nodes()
    for child in children:
        leaves = [l.label for l in child.traverse_leaves()]
        if len(leaves) == 1:
            return leaves[0]
    return None


def main(args):
    true_tree = treeswift.read_tree_newick(args.true)

    true_og = is_outgroup_rooted(true_tree)
    if true_og is None:
        sys.exit("True tree is NOT rooted at an outgroup!\n")

    true_clades = extract_clade_set(true_tree)

    supps =[]
    brlns = []
    for et in args.esti:
        esti_tree = treeswift.read_tree_newick(et)

        esti_og = is_outgroup_rooted(esti_tree)
        if esti_og is None:
            sys.exit("Estimated tree is NOT rooted at an outgroup!\n")
        if esti_og != true_og:
            sys.exit("True and estimated trees have different outgroups!\n")

        for node in esti_tree.traverse_preorder():
            if (node.label is not None) and (not node.is_leaf()):
                supp = float(node.label)
                brln = node.edge_length

                leaf_labels = sorted([l.label for l in node.traverse_leaves()])
                esti_clade = ",".join(leaf_labels)

                if esti_clade in true_clades:
                    pass
                else:
                    supps.append(supp)
                    brlns.append(brln)

    # Find max local pp
    supps = numpy.array(supps)
    sys.stdout.write("Local PP values:\n")
    sys.stdout.write("  Avg %f\n" % numpy.mean(supps))
    sys.stdout.write("  Std %f\n" % numpy.std(supps))
    sys.stdout.write("  Min %f\n" % numpy.min(supps))
    sys.stdout.write("  Max %f\n\n" % numpy.max(supps))

    # Find max branch length
    brlns = numpy.array(brlns)
    sys.stdout.write("Branch Lengths:\n")
    sys.stdout.write("  Avg %f\n" % numpy.mean(brlns))
    sys.stdout.write("  Std %f\n" % numpy.std(brlns))
    sys.stdout.write("  Min %f\n" % numpy.min(brlns))
    sys.stdout.write("  Max %f\n" % numpy.max(brlns))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--true", type=str,
                        help="Input true tree",
                        required=True)

    parser.add_argument("-e", "--esti", type=str, nargs='+',
                        help="Input list of estimated trees",
                        required=True)

    main(parser.parse_args())
