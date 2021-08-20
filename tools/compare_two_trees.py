"""
This file is used compare two phylogenetic trees (newick strings).

Copyright (c) 2020 Erin K. Molloy
All rights reserved.

License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives
import os
import sys


def compare_trees(tr1, tr2):
    """
    Compares two trees

    Parameters
    ----------
    tr1 : dendropy tree object
            First tree (typically the model tree)
    tr2 : dendropy tree object
            Second tree (typically the estimated tree)

    Returns
    -------
    nl : int
         Size of the shared leaf set, i.e., the number of leaves in both trees
    i1 : int
          Number of internal edges in tree 1 after restriction to shared leaves
    i2 : int
          Number of internal edges in tree 2 after restriction to shared leaves
    fn : int
         Number of edges in tree 1 that are not in tree 2
    fp : int
         Number of edges in tree 2 that are not in tree 1
    rf : float
         Normalized Robinson-Foulds (RF) distance between tree 1 and 2

    Example
    -------
    If tree 1 corresponds to "(((A,B,C),D),E);"
    and tree 2 corresponds to "((((A,B),C),D),E);",
    then the output is "5 1 2 0 1 0.25". In this example,
      + tree 1 and tree 2 share five leaves (A, B, C, D, E).
      + tree 1 has one internal edge "A,B,C|D,E"
      + tree 2 has two internal edges "A,B|C,D,E" and "A,B,C|D,E"
      + no edges in the tree 1 are missing from tree 2
      + one edge in the tree 2 is missing from the tree 1
      + normalized RF distance is (FP+FN)/(2*NL-6) = (1+0)/(2*5-6) = 0.25
    """

    # Unroot the two trees!
    tr1.is_rooted = False
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    tr2.is_rooted = False
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    # Restrict the two trees to the same leaf set if necessary!
    lb1 = set([leaf.taxon.label for leaf in tr1.leaf_nodes()])
    lb2 = set([leaf.taxon.label for leaf in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    # Compare trees!
    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    i1 = len(tr1.internal_edges(exclude_seed_edge=True))
    i2 = len(tr2.internal_edges(exclude_seed_edge=True))
    [fp, fn] = false_positives_and_negatives(tr1, tr2)

    # Compute normalized values
    if i1 != nl - 3:
        print("WARNING: True tree is not binary!")

    fpr = fp / i1  # Divided by # of internal branches in true species tree!
    fnr = fn / i1  # Divided by # of internal branches in true species tree!
    rf = (fn + fp) / (i1 + i1)

    return(nl, i1, i2, fn, fp, fnr, fpr, rf)


def main(args):
    taxa = dendropy.TaxonNamespace()

    truetree = dendropy.Tree.get(path=args.truetree,
                                 schema='newick',
                                 rooting='force-unrooted',
                                 taxon_namespace=taxa)

    estitree = dendropy.Tree.get(path=args.estitree,
                                 schema='newick',
                                 rooting='force-unrooted',
                                 taxon_namespace=taxa)

    if args.threshold is not None:
        for node in estitree.postorder_node_iter():
            if node.edge.length is None:
                node.edge.length = 1.0
            else:
                node.edge.length = 1.0
            if node.label is not None:
                support = float(node.label)
                if support < args.threshold:
                    node.edge.length = None
        estitree.collapse_unweighted_edges()

    [nl, i1, i2, fn, fp, fnr, fpr, rf] = compare_trees(truetree, estitree)
    sys.stdout.write('%d,%d,%d,%d,%d,%1.6f,%1.6f,%1.6f\n' %
                     (nl, i1, i2, fn, fp, fnr, fpr, rf))
    sys.stdout.flush()
    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--truetree", type=str,
                        help="Input file containing true tree",
                        required=True)

    parser.add_argument("-e", "--estitree", type=str,
                        help="Input file containing estimated tree",
                        required=True)

    parser.add_argument("-x", "--threshold", type=float,
                        help="Threshold for collapsing branches "
                              "in estimated tree",
                        required=False)

    main(parser.parse_args())
