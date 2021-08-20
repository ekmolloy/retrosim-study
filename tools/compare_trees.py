"""
This file is used compare phylogenetic trees (newick strings).

Copyright (c) 2021 Erin K. Molloy
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


def main(args):
    taxa = dendropy.TaxonNamespace()

    trees = []
    for tf in args.trees:
        tree = dendropy.Tree.get(path=tf,
                                 schema='newick',
                                 rooting='force-unrooted',
                                 taxon_namespace=taxa)

        tree.is_rooted = False
        tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
        tree.update_bipartitions()

        trees.append(tree)

    nt = len(trees)

    found = set([])
    sames = {}

    for i in range(0, nt - 1):
        ii = i + 1
        if ii not in found:
            found.add(ii)
            same = []
            for j in range(ii, nt):
                jj = j + 1

                ti = trees[i]
                tj = trees[j]

                [fp, fn] = false_positives_and_negatives(ti, tj)

                if fp + fn == 0:
                    same.append(str(jj))
                    found.add(jj)

            sames[str(ii)] = same

    for key, val in sames.items():
        sys.stdout.write("%s : %s\n" % (key, ' '.join(val)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--trees", type=str, nargs='+',
                        help="Input file containing trees",
                        required=True)

    main(parser.parse_args())
