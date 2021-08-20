"""
Grab branch info

Written by Erin K. Molloy (github: ekmolloy)
"""
import argparse
import dendropy
import sys


def parse_branch_label(label):
    data = {}

    try:
        xs = label.split(';')
    except TypeError:
        return data

    for x in xs:
        try:
            key, value = x.split('=')
            data[key] = float(value)
        except ValueError:
            return data

    return data


def main(args):
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=args.input,
                             schema='newick',
                             rooting='force-unrooted',
                             taxon_namespace=taxa)

    taxa = set([x.label for x in taxa])

    for node in tree.preorder_node_iter():
        bipA = set([x.taxon.label for x in node.leaf_nodes()])
        bipB = taxa.difference(bipA)

        if (len(bipA) > 1) and (len(bipB) > 1):
            bipA = ','.join(sorted(list(bipA)))
            bipB = ','.join(sorted(list(bipB)))
            tmp = sorted([bipA, bipB])

            if node.edge.length is not None:
                brln = str("%f" % float(node.edge.length))
            else:
                brln = "NA"

            if node.label is not None:
                label = parse_branch_label(node.label)

                if len(label) > 0:
                    f1 = str("%f" % label["AA_q1"])
                    f2 = str("%f" % label["AA_q2"])
                    f3 = str("%f" % label["AA_q3"])
                    en = str("%f" % label["AA_EN"])
                    pp = str("%f" % label["AA_pp1"])
                    map_gt = str("%f" % label["AA_GT_MAPBL"])
                    mle_gt = str("%f" % label["AA_GT_MLEBL"])
                    mle_ri = str("%f" % label["AA_RI_MLEBL"])
                    info = ','.join([map_gt, mle_gt, mle_ri, f1, f2, f3, en, pp])
                else:
                    info = "NA,NA,NA,NA,NA,NA,NA," + label

            else:
                info = "NA,NA,NA,NA,NA,NA,NA,NA"

            sys.stdout.write("%s,\"%s\",\"%s\",%s,%s\n" %
                             (args.prefix, tmp[0], tmp[1], brln, info))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input newick string",
                        required=True)

    parser.add_argument("-p", "--prefix", type=str,
                        help="Prefix",
                        required=True)

    main(parser.parse_args())
