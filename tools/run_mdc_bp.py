"""
Run MDC
"""
import argparse
import os
import sys


def main(args):
    tmpfile = args.output + ".tmp"
    nexfile = args.output + ".nex"
    trefile = args.output + ".tre"

    with open(args.input, 'r') as fi, \
         open(nexfile, 'w') as fo:

        genes = []
        fo.write("#NEXUS\n\n")
        fo.write("BEGIN NETWORKS;\n\n")
        for i, line in enumerate(fi):
            gene = str("gt%d" % (i + 1))
            genes.append(gene)
            xline = line.strip()
            if xline[-1] == ';':
                fo.write("Network %s = %s\n" % (gene, xline))
            else:
                fo.write("Network %s = %s;\n" % (gene, xline))
        fo.write("END;\n\nBEGIN PHYLONET;\n\nInfer_ST_MDC_UR (")
        fo.write(", ".join(genes))
        fo.write(") ")
        if args.unresolvedst:
            fo.write(" -ur ")
        if args.exact:
            fo.write(" -x ")
        fo.write(tmpfile + ";\n\nEND;")

    os.system("java -jar " + args.phylonet + " " + nexfile)
    os.system("tail -n1 " + tmpfile + "| awk \'{print $1}\' > " + trefile)
    os.system("rm " + nexfile)
    os.system("rm " + tmpfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--phylonet", type=str,
                        help="Path to PhyloNet",
                        required=True)

    parser.add_argument("-i", "--input", type=str,
                        help="Input file containing gene trees "
                             "(one newick string per line)",
                        required=True)

    parser.add_argument("-u", "--unresolvedst", action="store_true",
                        help="Allow non-binary species tree generation")

    parser.add_argument("-x", "--exact", action="store_true",
                        help="Run MDC in exact mode")

    parser.add_argument("-o", "--output", type=str,
                        help="Output file containing species tree "
                             "(newick string)",
                        required=True)

    main(parser.parse_args())
