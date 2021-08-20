"""
Run MDC
"""
import argparse
import os
import sys


def main(args):
    tmp = args.input
    tmp = tmp.replace(".nex", "").rsplit("/", 1)
    if len(tmp) == 1:
        name = tmp[0]
    else:
        name = tmp[-1]

    if args.outpath is None:
        outpath = ""
    else:
        outpath = args.outpath + '/'

    ogs = " ".join(args.outgroups)

    if args.parsimony is None:
        method = "parsimony"
    else:
        method = args.parsimony

    if args.keep < 100000:
        name = name + "-" + str(args.keep)
        nexfile = outpath + method + "-" + name + ".nex"
        logfile = outpath + method + "-" + name + ".log"
        output1 = outpath + method + "-strict-" + name + ".tre"
        output2 = outpath + method + "-all-" + name + ".trees"
    else:
        nexfile = outpath + method + "-" + name + ".nex"
        logfile = outpath + method + "-" + name + ".log"
        output1 = outpath + method + "-strict-" + name + ".tre"
        output2 = outpath + method + "-all-" + name + ".trees"

    with open(nexfile, 'w') as fo:
        fo.write("#NEXUS\n")
        fo.write("BEGIN PAUP;\n")
        fo.write("set autoclose=yes warntree=no warnreset=no;\n")
        fo.write("execute " + args.input + ";\n")
        fo.write("outgroup " + ogs + ";\n")
        if args.parsimony == "dollo":
            fo.write("ctype dollo:1-100000;")
        elif args.parsimony == "camin-sokal":
            fo.write("ctype irrev:1-100000;[irrev=Camin-Sokal]\n")
        if args.keep < 100000:
            fo.write("exclude " + str(args.keep + 1) + "-100000;\n")
        if args.bband:
            fo.write("bandb;\n")
        else:
            fo.write("hsearch addSeq=random nreps=100 swap=TBR;\n")
        fo.write("contree all/strict=yes treefile=" + output1)
        fo.write(" format=newick;\n")
        fo.write("savetrees File=" + output2)
        fo.write(" root=yes trees=all format=newick;\n")
        fo.write("END;\n")

    os.system(args.paup + " -n " + nexfile + " &> " + logfile)
    os.system("rm " + nexfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--paup", type=str,
                        help="Path to PAUP",
                        required=True)

    parser.add_argument("-p", "--parsimony", type=str,
                        help="Parsimony method [None, dollo, camin-sokal]")

    parser.add_argument("-b", "--bband", action="store_true",
                        help="Do branch and bound")

    parser.add_argument("-k", "--keep", type=int,
                        help="Number of retroelements "
                             "to keep e.g. 10, 50, 100, etc.",
                        required=True)

    parser.add_argument("-i", "--input", type=str,
                        help="Input nexus file",
                        required=True)

    parser.add_argument("-g", "--outgroups", type=str, nargs='+',
                        help="Outgroups e.g. "
                             "\'TaxonU TaxonV TaxonW TaxonX TaxonY TaxonZ\'",
                        required=True)

    parser.add_argument("-o", "--outpath", type=str,
                        help="Output path", required=False)

    main(parser.parse_args())
