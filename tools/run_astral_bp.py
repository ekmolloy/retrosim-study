"""
Script to run ASTRAL_BP

Written by Erin K. Molloy (Github: ekmolloy)
"""
import argparse
import numpy
import os
from scipy import special
import sys
import treeswift


class BinarySitePatternMatrix:
    def __init__(self, myfp, mykeep=None):
        self.ntaxa = -1
        self.nsites = -1
        self.keep = mykeep
        self.i2xmap = []
        self.x2imap = {}
        self.read_input(myfp)

    def read_input(self, fp):
        line = fp.readline()
        if line[0] == '#':
            self.read_nexus(fp)
        elif line[0] == '>':
            raise Exception("Bad Input: No support for FASTA files; "
                            "convert to PHYLIP format!")
        else:
            temp = line.split()
            if len(temp) != 2:
                raise Exception("Bad Input: Unable to read as PHYLIP file!")
            self.ntaxa = int(temp[0])
            self.nsites = int(temp[1])
            read_binary_matrix_block(fp)

    def read_nexus(self, fp):
        for line in fp:
            temp = line.lower()

            word = "ntax="
            sind = temp.find(word)
            if sind > -1:
                eind = sind + len(word)
                ntaxa = temp[eind:].split()[0].replace(';', '')
                self.ntaxa = int(ntaxa)

            word = "nchar="
            sind = temp.find(word)
            if sind > -1:
                eind = sind + len(word)
                nsites = temp[eind:].split()[0].replace(';', '')
                self.nsites = int(nsites)

            word = "matrix"
            sind = temp.find(word)
            if sind > -1:
                if (self.ntaxa == -1) or (self.nsites == -1):
                    raise Exception("Bad Input: Matrix block appears "
                                    "before definition of ntax or nchar!")
                else:
                    if self.keep is None:
                        self.keep = self.nsites
                    else:
                        if self.keep > self.nsites:
                            raise("Not enough sites in data set!")

                    self.read_binary_matrix_block(fp)

    def read_binary_matrix_block(self, fp):
        self.matrix = numpy.chararray((self.keep, self.ntaxa))

        j = 0
        name = ""

        for line in fp:
            if line[:100].find(";") > -1:
                break

            temp = line.split()
            if len(temp) != 2:
                continue

            self.i2xmap.append(temp[0])
            self.x2imap[temp[0]] = j

            for i, c in enumerate(temp[1]):
                if i == self.keep:
                    break
                self.matrix[i, j] = c

            j += 1

    def convert_site_to_newick(self, rowind):
        str0 = "("
        str1 = "("
        num0 = 0
        num1 = 0

        for j, label in enumerate(self.i2xmap):
            c = self.matrix[rowind, j]
            if c == b'0':
                if num0 == 0:
                    str0 += label
                else:
                    str0 += "," + label
                num0 += 1
            if c == b'1':
                if num1 == 0:
                    str1 += label
                else:
                    str1 += "," + label
                num1 += 1

        if (num0 > 1) and (num1 > 1):
            return str0 + "," + str1 + "))\n"
        else:
            return ""

    def write_newick(self, fp):
        for i in range(self.keep):
            fp.write(self.convert_site_to_newick(i))


class AstralBPtree:
    def __init__(self, myinput, myoutput, myargs, myspmat=None):
        # Process ASTRAL arguments
        ajar = myargs.astraljar
        aargs = []
        aargs.append("-c 0.5")
        aargs.append("-t 2")
        if myargs.exact:
            aargs.append("-x")
        if myargs.extratrees is not None:
            aargs.append(str("-e %s" % myargs.extratrees))
        if myargs.mintaxop is not None:
            aargs.append(str("-m %d" % myargs.mintaxop))
        if myargs.outgroup is not None:
            aargs.append(str("--outgroup %s" % myargs.outgroup))
        if myargs.savecompleted is not None:
            aargs.append(str("-k %s" % myargs.savecompleted))
        if myargs.tree is not None:
            aargs.append("-q %s" % myargs.tree)
        aargs.append(str("-i %s" % myinput))
        aargs.append(str("-o %s" % myoutput))

        # Run ASTRAL
        cmd = "java -jar " + ajar + " " + ' '.join(aargs)
        os.system(cmd)

        self.tree = treeswift.read_tree_newick(myoutput)
        os.remove(myoutput)

        # Store branch length information
        self.ntaxa = 0
        self.x2imap = {}
        self.i2xmap = []
        xind = 0
        for node in self.tree.traverse_postorder():
            if node.is_leaf():
                node.brinfo = None
                self.ntaxa += 1
                self.i2xmap.append(node.label)
                self.x2imap[node.label] = xind
                xind += 1
            elif node.label is not None:
                node.brinfo = self.parse_branch_label(node.label, prefix="AA")
                z1 = node.brinfo["AA_f1"]
                n = node.brinfo["AA_EN"]
                map_brln_gt = self.map_brln_gt(z1, n, 0.5)
                if round(node.edge_length, 6) != round(map_brln_gt, 6):
                    raise("Unable to verify ASTRAL's branch length!")
                node.brinfo["AA_GT_MAPBL"] = map_brln_gt
                node.brinfo["AA_GT_MLEBL"] = self.mle_brln_gt(z1, n)
                node.brinfo["AA_RI_MLEBL"] = self.mle_brln_ri(z1, n)
            else:
                node.brinfo = None

        self.found_short_quartets = False
        self.computed_freqs_for_short_quartets = False
        self.spmat = myspmat

    def parse_branch_label(self, label, prefix=None):
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
        for x in label[2:-2].split(';'):
            key, value = x.split('=')
            if prefix is not None:
                key = prefix + '_' + key
            data[key] = float(value)
        return data

    def mle_brln_gt(self, z, n):
        if n == 0:
            return numpy.nan

        t = float(z) / n

        if t <= 1.0 / 3:
            return 0.0
        if t == 1:
            return numpy.inf

        return -1.0 * numpy.log(1.5 * (1.0 - t))

    def map_brln_gt(self, z, n, c):
        x = z / (n + (2 * c))
        toreturn = -1.0 * numpy.log(1.5 * (1. - x))

        if toreturn < 0:
            return 0.0

        return toreturn

    def mle_brln_ri(self, z, n):
        if n == 0:
            return numpy.nan

        t = float(z) / n

        if t <= 1. / 3:
            return 0
        elif t == 1.:
            return numpy.inf

        y = ((2.0 / 3) / (1.0 - t)) - 1.0
        x = special.lambertw(y, k=0, tol=1e-08)

        if x < 0.0:
            return 0.0

        return x.real

    def localpp(self, z1, z2, z3, n, c):
        """
        Computes local posterior probability (local pp) defined in
        Equation (S4) of Sayyari and Mirarab (2016) as

            y = h(z1) / [h(z1) + 2^(z2-z1)h(z2) + 2^(z3-z1)h(z3)]

        where

            h(zi) = B(zi+1, n-zi+2c) * (1 - I_{1/3}(zi+1, n-zi+2c)).

        Recall that B is the beta function

            B(a, b) = (Gamma(a) * Gamma(b)) / Gamma(a+b)

        and I is the regularized incomplete beta function.

        Let ai = zi+1 and bi = n - zi + 2c.

        Then, we can rewrite y as

            1 - I_{1/3}(a1, b1)

        divided by

            [2^(zj-z1) * B(aj, bj)/B(a1, a2)] * (1 - I_{1/3}(aj, bj))

        summed for j = 1, 2, 3.

        In the ASTRAL code, the bracketted term is computed by first
        taking its log, presumably to avoid numerical issues, i.e.

            (zj-z1) * log(2) + log(B(aj, bj)) - log(B(a1, a2))

        becomes

            (zj-z1) * log(2) + [logG(aj) + logG(bj)] - [logG(a1) + logG(b1)]

        becase logG(aj+bj) = logG(a1+b1) = logG(n + 1 + 2c).
        """
        def ffunc(zi, n, c):
            toreturn = 1.0 - special.betainc(zi + 1, n - zi + 2 * c, 1.0 / 3)
            if toreturn <= 1e-15:
                return 0.0
            return toreturn

        def gfunc(zi, z1, n, c):
            ai = zi + 1
            bi = n - zi + (2 * c)
            a1 = z1 + 1
            b1 = n - z1 + (2 * c)

            x = numpy.log(2) * (zi - z1) \
                + (special.loggamma(ai) + special.loggamma(bi)) \
                - (special.loggamma(a1) + special.loggamma(b1))

            return numpy.exp(x)

        def fgfunc(zi, z1, n, c):
            toreturn = ffunc(zi, n, c) * gfunc(zi, z1, n, c)
            if (toreturn == numpy.nan) or (toreturn == numpy.inf):
                if (zi * 3 < n):
                    return 0.0
                elif (zi * 3 >= n):
                    if (z1 > zi):
                        return 0.0
                    elif (z1 < zi):
                        return numpy.inf
                    else:
                        raise("Something bad happened computing localpp for "
                              "zi=%d, z1=%d, n=%d, c=%f\n" % (zi, z1, n, c))
                else:
                    return numpy.inf
            return toreturn

        x1 = ffunc(z1, n, c)
        x2 = fgfunc(z2, z1, n, c)
        x3 = fgfunc(z3, z1, n, c)

        toreturn = x1 / (x1 + x2 + x3)

        if toreturn == numpy.inf:
            raise("Something bad happened computing localpp for "
                  "zi=%d, z1=%d, n=%d, c=%f\n" % (zi, z1, n, c))

        if toreturn == numpy.nan:
            if (z1 * 3 < n):
                return 0.0
            else:
                raise("Something bad happened computing localpp for "
                      "zi=%d, z1=%d, n=%d, c=%f\n" % (zi, z1, n, c))

        return toreturn

    def find_short_quartet_for_each_internal_branch(self):
        if self.found_short_quartets:
            return

        # Prepare intermediate data structures
        for node in self.tree.traverse_postorder():
            node.clade = numpy.zeros(self.ntaxa)
            node.dist2 = numpy.zeros(self.ntaxa)
            if node.is_leaf():
                node.edge_length = 0.0
            elif node.is_root():
                children = node.child_nodes()
                children[0].edge_length = 0.0
                children[1].edge_length = 0.0

        # Get distance to taxa in clade below edge
        for node in self.tree.traverse_postorder():
            if node.is_leaf():
                xind = self.x2imap[node.label]
                node.clade[xind] = 1.0
            else:
                children = node.child_nodes()
                node.clade = children[0].clade + children[1].clade
                if not node.is_root():
                    parent = node.get_parent()
                    for i in range(self.ntaxa):
                        if node.clade[i] != 0:
                            parent.dist2[i] = node.dist2[i] + node.edge_length

        # Get distance to taxa outside of clade below edge
        for node in self.tree.traverse_preorder():
            if node.brinfo is not None:
                parent = node.get_parent()
                for i in range(self.ntaxa):
                    if node.clade[i] == 0:
                        node.dist2[i] = parent.dist2[i] + node.edge_length

        # Find quartet with shortest terminal branches
        for node in self.tree.traverse_postorder():
            if node.brinfo is not None:
                inds = list(numpy.argsort(node.dist2))

                min0 = []
                min1 = []
                for i in inds:
                    if node.clade[i] == 0:
                        min0.append(i)
                    else:
                        min1.append(i)
                    if (len(min0) == 2) and (len(min1) == 2):
                        break

                node.sq = []
                node.sq.append(self.i2xmap[min0[0]])
                node.sq.append(self.i2xmap[min0[1]])
                node.sq.append(self.i2xmap[min1[0]])
                node.sq.append(self.i2xmap[min1[1]])

            node.clade = None
            node.dist2 = None

        self.i2xmap = None
        self.x2imap = None
        self.found_short_quartets = True

    def compute_info_for_short_quartets(self):
        if self.computed_freqs_for_short_quartets:
            return

        if self.spmat is None:
            raise("Need to provide site pattern matrix!")

        self.find_short_quartet_for_each_internal_branch()

        for node in self.tree.traverse_postorder():
            if node.brinfo is not None:
                # Frequentist estimate
                ABvsCD = 0
                ADvsBC = 0
                ACvsBD = 0

                iA = self.spmat.x2imap[node.sq[0]]
                iB = self.spmat.x2imap[node.sq[1]]
                iC = self.spmat.x2imap[node.sq[2]]
                iD = self.spmat.x2imap[node.sq[3]]

                for i in range(self.spmat.keep):
                    dA = self.spmat.matrix[i, iA]
                    dB = self.spmat.matrix[i, iB]
                    dC = self.spmat.matrix[i, iC]
                    dD = self.spmat.matrix[i, iD]

                    if (dA <= b'1') and (dB <= b'1') and \
                       (dC <= b'1') and (dD <= b'1'):
                        if (dA != dC):
                            if (dA == dB) and (dC == dD):
                                ABvsCD += 1.0
                            elif (dA == dD) and (dB == dC):
                                ADvsBC += 1.0
                        elif (dA != dB):
                            if (dA == dC) and (dB == dD):
                                ACvsBD += 1.0

                z1 = ABvsCD
                z2 = ADvsBC
                z3 = ACvsBD
                n = z1 + z2 + z3

                node.brinfo["SQ_f1"] = z1
                node.brinfo["SQ_f2"] = z2
                node.brinfo["SQ_f3"] = z3
                node.brinfo["SQ_N"] = n

                if n > 0:
                    node.brinfo["SQ_q1"] = z1 / n
                    node.brinfo["SQ_q2"] = z2 / n
                    node.brinfo["SQ_q3"] = z3 / n
                    
                else:
                    node.brinfo["SQ_q1"] = numpy.nan
                    node.brinfo["SQ_q2"] = numpy.nan
                    node.brinfo["SQ_q3"] = numpy.nan

                node.brinfo["SQ_pp1"] = self.localpp(z1, z2, z3, n, 0.5)
                node.brinfo["SQ_pp2"] = self.localpp(z2, z1, z3, n, 0.5)
                node.brinfo["SQ_pp3"] = self.localpp(z3, z1, z2, n, 0.5)
                node.brinfo["SQ_RI_MLEBL"] = self.mle_brln_ri(z1, n)

        self.computed_freqs_for_short_quartets = True

    def get_newick(self, blkey="AA_RI_MLEBL", bskey="AA_pp1"):
        for node in self.tree.traverse_preorder():
            if node.brinfo is not None:
                if bskey == "t2":
                    label = ""
                    for k, v in node.brinfo.items():
                        label = label + ";" + k + "=" + str("%1.6f" % v)
                    label = '\'' + label[1:] + '\''
                else:
                    label = node.brinfo[bskey]

                node.label = label
                node.edge_length = node.brinfo[blkey]

        return self.tree.newick()


def main(args):
    if args.newick and args.dosq:
        sys.exit("Unable to use options --newick and --dosq together!\n")

    # Process output arguments
    [outpath, outfile] = os.path.split(args.output_prefix)
    trfile1 = os.path.join(outpath, outfile + "-t2.tre")
    trfile2 = os.path.join(outpath, outfile + ".tre")

    # Process input arguments
    if args.newick:
        bpfile = args.input
        spmat = None
    else:
        # Transform 0/1 site patterns to newick strings
        [inpath, infile] = os.path.split(args.input)
        spfile = args.input
        bpfile = os.path.join(inpath, infile.rsplit('.', 1)[0] + '.bptrees')
        with open(spfile, 'r') as fp:
            spmat = BinarySitePatternMatrix(fp, mykeep=args.nsites)
        with open(bpfile, 'w') as fp:
            spmat.write_newick(fp)

    # Run ASTRAL
    astralbp = AstralBPtree(bpfile, trfile1, args, myspmat=spmat)

    if args.dosq:
        astralbp.compute_info_for_short_quartets()

    # Estimate branch lengths and branch support
    with open(trfile1, 'w') as fp:
        fp.write(astralbp.get_newick(bskey="t2") + '\n')

    with open(trfile2, 'w') as fp:
        fp.write(astralbp.get_newick() + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-j", "--astraljar", type=str,
                        help="ASTRAL .jar file, including path",
                        required=True)

    parser.add_argument("-x", "--exact", action="store_true",
                        help="Run ASTRAL in exact mode")

    parser.add_argument("-e", "--extratrees", type=str,
                        help="File with trees to expand ASTRAL's search space",
                        required=False)

    parser.add_argument("-k", "--savecompleted", type=str,
                        help="Save completed trees in ASTRAL's search space",
                        required=False)

    parser.add_argument("-m", "--mintaxop", type=int,
                        help="Minimum taxon occupancy to include RI",
                        required=False)

    parser.add_argument("-q", "--tree", type=str,
                        help="Use ASTRAL to score tree",
                        required=False)

    parser.add_argument("-g", "--outgroup", type=str,
                        help="Outgroup to root ASTRAL tree",
                        required=False)

    parser.add_argument("-dosq", "--dosq", action="store_true",
                        help="Use shortest quartet around branch")

    parser.add_argument("-n", "--nsites", type=int,
                        help="Run ASTRAL_BP on the first n sites",
                        required=False)

    parser.add_argument("-newick", "--newick", action="store_true",
                        help="Input is newick strings")

    parser.add_argument("-i", "--input", type=str,
                        help="Input file with RI (1/0) matrix (.nex or .phy)",
                        required=True)

    parser.add_argument("-o", "--output_prefix", type=str,
                        help="Output prefix",
                        required=True)

    main(parser.parse_args())
