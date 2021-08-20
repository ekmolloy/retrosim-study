"""
Script to convert alignment in NEXUS format to PHYLIP format

Written by Erin K. Molloy (Github: ekmolloy)
"""
import argparse
import sys


def convert_nexus_to_phylip(infp, outfp, saln, ealn):
    ntaxa = -1
    nsites = -1

    # Check if NEXUS file
    line = infp.readline()
    if line[0] != '#':
        raise Exception("Input is not a NEXUS file!")

    # Process header
    for line in infp:
        temp = line.lower()

        word = "ntax="
        sind = temp.find(word)
        if sind > -1:
            eind = sind + len(word)
            ntaxa = int(temp[eind:].split()[0].replace(';', ''))

        word = "nchar="
        sind = temp.find(word)
        if sind > -1:
            eind = sind + len(word)
            nsites = int(temp[eind:].split()[0].replace(';', ''))

        word = "matrix"
        sind = temp.find(word)
        if sind > -1:
            if (ntaxa == -1) or (nsites == -1):
                raise Exception("Bad Input: Matrix block appears "
                                    "before definition of ntax or nchar!")
            break
        else:
            continue


    if saln is None:
        saln = 0

    if ealn is None:
        ealn = nsites
    else:
        if (saln > ealn) or (ealn > nsites):
            raise("Unable to index alignment properly!")
        nsites = ealn - saln

    outfp.write("%d %d\n" % (ntaxa, nsites))

    # Process matrix
    for line in infp:
        if line[:100].find(";") > -1:
            break

        temp = line.split()
        if len(temp) != 2:
            continue

        name = temp[0]
        data = temp[1].strip()
        data = data[saln:ealn]
        data = data.replace('0', 'A')
        data = data.replace('1', 'T')

        outfp.write("%s %s\n" % (name, data))


def main(args):
    with open(args.input, 'r') as infp, open(args.output, 'w') as outfp:
        convert_nexus_to_phylip(infp, outfp, args.saln, args.ealn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Alignment in NEXUS format",
                        required=True)

    parser.add_argument("-s", "--saln", type=int,
                        help="Start of alignment (inclusive; default 0)",
                        required=False)

    parser.add_argument("-e", "--ealn", type=int,
                        help="End of alignment (exclusive; default L)",
                        required=False)

    parser.add_argument("-o", "--output", type=str,
                        help="Output file in PHYLIP format",
                        required=True)

    main(parser.parse_args())