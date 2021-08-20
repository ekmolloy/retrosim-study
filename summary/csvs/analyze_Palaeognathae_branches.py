import pandas
import numpy
import sys

#dofp = 0
dofp = 1

nret = 5000

df = pandas.read_csv("../../data/data-branch-info.csv")
df = df[df["DATA"] == "Palaeognathae"]

# Find clades to study
tdf = df[df["MTHD"] == "true"]  # true
edf = df[df["NRET"] == nret]    # estimated

true_bipAs = set(tdf.BIPA.values)
esti_bipAs = set(edf.BIPA.values)

if dofp:
    bipAs = list(esti_bipAs.difference(true_bipAs))
    outf = "Palaeognathae-5000-false-positive-branches.csv"
else:
    bipAs = list(esti_bipAs.intersection(true_bipAs))
    outf = "Palaeognathae-5000-true-positives-branches.csv"

# Sort clades
clades = []
clade_sizes = []
for bipA in bipAs:
    bipB = edf[edf["BIPA"] == bipA].BIPB.values[0]

    if bipA.find("galGal") == -1:
        clade = bipA
    else:
        clade = bipB

    clades.append(clade)
    clade_sizes.append(len(clade.split(",")))

inds = numpy.argsort(clade_sizes)
inds = inds[::-1]

bipAs = [bipAs[i] for i in inds]
clades = [clades[i] for i in inds]

# Check if clades are in estimated species tree for each replicate data set
new_cols = ["CLADE"]
for repl in range(1, 26):
    new_cols.append(str("REPL%d" % repl))
new_rows = []

for clade, bipA in zip(clades, bipAs):
    xdf = edf[edf["BIPA"] == bipA]

    new_row = {}
    new_row["CLADE"] = clade

    for repl in range(1, 26):
        ydf = xdf[xdf["REPL"] == repl]

        key = str("REPL%d" % repl)

        if ydf.shape[0] > 1:
            print(ydf)
            sys.exit("Found too many rows")

        if ydf.shape[0] == 0:
            new_row[key] = "NA"
        else:
            localpp = ydf.PP.values[0]
            brlen = ydf.BRLN_MLE_RI.values[0]
            en = ydf.EN.values[0]

            new_row[key] = str("%1.4f" % (localpp))

            #str("%1.4f/%1.4f" % (brlen, localpp))

    new_rows.append(new_row)

new_df = pandas.DataFrame(new_rows, columns=new_cols)
new_df.to_csv(outf, index=False)
