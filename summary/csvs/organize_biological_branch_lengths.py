import pandas
import sys

df = pandas.read_csv("../../data/data-branch-info.csv")

new_cols = ["BIPA", "BIPB", "BRLN_TENT", "BRLN_MAP_GT", "BRLN_MLE_GT", "BRLN_MLE_RI", "Q1", "Q2", "Q3", "EN", "PP"]

df1 = df[df["DATA"] == "Palaeognathae"]
df1 = df1[df1["NRET"] != 10]
df1 = df1[df1["NRET"] != 50]
df1 = df1[df1["NRET"] != 100]
df1 = df1[df1["NRET"] != 500]
df1 = df1[df1["NRET"] != 1000]
df1 = df1[df1["NRET"] != 5000]
df1 = df1[df1["NRET"] != 10000]
df1 = df1[df1["NRET"] != 50000]

df2 = df[df["DATA"] == "biological"]

new_rows = []
bipAs = list(set(list(df1["BIPA"].values) + list(df2["BIPA"].values)))

for bipA in bipAs:
	bip1 = set(bipA.split(','))
	bip2 = set(bipA.split(','))
	b1mb2 = bip1.difference(bip2)
	b2mb1 = bip2.difference(bip1)
	if len(b1mb2) != 0:
		sys.exit("Problem with taxon labeling: %s\n" % b1mb2)
	if len(b2mb1) != 0:
		sys.exit("Problem with taxon labeling: %s\n" % b2mb1)

	xdf1 = df1[df1["BIPA"] == bipA]
	xdf2 = df2[df2["BIPA"] == bipA]

	if xdf2.shape[0] == 0:
		sys.exit("Different bipartition sets!\n")

	row = {}
	row["BIPA"] = bipA

	row["BIPB"] = xdf1["BIPB"].values[0]
	row["BRLN_TENT"] = xdf1.BRLN.values[0]

	row["BRLN_MAP_GT"] = xdf2.BRLN_MAP_GT.values[0]
	row["BRLN_MLE_GT"] = xdf2.BRLN_MLE_GT.values[0]
	row["BRLN_MLE_RI"] = xdf2.BRLN_MLE_RI.values[0]
	row["Q1"] = xdf2.Q1.values[0]
	row["Q2"] = xdf2.Q2.values[0]
	row["Q3"] = xdf2.Q3.values[0]
	row["EN"] = xdf2.EN.values[0]
	row["PP"] = xdf2.PP.values[0]

	new_rows.append(row)

new_df = pandas.DataFrame(new_rows, columns=new_cols)
new_df.to_csv("biological-branch-lengths.csv", index=False)
