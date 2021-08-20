import pandas
import numpy
import sys

df = pandas.read_csv("../../data/data-species-tree-error.csv")

mthds = ["parsimony-strict",
         "camin-sokal-strict",
         "dollo-strict",
         "astrid-bp",
         "astral-bp",
         "mdc-bp",
         "mdc-bp-ur"]

cols = ["MODL", "NRET", "MTHD", "THRC",
        "AVG_FNR", "AVG_FPR",
        "STD_FNR", "STD_FPR",
        "SE_FNR", "SE_FPR",
        "NREP"]
rows = []

for modl in ["5taxa", "6taxa", "26taxa", "Palaeognathae"]:
	wdf = df[df["MODL"] == modl]

	if (modl == "5taxa") or (modl == "6taxa"):
		nrets = [100000]
	else:
		nrets = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]

	for nret in nrets:
		xdf = wdf[wdf["NRET"] == nret]
		for mthd in mthds:
			ydf = xdf[xdf["MTHD"] == mthd]
			
			if ydf.shape[0] != 25:
				print(ydf)
				sys.exit("Missing replicates!\n")

			row = {}
			row["MODL"] = modl
			row["NRET"] = nret
			row["MTHD"] = mthd
			row["THRC"] = 0.0
			row["AVG_FNR"] = numpy.mean(ydf.FNR.values)
			row["AVG_FPR"] = numpy.mean(ydf.FPR.values)
			row["STD_FNR"] = numpy.std(ydf.FNR.values)
			row["STD_FPR"] = numpy.std(ydf.FPR.values)
			row["SE_FNR"] = numpy.std(ydf.FNR.values) / numpy.sqrt(25)
			row["SE_FPR"] = numpy.std(ydf.FPR.values) / numpy.sqrt(25)
			row["NREP"] = 25
			rows.append(row)

avg_df = pandas.DataFrame(rows, columns=cols)
avg_df.to_csv("mean-species-tree-error.csv", index=False)
