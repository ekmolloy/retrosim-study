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
        "AVG_FN", "AVG_FP", "AVG_FNR", "AVG_FPR",
        "STD_FN", "STD_FP", "STD_FNR", "STD_FPR",
        "SE_FN", "SE_FP", "SE_FNR", "SE_FPR",
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

			fnr = ydf["FN"] / ydf["NINT_TRUE"]
			fpr = ydf["FP"] / ydf["NINT_ESTI"]  # Correcting formula

			row = {}
			row["MODL"] = modl
			row["NRET"] = nret
			row["MTHD"] = mthd
			row["THRC"] = 0.0

			row["AVG_FN"] = numpy.mean(ydf["FN"])
			row["AVG_FP"] = numpy.mean(ydf["FP"])

			row["STD_FN"] = numpy.std(ydf["FN"])
			row["STD_FP"] = numpy.std(ydf["FP"])

			row["SE_FN"] = row["STD_FN"] / numpy.sqrt(25)
			row["SE_FP"] = row["STD_FP"] / numpy.sqrt(25)

			row["AVG_FNR"] = numpy.mean(fnr)
			row["AVG_FPR"] = numpy.mean(fpr)

			row["STD_FNR"] = numpy.std(fnr)
			row["STD_FPR"] = numpy.std(fpr)

			row["SE_FNR"] = row["STD_FNR"] / numpy.sqrt(25)
			row["SE_FPR"] = row["STD_FPR"] / numpy.sqrt(25)

			row["NREP"] = 25
			rows.append(row)

avg_df = pandas.DataFrame(rows, columns=cols)
avg_df.to_csv("mean-species-tree-error.csv", index=False)
