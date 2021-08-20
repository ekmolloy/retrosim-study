import pandas
import numpy
import sys

df = pandas.read_csv("../../data/data-species-tree-error-collapsed.csv")

mthds = ["astral-bp"]

cols = ["MODL", "NRET", "MTHD", "THRC",
        "AVG_FNR", "AVG_FPR",
        "STD_FNR", "STD_FPR",
        "SE_FNR", "SE_FPR",
        "NREP"]
rows = []

for modl in ["26taxa", "Palaeognathae"]:
	wdf = df[df["MODL"] == modl]
	for nret in [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]:
		xdf = wdf[wdf["NRET"] == nret]
		for mthd in mthds:
			ydf = xdf[xdf["MTHD"] == mthd]
			for thrc in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
				zdf = ydf[ydf["THRC"] == thrc]

				if zdf.shape[0] != 25:
					sys.exit("Missing replicates!\n")

				row = {}
				row["MODL"] = modl
				row["NRET"] = nret
				row["MTHD"] = mthd
				row["THRC"] = thrc
				row["AVG_FNR"] = numpy.mean(zdf.FNR.values)
				row["AVG_FPR"] = numpy.mean(zdf.FPR.values)
				row["STD_FNR"] = numpy.std(zdf.FNR.values)
				row["STD_FPR"] = numpy.std(zdf.FPR.values)
				row["SE_FNR"] = numpy.std(zdf.FNR.values) / numpy.sqrt(25)
				row["SE_FPR"] = numpy.std(zdf.FPR.values) / numpy.sqrt(25)
				row["NREP"] = 25
				rows.append(row)

avg_df = pandas.DataFrame(rows, columns=cols)
avg_df.to_csv("mean-species-tree-error-collapsed.csv", index=False)
