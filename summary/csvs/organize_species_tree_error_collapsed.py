import pandas
import numpy
import sys

df = pandas.read_csv("../../data/data-species-tree-error-collapsed.csv")

mthds = ["astral-bp"]

cols = ["MODL", "NRET", "MTHD", "THRC",
        "AVG_FN", "AVG_FP", "AVG_FNR", "AVG_FPR",
        "STD_FN", "STD_FP", "STD_FNR", "STD_FPR",
        "SE_FN", "SE_FP", "SE_FNR", "SE_FPR",
        "NREP"]

rows = []

for modl in ["26taxa", "Palaeognathae"]:
	for mthd in mthds:
		for thrc in [0.5, 0.6, 0.7, 0.8, 0.9]:
			for nret in [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]:
				zdf = df[(df["NRET"] == nret) & 
				         (df["MTHD"] == mthd) &
				         (df["MODL"] == modl) &
				         (df["THRC"] == thrc)]				

				if zdf.shape[0] != 25:
					sys.exit("Missing replicates!\n")

				fnr = zdf["FN"] / zdf["NINT_TRUE"]
				fpr = zdf["FP"] / zdf["NINT_ESTI"]  # Correcting formula

				#x = numpy.mean(zdf["FPR"])
				#y = numpy.mean(fpr)
				#print("%s model, %d nret, %f - %f (old) vs. %f (new)" % (modl, nret, thrc, x, y))
				#print(zdf["NINT_TRUE"].values)
				#print(zdf["NINT_ESTI"].values)
				#for x,y in zip(zdf["FPR"].values, fpr.values):
				#	print("%f vs. %f" % (x,y))

				row = {}
				row["MODL"] = modl
				row["NRET"] = nret
				row["MTHD"] = mthd
				row["THRC"] = thrc

				row["AVG_FN"] = numpy.mean(zdf["FN"])
				row["AVG_FP"] = numpy.mean(zdf["FP"])

				row["STD_FN"] = numpy.std(zdf["FN"])
				row["STD_FP"] = numpy.std(zdf["FP"])

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
avg_df.to_csv("mean-species-tree-error-collapsed.csv", index=False)
