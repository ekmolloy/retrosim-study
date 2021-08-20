import pandas
import numpy

modls = ["4 taxa anomaly zone",
		 "5 taxa anomaly zone",
		 "Palaeognathae",
		 "26 taxa"]


df = pandas.read_csv("../../csvs/biological-branch-lengths.csv")
inds = numpy.argsort(df.BRLN_TENT.values)

with open("bio-brln-table.txt", 'w') as f:
	f.write("\\begin{landscape}\n")
	f.write("\\begin{table}\n")
	f.write("\\caption{Add caption...}\n")
	f.write("\\begin{tabular}{l c cc r}\n")
	f.write("\\toprule\n")
	f.write("Clade & ASTRAL TENT Analysis & \\multicolumn{3}{c}{ASTRAL\\_BP Analysis} \\\\\n")
	f.write("& Branch Length  & Branch Length & Quartet Support & EN \\\\\n")
	f.write("&  MAP-GT & MAP-GT / MLE-GT / MLE-RI & & \\\\\n")
	f.write("\\midrule\n")
	for j in inds:
		xdf = df.iloc[j]
		bipA = xdf.BIPA
		bipB = xdf.BIPB
		brln = round(round(xdf.BRLN_TENT, 5), 4)
		map_gt = round(round(xdf.BRLN_MAP_GT, 5), 4)
		mle_gt = round(round(xdf.BRLN_MLE_GT, 5), 4)
		mle_ri = round(round(xdf.BRLN_MLE_RI, 5), 4)
		q1 = round(round(xdf.Q1, 5), 4)
		q2 = round(round(xdf.Q2, 5), 4)
		q3 = round(round(xdf.Q3, 5), 4)
		en = round(round(xdf.EN, 3), 2)
		f.write("%s vs. %s & %1.4f & %1.4f / %1.4f / %1.4f & %1.4f / %1.4f / %1.4f & %1.2f \\\\\n" \
					% (bipA, bipB, brln, map_gt, mle_gt, mle_ri, q1, q2, q3, en))
	f.write("\\bottomrule\n")
	f.write("\\end{tabular}\n")
	f.write("\\end{table}\n")
	f.write("\\end{landscape}\n")




