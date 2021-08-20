import pandas
import math
import numpy

modls = ["5-taxon-AZ-Sim",
		 "6-taxon-AZ-Sim",
		 "Palaeognathae-AZ-Sim",
		 "26-taxon-Sim"]

dfs = []
dfs.append(pandas.read_csv("../../csvs/5taxa-mean-simulated-branch-lengths.csv"))
dfs.append(pandas.read_csv("../../csvs/6taxa-mean-simulated-branch-lengths.csv"))
dfs.append(pandas.read_csv("../../csvs/Palaeognathae-mean-simulated-branch-lengths.csv"))
dfs.append(pandas.read_csv("../../csvs/26taxa-mean-simulated-branch-lengths.csv"))

with open("sim-brln-table.txt", 'w') as f:
	f.write("\\begin{table}\n")
	f.write("\\caption{Add caption...}\n")
	f.write("\\begin{tabular}{cccr}\n")
	f.write("\\toprule\n")
	f.write("True Branch Length & Estimated Branch Length & Percent Error & EN \\\\\n")
	f.write(" & MAP-GT / MLE-GT / MLE-RI & MAP-GT / MLE-GT / MLE-RI & \\\\\n")
	f.write("\\midrule\n")
	for i in range(len(modls)):
		f.write("\\multicolumn{4}{l}{\\em %s}\\\\\n" % modls[i])
		df = dfs[i]
		inds = numpy.argsort(df.BRLN_TRUE.values)
		for j in inds:
			xdf = df.iloc[j]
			brln_true = round(float(xdf.BRLN_TRUE), 4)

			brln_map_gt = round(round(xdf.AVG_BRLN_MAP_GT, 5), 4)
			brln_mle_gt = round(round(xdf.AVG_BRLN_MLE_GT, 5), 4)
			brln_mle_ri = round(round(xdf.AVG_BRLN_MLE_RI, 5), 4)

			nrepl_egtt_map_gt = xdf.NREPL_EGTT_MAP_GT
			nrepl_egtt_mle_gt = xdf.NREPL_EGTT_MLE_GT
			nrepl_egtt_mle_ri = xdf.NREPL_EGTT_MLE_RI

			perc_map_gt = round(round(xdf.AVG_PERR_MAP_GT, 1), 0)
			perc_mle_gt = round(round(xdf.AVG_PERR_MLE_GT, 1), 0)
			perc_mle_ri = round(round(xdf.AVG_PERR_MLE_RI, 1), 0)
			en = round(round(xdf.AVG_EN, 3), 2)

			if perc_mle_ri < perc_map_gt:
				f.write("%1.4f & %1.4f / %1.4f / %1.4f & %1.2f (%d) / %1.2f (%d) / {\\bf %1.2f} (%d) & %1.2f \\\\\n" \
					% (brln_true, brln_map_gt, brln_mle_gt, brln_mle_ri, perc_map_gt, nrepl_egtt_map_gt, perc_mle_gt, nrepl_egtt_mle_gt, perc_mle_ri, nrepl_egtt_mle_ri, en))
			else:
				f.write("%1.4f & %1.4f / %1.4f / %1.4f & %1.2f (%d) / %1.2f (%d) / %1.2f (%d) & %1.2f \\\\\n" \
					% (brln_true, brln_map_gt, brln_mle_gt, brln_mle_ri, perc_map_gt, nrepl_egtt_map_gt, perc_mle_gt, nrepl_egtt_mle_gt, perc_mle_ri, nrepl_egtt_mle_ri, en))

	f.write("\\bottomrule\n")
	f.write("\\end{tabular}\n")
	f.write("\\end{table}\n")
