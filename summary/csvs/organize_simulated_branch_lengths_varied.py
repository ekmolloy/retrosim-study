import pandas
import numpy
import sys

df = pandas.read_csv("../../data/data-branch-info.csv")

new_cols = ["NRET", "BIPA", "BIPB", "REPL",
			"BRLN_TRUE",
			"BRLN_MAP_GT", "EGTT_MAP_GT", "AERR_MAP_GT", "PERR_MAP_GT",
			"BRLN_MLE_GT", "EGTT_MLE_GT", "AERR_MLE_GT", "PERR_MLE_GT",
			"BRLN_MLE_RI", "EGTT_MLE_RI", "AERR_MLE_RI", "PERR_MLE_RI",
			"Q1", "Q2", "Q3", "EN", "PP"]

avg_cols = ["NRET", "NREP", "BIPA", "BIPB",
			"BRLN_TRUE",
			"AVG_BRLN_MAP_GT", "STD_BRLN_MAP_GT", "NREPL_EGTT_MAP_GT",
			"AVG_AERR_MAP_GT", "STD_AERR_MAP_GT",
			"AVG_PERR_MAP_GT", "STD_PERR_MAP_GT",
			"AVG_BRLN_MLE_GT", "STD_BRLN_MLE_GT", "NREPL_EGTT_MLE_GT",
			"AVG_AERR_MLE_GT", "STD_AERR_MLE_GT",
			"AVG_PERR_MLE_GT", "STD_PERR_MLE_GT",
			"AVG_BRLN_MLE_RI", "STD_BRLN_MLE_RI", "NREPL_EGTT_MLE_RI",
			"AVG_AERR_MLE_RI", "STD_AERR_MLE_RI",
			"AVG_PERR_MLE_RI", "STD_PERR_MLE_RI",
			"AVG_Q1", "AVG_Q2", "AVG_Q3", "AVG_EN", "AVG_PP",
			"STD_Q1", "STD_Q2", "STD_Q3", "STD_EN", "STD_PP"]

for modl in ["26taxa", "Palaeognathae"]:
	wdf = df[df["DATA"] == modl]

	tdf = wdf[wdf["MTHD"] == "true"]
	bipAs = list(set(tdf["BIPA"].values))

	new_rows = []
	avg_rows = []

	for bipA in bipAs:
		bipB = tdf[tdf["BIPA"] == bipA].BIPB.values[0]
		truebrln = tdf[tdf["BIPA"] == bipA].BRLN.values[0]
		xdf = wdf[wdf["BIPA"] == bipA]

		for nret in [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]:
			ydf = xdf[xdf["NRET"] == nret]  # Grabs ASTRAL-BP for NRET retroelement insertions
			nrep = ydf.shape[0]

			if nrep == 0:
				for repl in range(1, 26):
					row = {}
					row["NRET"] = nret
					row["BIPA"] = bipA
					row["BIPB"] = bipB
					row["REPL"] = repl
					row["BRLN_TRUE"] = truebrln
					new_rows.append(row)

				row = {}
				row["NRET"] = nret
				row["NREP"] = nrep
				row["BIPA"] = bipA
				row["BIPB"] = bipB
				row["BRLN_TRUE"] = truebrln
				avg_rows.append(row)
			else:
				brln_map_gt = []
				brln_mle_gt = []
				brln_mle_ri = []
				aerr_map_gt = []
				aerr_mle_gt = []
				aerr_mle_ri = []
				perr_map_gt = []
				perr_mle_gt = []
				perr_mle_ri = []
				nrepl_egtt_map_gt = 0
				nrepl_egtt_mle_gt = 0
				nrepl_egtt_mle_ri = 0
				q1 = []
				q2 = []
				q3 = []
				en = []
				pp = []

				for repl in range(1, 26):
					zdf = ydf[ydf["REPL"] == repl]

					if zdf.shape[0] > 1:
						print(zdf)
						sys.exit("Found too many rows")
					elif zdf.shape[0] == 0:
						row = {}
						row["NRET"] = nret
						row["BIPA"] = bipA
						row["BIPB"] = bipB
						row["REPL"] = repl
						row["BRLN_TRUE"] = truebrln
						new_rows.append(row)
					else:
						#if zdf.BRLN_MAP_GT.values[0] == numpy.inf:
						#	sys.stdout.write("WARNING: branch for %s replicate %d -- branch has infinite length!\n" % (modl, repl))
						#	sys.stdout.write("\t" + bipA + " vs. " + bipB + "\n")

						row = {}
						row["NRET"] = nret
						row["BIPA"] = bipA
						row["BIPB"] = bipB
						row["REPL"] = repl
						row["BRLN_TRUE"] = truebrln
						row["BRLN_MAP_GT"] = zdf.BRLN_MAP_GT.values[0]
						row["BRLN_MLE_GT"] = zdf.BRLN_MLE_GT.values[0]
						row["BRLN_MLE_RI"] = zdf.BRLN_MLE_RI.values[0]
						row["AERR_MAP_GT"] = abs(row["BRLN_TRUE"] - row["BRLN_MAP_GT"])
						row["AERR_MLE_GT"] = abs(row["BRLN_TRUE"] - row["BRLN_MLE_GT"])
						row["AERR_MLE_RI"] = abs(row["BRLN_TRUE"] - row["BRLN_MLE_RI"])
						row["PERR_MAP_GT"] = (row["AERR_MAP_GT"] / row["BRLN_TRUE"]) * 100.0
						row["PERR_MLE_GT"] = (row["AERR_MLE_GT"] / row["BRLN_TRUE"]) * 100.0
						row["PERR_MLE_RI"] = (row["AERR_MLE_RI"] / row["BRLN_TRUE"]) * 100.0

						# Check if estimated value is greater than true value!
						if row["BRLN_MAP_GT"] > row["BRLN_TRUE"]:
							row["EGTT_MAP_GT"] = 1
						else:
							row["EGTT_MAP_GT"] = 0

						if row["BRLN_MLE_GT"] > row["BRLN_TRUE"]:
							row["EGTT_MLE_GT"] = 1
						else:
							row["EGTT_MLE_GT"] = 0

						if row["BRLN_MLE_RI"] > row["BRLN_TRUE"]:
							row["EGTT_MLE_RI"] = 1
						else:
							row["EGTT_MLE_RI"] = 0

						row["Q1"] = zdf.Q1.values[0]
						row["Q2"] = zdf.Q2.values[0]
						row["Q3"] = zdf.Q3.values[0]
						row["EN"] = zdf.EN.values[0]
						row["PP"] = zdf.PP.values[0]
						new_rows.append(row)

						brln_map_gt.append(row["BRLN_MAP_GT"])
						brln_mle_gt.append(row["BRLN_MLE_GT"])
						brln_mle_ri.append(row["BRLN_MLE_RI"])
						aerr_map_gt.append(row["AERR_MAP_GT"])
						aerr_mle_gt.append(row["AERR_MLE_GT"])
						aerr_mle_ri.append(row["AERR_MLE_RI"])
						perr_map_gt.append(row["PERR_MAP_GT"])
						perr_mle_gt.append(row["PERR_MLE_GT"])
						perr_mle_ri.append(row["PERR_MLE_RI"])
						nrepl_egtt_map_gt += row["EGTT_MAP_GT"]
						nrepl_egtt_mle_gt += row["EGTT_MLE_GT"]
						nrepl_egtt_mle_ri += row["EGTT_MLE_RI"]
						q1.append(row["Q1"])
						q2.append(row["Q2"])
						q3.append(row["Q3"])
						en.append(row["EN"])
						pp.append(row["PP"])

				row = {}
				row["NRET"] = nret
				row["NREP"] = nrep
				row["BIPA"] = bipA
				row["BIPB"] = bipB
				row["BRLN_TRUE"] = truebrln

				inf_gt = ~numpy.isfinite(brln_mle_gt)
				inf_ri = ~numpy.isfinite(brln_mle_ri)
				if not numpy.array_equal(inf_gt, inf_ri):
					sys.exit("Issue pulling out infs!")
				keep = set(range(0, nrep))
				keep = numpy.array(list(keep.difference(set(numpy.where(inf_gt)[0]))))

				if len(keep) > 0:
					brln_map_gt = numpy.array(brln_map_gt)
					brln_mle_gt = numpy.array(brln_mle_gt)
					brln_mle_ri = numpy.array(brln_mle_ri)
					aerr_map_gt = numpy.array(aerr_map_gt)
					aerr_mle_gt = numpy.array(aerr_mle_gt)
					aerr_mle_ri = numpy.array(aerr_mle_ri)
					perr_map_gt = numpy.array(perr_map_gt)
					perr_mle_gt = numpy.array(perr_mle_gt)
					perr_mle_ri = numpy.array(perr_mle_ri)
					row["AVG_BRLN_MAP_GT"] = numpy.mean(brln_map_gt[keep])
					row["AVG_BRLN_MLE_GT"] = numpy.mean(brln_mle_gt[keep])
					row["AVG_BRLN_MLE_RI"] = numpy.mean(brln_mle_ri[keep])
					row["AVG_AERR_MAP_GT"] = numpy.mean(aerr_map_gt[keep])
					row["AVG_AERR_MLE_GT"] = numpy.mean(aerr_mle_gt[keep])
					row["AVG_AERR_MLE_RI"] = numpy.mean(aerr_mle_ri[keep])
					row["AVG_PERR_MAP_GT"] = numpy.mean(perr_map_gt[keep])
					row["AVG_PERR_MLE_GT"] = numpy.mean(perr_mle_gt[keep])
					row["AVG_PERR_MLE_RI"] = numpy.mean(perr_mle_ri[keep])
					row["STD_BRLN_MAP_GT"] = numpy.std(brln_map_gt[keep])
					row["STD_BRLN_MLE_GT"] = numpy.std(brln_mle_gt[keep])
					row["STD_BRLN_MLE_RI"] = numpy.std(brln_mle_ri[keep])
					row["STD_AERR_MAP_GT"] = numpy.std(aerr_map_gt[keep])
					row["STD_AERR_MLE_GT"] = numpy.std(aerr_mle_gt[keep])
					row["STD_AERR_MLE_RI"] = numpy.std(aerr_mle_ri[keep])
					row["STD_PERR_MAP_GT"] = numpy.std(perr_map_gt[keep])
					row["STD_PERR_MLE_GT"] = numpy.std(perr_mle_gt[keep])
					row["STD_PERR_MLE_RI"] = numpy.std(perr_mle_ri[keep])

				row["NREPL_EGTT_MAP_GT"] = nrepl_egtt_map_gt
				row["NREPL_EGTT_MLE_GT"] = nrepl_egtt_mle_gt
				row["NREPL_EGTT_MLE_RI"] = nrepl_egtt_mle_ri
				row["AVG_Q1"] = numpy.mean(q1)
				row["AVG_Q2"] = numpy.mean(q2)
				row["AVG_Q3"] = numpy.mean(q3)
				row["AVG_EN"] = numpy.mean(en)
				row["AVG_PP"] = numpy.mean(pp)
				row["STD_Q1"] = numpy.std(q1) 
				row["STD_Q2"] = numpy.std(q2)
				row["STD_Q3"] = numpy.std(q3)
				row["STD_EN"] = numpy.std(en)
				row["STD_PP"] = numpy.std(pp)
				avg_rows.append(row)

	new_df = pandas.DataFrame(new_rows, columns=new_cols)
	new_df.to_csv(modl + "-simulated-branch-lengths-varied.csv", index=False)

	avg_df = pandas.DataFrame(avg_rows, columns=avg_cols)
	avg_df.to_csv(modl + "-mean-simulated-branch-lengths-varied.csv", index=False)
