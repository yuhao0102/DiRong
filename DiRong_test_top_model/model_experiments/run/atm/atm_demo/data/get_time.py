import csv

result = [["decomp_type_id", "atm_point_num", "ocn_grid_name", "proc_num", "overlapping_rate", "new algorithm", "new algorithm", "new algorithm sum", "old algorithm", "old algorithm", "old algorithm sum", "coupling time"]]
with open("record_time", "r") as f:
	lines = f.readlines()
	i = 0
	while(i < len(lines)):
		line22 = []
		for t in range(0, 10):
			if(i >= len(lines)):
				break	
			line1 = lines[i].strip().split(" ")
			for ii in range(0,len(line1)):
				if(line1[ii] <> ''):
					line22.append(line1[ii])
			i = i + 1
		i += 1
		result.append(line22)
csv_writer = csv.writer(open("atm_time_for_cc_internal.csv","w"))
csv_writer.writerows(result)
