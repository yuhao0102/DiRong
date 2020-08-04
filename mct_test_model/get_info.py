import csv

result = [["decomp_type_id", "atm_point_num", "proc_num", "overlapping_rate", ]]
with open("record_time_atm", "r") as f:
	lines = f.readlines()
	i = 0
	while(i < len(lines)):
		line22 = []
		line1 = lines[i].strip().split(" ")
		for ii in line1:
			if (ii <> ''):
				line22.append(ii)
		line1 = lines[i+1].strip().split(" ")
                for ii in line1:
                        if (ii <> ''):
                                line22.append(ii)
		line1 = lines[i+2].strip().split(" ")
		for ii in line1:
                        if (ii <> ''):
                                line22.append(ii)
		i=i+3
		result.append(line22)
csv_writer = csv.writer(open("mct_atm.csv","w"))
csv_writer.writerows(result)
result = []
with open("record_time_ocn", "r") as f:
        lines = f.readlines()
        i = 0
        while(i < len(lines)):
                line22 = []
                line1 = lines[i].strip().split(" ")
                for ii in line1:
                        if (ii <> ''):
                                line22.append(ii)
                line1 = lines[i+1].strip().split(" ")
                for ii in line1:
                        if (ii <> ''):
                                line22.append(ii)
                line1 = lines[i+2].strip().split(" ")
                for ii in line1:
                        if (ii <> ''):
                                line22.append(ii)
                i=i+3
                result.append(line22)
csv_writer = csv.writer(open("mct_ocn.csv","w"))
csv_writer.writerows(result)

