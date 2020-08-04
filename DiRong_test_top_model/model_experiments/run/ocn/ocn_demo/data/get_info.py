import csv

result = [["decomp_type_id", "atm_point_num", "ocn_grid_name", "proc_num", "overlapping_rate", "new algorithm sum", "old algorithm sum", "coupling time"]]
with open("record_time", "r") as f:
	lines = f.readlines()
	i = 0
	while(i < len(lines)):
		line1 = lines[i].strip().split(" ")
		line22 = []
		for ii in range(0,len(line1)):
			if(line1[ii] <> ''):
				line22.append(line1[ii])

		temp1 = float(lines[i+1].strip().split(" ")[5])
		temp2 = float(lines[i+3].strip().split(" ")[5])
		line22.append(temp1+temp2)		

                temp1 = float(lines[i+2].strip().split(" ")[5])
                temp2 = float(lines[i+4].strip().split(" ")[5])
		line22.append(temp1+temp2)

		temp1 = float(lines[i+5].strip().split(" ")[6])
		line22.append(temp1)

		i=i+6
		result.append(line22)
csv_writer = csv.writer(open("ocn_time_cc.csv","w"))
csv_writer.writerows(result)
