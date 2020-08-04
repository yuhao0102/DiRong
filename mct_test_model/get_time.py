import csv

result = [["decomp", "atm_point_num", "ocn_grid_name", "proc_num", "total time", "router", "GlobalSegMap_init", "SparseMatrix(include read nc)", "GlobalSegMap_init(for interp_O)", "scatter in SparseMatrixPlus", "2function in SparseMatrixPlus", "SparseMatrixPlus"]]
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
                line22.append(float(ii))
        line1 = lines[i+2].strip().split(" ")
        for ii in line1:
            if (ii <> ''):
                line22.append(float(ii))
        line1 = lines[i+3].strip().split(" ")
	for ii in line1:
            if (ii <> ''):
                line22.append(float(ii))

        line1 = lines[i+4].strip().split(" ")
        line22.insert(2, line1[1])       
	i=i+5
        result.append(line22)
csv_writer = csv.writer(open("mct_atm_bscc.csv","w"))
csv_writer.writerows(result)


result = [["decomp", "atm_point_num", "ocn_grid_name", "proc_num", "total time", "router", "GlobalSegMap_init", "SparseMatrix(include read nc)", "GlobalSegMap_init", "function in SparseMatrixPlus", "SparseMatrixPlus", "scatter in SparseMatrixPlus"]]
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
                line22.append(float(ii))
        line1 = lines[i+2].strip().split(" ")
        for ii in line1:
            if (ii <> ''):
                line22.append(float(ii))
        line1 = lines[i+3].strip().split(" ")
        for ii in line1:
            if (ii <> ''):
                line22.append(float(ii))

	line1 = lines[i+4].strip().split(" ")
        line22.insert(2, line1[0])
        i=i+5    
        result.append(line22)
csv_writer = csv.writer(open("mct_ocn_bscc.csv","w"))
csv_writer.writerows(result)
