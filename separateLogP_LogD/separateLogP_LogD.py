import csv

def separateLogD_LogP(input_file,logP_file,logD_file):

    with open (input_file,'r') as f:
        next(f)
        data_lst = []
        for line in f:
            data_lst.append(line.rstrip('\n').split(','))

        RECORDID_logD= []
        RECORDID_logP= []

        """the input file has only the data with annotated pH, 19 collumns in the format 
           below, with *logd_pH corresponding to calculated logD to pH from 0 to 14"""
        for SMILES,RECORDID,pH,logd_74,*logd_pH in data_lst:
            appended_logP = False
            logP = max(float(logd_74), *[float(logd) for logd in logd_pH])

            if float(pH) == 7.4: 
                if (float(logP) - float(logd_74)) <= 0.3:
                    RECORDID_logP.append(RECORDID)
                else:
                    RECORDID_logD.append(RECORDID)          
            else:
                for pH_value, logd_value in enumerate(logd_pH):
                    if float(pH) == pH_value:
                        if (float(logP) - float(logd_value)) <= 0.3:
                            RECORDID_logP.append(RECORDID)
                        else:
                            RECORDID_logD.append(RECORDID)

    with open(logP_file, 'w', newline='') as fileP, open(logD_file, 'w', newline='') as fileD:
        writerP = csv.writer(fileP)
        writerP.writerow(["RECORDID"])
        for RID_P in RECORDID_logP:
            writerP.writerow([RID_P])

        writerD = csv.writer(fileD)
        writerD.writerow(["RECORDID"])
        for RID_D in RECORDID_logD:
            writerD.writerow([RID_D])


""" Although the data is not publicly available yet, the 'example_input.csv' 
has three lines of examples calculated by ChemAxon. """   

separateLogD_LogP("example_input.csv","logP.csv","logD.csv")