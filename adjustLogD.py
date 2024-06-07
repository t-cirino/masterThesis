import csv

def adjustLogD(input_file, output_file):

    with open (input_file,'r') as infile, open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        next(infile)
        data_lst = []
        for line in infile:
            data_lst.append(line.rstrip('\n').split(','))

        writer_logP = csv.writer(outfile)
        writer_logP.writerow(["SMILES","RECORDID","pH","logPow"])

        # the .csv file must have 20 collumns, including logD7.4 and logD0 to logD14
        for SMILES,RECORDID,exp_logD,pH,logD_74,*logD_pH in data_lst:

            max_logD = max(float(logD_74), *[float(logd) for logd in logD_pH])
            maxLogD_pH = [7.4] if float(logD_74) == max_logD else []
            maxLogD_pH.extend([float(pH_value) for pH_value, logd_value in enumerate(logD_pH) if float(logd_value) == max_logD])
        
            if float(pH) == 7.4:
                adjusted_logP = float(exp_logD) + float(max_logD) - float(logD_74)
                adjusted_logP = round(adjusted_logP, 2)
                pH_theo = round(sum(maxLogD_pH)/len(maxLogD_pH),1)
                writer_logP.writerow([SMILES,RECORDID,pH_theo,adjusted_logP])

            else:
                for pH_value, logD_value in enumerate(logD_pH):
                    if float(pH) == pH_value:
                        adjusted_logP = float(exp_logD) + float(max_logD) - float(logD_value)
                        adjusted_logP = round(adjusted_logP, 2)
                        pH_theo = round(sum(maxLogD_pH)/len(maxLogD_pH),1)
                        writer_logP.writerow([SMILES,RECORDID,pH_theo,adjusted_logP])

adjustLogD("example_adjustLogD.csv", "example_output.csv")