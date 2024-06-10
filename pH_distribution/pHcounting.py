import csv

def createDistribution(input_file):

    pH_count = {f'pH{i}': 0 for i in range(0, 14, 1)}
    pH_count["pH7.4"] = 0

    with open ('measured_pH.csv','r', encoding='utf-8') as f:
        next(f)
        data = []
        for line in f:
            data.append(line.rstrip('\n').split(','))

    for SMILES,ID,pH in data:
        pH_ = float(pH)
        if pH_ == 7.4:
            pH_count["pH7.4"] +=1
        else:
            for i in range(0, 14, 1):
                if i <= pH_ < i + 1:
                    pH_count[f'pH{i}'] += 1
                    break

    print(pH_count)

createDistribution("measured_pH.csv")

