import sys
import pandas as pd

def to_jskycalc(target_list):
    df = pd.read_csv(target_list)
    jskycalc = open("jskycalc_input.txt","w") 
    for i in range(df.shape[0]):
        ra = str(df['RA'].values[i])
        dec = str(df['Dec'].values[i])
        if "T mag" in list(df):
            ID = str(df['TIC ID'].values[i])
            mag = str(df['T mag'].values[i])
            name = "TESS"
        elif "Kep mag" in list(df):
            ID = str(df['EPIC ID'].values[i])
            mag = str(df['Kep mag'].values[i])
            name = "K2"        
        if len(list(ID)) > 7:
            jskycalc.write(ID+"\t"+ra+"\t"+dec+"\t2000.0\n")
        else:
            jskycalc.write(ID+"\t\t"+ra+"\t"+dec+"\t2000.0\n")
    jskycalc.close()

    return

if __name__ == "__main__":
    target_list = str(sys.argv[1])
    to_jskycalc(target_list)