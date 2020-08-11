import pandas as pd



if __name__ == "__main__":
   path = "/Users/Will/Desktop/TA sequences.xlsx"
   df = pd.read_excel(path)

   with open('ta_seqs.fasta', 'w') as file:
        for index, row in df.iterrows():
             file.write('>')
             file.write(row['enzyme_name'])
             file.write('\n')
             file.write(row['sequence'])
             file.write('\n')