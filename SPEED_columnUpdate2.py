import pandas as pd
import time
import math

myinput = '2018_Speed_QCpassed.tsv'
print("loading csv in to dataframe ... ")
# load csv into pandas dataframe
df1 = pd.read_csv(myinput, header=0, delimiter='\t')
# count the number of times a gene is reported per patient and generate a df to hold the tally
print("counting variants per gene per patient ... ")
grouped = df1.groupby(['bridge_id', 'symbol'])
df2 = pd.DataFrame(grouped.size().reset_index(name = "Group_Count"))

# merge df1 and df2 to permit looping through and updating of model if variant have not passed QC
print("merging dataframes ... ")
df3 = df1.merge(df2, how='left')
#df4 = df3.head(5000)

# remove df1 and df2 to free up memory
del df1
del df2
# del df3

# loop through df3 and update post QC model if "Multi" model is recorded but only 1 variant passed QC
print("Updating model flag to reflect variants that have pass QC ...")

def test_get_value(df):
    """First attempt: To update models, runs really slowly on large data, on 90K dataframe this did not compleate after 12 hours"""
    for i in df.index:
        print(i)
        model = df.loc[i, 'model']
        if "MULT" in model and df.loc[i, 'Group_Count'] < 2:
            model = model.replace('MULT', 'SING')
            df.loc[i, 'model_postqc'] = model
        else:
            df.loc[i, 'model_postqc'] = df.loc[i,'model']

def two(df3):
    """Attempt two: To update models, Also runs slowly on large data, not much difference from attempt 1"""
    for index, row in df3.iterrows():
        print(index)
        if "MULT" in df3.loc[index, 'model'] and df3.loc[index, 'Group_Count'] < 2:
            flag = row['model']
            flag = flag.replace('MULT', 'SING')
            df3.loc[index, 'model_postqc'] = flag
        else:
            df3.loc[index, 'model_postqc'] = df3.loc[index, 'model']

def update_model (df):
    """take 3: Update model on whole dataframe - vectorized; so does not loop through. Much Much faster! ran 90K rows in 60s"""
    # df.loc[(df['A'] == 'blue') & (df['B'] == 'red') & (df['B'] == 'square'), 'label'] = 'M5'
    df['model_postqc'] = df['model']
    df.loc[(df['model'].str.contains("MULT")) & (df['Group_Count'] < 2), 'model_postqc'] = df.model.replace({'MULT':'SING'},regex=True)
    # print (df['model'], df['model_postqc'])
    # export processed df to csv
    print("exporting to csv ...")
    df.to_csv('SPEED_coding_snvs_processed5.csv', sep='\t' , index = False)

start = time.time()
# test_get_value(df4)
# two(df4)
update_model(df3)
end = time. time()
print(end - start)

exit()




