## imports
from EnzymePh import *
import pandas as pd
import argparse
import string
import numpy as np
import pickle
import joblib
import string
from sklearn.ensemble import RandomForestRegressor,ExtraTreesRegressor


## argumnet parser 
parser = argparse.ArgumentParser()
parser.add_argument('--input_csv', default='input.csv')
args = parser.parse_args()

## output file
output_file_path = 'predictions.csv'


## take input as sequence and convert it to descriptor dataframe
composition = PredictFrame(args.input_csv,descriptors_list= ["CC",])
Biopyth = PredictFrame(args.input_csv,descriptors_list= ["BioPythDescriptor",])
AminoAcid = PredictFrame(args.input_csv,descriptors_list = ["AAC"])

## data saved in variable
comp = composition.show()
Bio = Biopyth.show()
aa = AminoAcid.show()

## columns to be removed from the 
columns =[*[alpha for alpha in string.ascii_uppercase if alpha not in ["B","J","O","U","X","Z"]],]

## dropping colunms Biopyth
df = Bio.drop(columns,axis=1)

##conacting the dataframe
dataf = pd.concat([df,comp,aa],axis=1,join="outer")

## loading  model
with open("model.pkl","rb") as predictive_model:
    model = joblib.load(predictive_model)



## transforming the data to array for prediction and writing output file
y_predictions = []
for i in range(len(df)):
    s = model.predict(np.array(dataf.iloc[i]).reshape(1,-1))
    y_predictions.append(s[0])


# Save predictions to file
df_predictions = pd.DataFrame({'prediction': y_predictions})
df_predictions.to_csv(output_file_path, index=False)


    



