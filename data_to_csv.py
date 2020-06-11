import signac
import pandas as pd

name = input('Provide file name (<file name>_data_raw.csv): ')

Data = signac.get_project()
df = Data.to_dataframe()

df.to_csv('{}_data_raw.csv'.format(name))
