import pandas as pd


coords_df = pd.read_csv('coords.csv', sep = ';')

print(coords_df['lon_esq'][0])