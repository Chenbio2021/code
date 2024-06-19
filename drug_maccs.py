# %%
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
import numpy as np
import pandas as pd

# %%
# 读取包含 name 和 SMILES 列的 CSV 文件
df = pd.read_csv('drug_smiles.csv')

# 创建一个空的 MACCS 指纹集合
maccs_df = pd.DataFrame()

# %%
# 计算每个分子的 MACCS 指纹
for index, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol is not None:
        maccs_key = MACCSkeys.GenMACCSKeys(mol)
        maccs_series = pd.Series(list(maccs_key.ToBitString()), name=row['Name'])
        maccs_df = maccs_df.append(maccs_series)


# %%
print(maccs_df)

# %%

# 将 DataFrame 保存为 CSV 文件
maccs_df.to_csv('drug_maccs.csv')


