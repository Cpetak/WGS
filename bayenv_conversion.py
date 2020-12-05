import pandas as pd
import os
import numpy as np

files = ['/users/c/p/cpetak/WGS/angsd_new_noout/LOM_angsd_polysites.mafs', '/users/c/p/cpetak/WGS/angsd_new_noout/TER_angsd_polysites.mafs', '/users/c/p/cpetak/WGS/angsd_new_noout/BOD_angsd_polysites.mafs', '/users/c/p/cpetak/WGS/angsd_new_noout/KIB_angsd_polysites.mafs', '/users/c/p/cpetak/WGS/angsd_new_noout/CAP_angsd_polysites.mafs', '/users/c/p/cpetak/WGS/angsd_new_noout/FOG_angsd_polysites.mafs']
pops = ["LOM","TER","BOD","KIB","CAP","FOG"]
data = []

for i in range(len(pops)):
    frame = pd.read_csv(files[i], sep="\t")
    frame['filename'] = pops[i]
    data.append(frame)

def init_proc_file(data):
	print(len(data))
	if data['filename'].iloc[1] == "CAP":
		data=data[data["nInd"]==19]
	elif data['filename'].iloc[1] == "FOG":
                data=data[data["nInd"]==18]
	else:
		data=data[data["nInd"]==20]
	print(len(data))
	data=data.iloc[::100, :]
	data['chromo'] = data['chromo'].apply(str)
	data['position'] = data['position'].apply(str)
	data['pos'] = data[['chromo', 'position']].apply(lambda x: ''.join(x), axis=1)
	df = data[["pos","knownEM", "filename", "nInd"]]
	df["al_count"] = df["knownEM"] * df["nInd"]
	df["al_count"] = df.al_count.round().astype(int)
	df.drop(df[df['al_count'] < 1].index, inplace = True)
	df["al_count2"] = df["nInd"] - df["al_count"]
	df = df.reset_index()
	df = df[["pos","filename", "al_count", "al_count2"]]
	print(len(df))
	pop = df['filename'].iloc[1]
	df = df.rename(columns={"al_count": "al_count"+pop, "al_count2": "al_count2"+pop})
	df = df.drop(['filename'], axis=1)
	return df

lom_df = init_proc_file(data[0])
ter_df = init_proc_file(data[1])
bod_df = init_proc_file(data[2])
kib_df = init_proc_file(data[3])
cap_df = init_proc_file(data[4])
fog_df = init_proc_file(data[5])

merged = pd.merge(lom_df, ter_df, on="pos", how="outer")
merged = pd.merge(merged, bod_df, on="pos", how="outer")
merged = pd.merge(merged, kib_df, on="pos", how="outer")
merged = pd.merge(merged, cap_df, on="pos", how="outer")
merged = pd.merge(merged, fog_df, on="pos", how="outer")

merged = merged.fillna(0)
merged['al_count2LOM'] = np.where(merged['al_countLOM'] == 0, 20, merged['al_count2LOM'])
merged['al_count2TER'] = np.where(merged['al_countTER'] == 0, 20, merged['al_count2TER'])
merged['al_count2BOD'] = np.where(merged['al_countBOD'] == 0, 20, merged['al_count2BOD'])
merged['al_count2KIB'] = np.where(merged['al_countKIB'] == 0, 20, merged['al_count2KIB'])
merged['al_count2CAP'] = np.where(merged['al_countCAP'] == 0, 20, merged['al_count2CAP'])
merged['al_count2FOG'] = np.where(merged['al_countFOG'] == 0, 20, merged['al_count2FOG'])

final_df = pd.DataFrame([el for el in merged["pos"].values for _ in range(2)])
final_df = final_df.rename(columns={0: "pos"})
final_df["LOM"] = [ el for pair in zip(merged["al_countLOM"].values, merged["al_count2LOM"].values) for el in pair]
final_df["TER"] = [ el for pair in zip(merged["al_countTER"].values, merged["al_count2TER"].values) for el in pair]
final_df["BOD"] = [ el for pair in zip(merged["al_countBOD"].values, merged["al_count2BOD"].values) for el in pair]
final_df["KIB"] = [ el for pair in zip(merged["al_countKIB"].values, merged["al_count2KIB"].values) for el in pair]
final_df["CAP"] = [ el for pair in zip(merged["al_countCAP"].values, merged["al_count2CAP"].values) for el in pair]
final_df["FOG"] = [ el for pair in zip(merged["al_countFOG"].values, merged["al_count2FOG"].values) for el in pair]

final_df["LOM"] = final_df.LOM.round().astype(int)
final_df["TER"] = final_df.TER.round().astype(int)
final_df["BOD"] = final_df.BOD.round().astype(int)
final_df["KIB"] = final_df.KIB.round().astype(int)
final_df["CAP"] = final_df.CAP.round().astype(int)
final_df["FOG"] = final_df.FOG.round().astype(int)

#final_df = final_df.drop(['pos'], axis=1)

final_df.to_csv('/users/c/p/cpetak/WGS/angsd_new_noout/snpsfile_100_pos', index=False,sep="\t",header=False)
