import importlib
import region_interaction_lib as lib

importlib.reload(lib)

import os
import pandas as pd
import itertools

pd.options.mode.chained_assignment = None

# %% Input data

# Path of the folder containing all the information on the atlases
atlas_folder_path = os.path.abspath(r'../atlas')

# Read Atlas folder
atlas_sub_files = os.listdir(atlas_folder_path)
atlas_names = sorted(set(map(lambda x: x.split('.')[0], atlas_sub_files)))

# Separation of csv and nii files then storage in dataframe
csv_files_path = sorted(list(filter(lambda x: '.csv' in x, atlas_sub_files)))
nii_files_path = sorted(list(filter(lambda x: '.nii' in x, atlas_sub_files)))

atlas_info = pd.DataFrame(data=[nii_files_path, csv_files_path],
                          index=['nii_path', 'csv_path'],
                          columns=atlas_names)

atlas_info = atlas_info.applymap(lambda x: os.path.join(atlas_folder_path, x))

# %% Creation of a dataframe containing all the labels of the original atlases
atlas_labels = pd.DataFrame()

for atlas_name in atlas_info.columns:
    csv_path = atlas_info.loc['csv_path', atlas_name]
    labels = lib.get_atlaslabels(csv_path)
    atlas_labels = pd.concat([atlas_labels, labels], axis=1)

atlas_labels.columns = atlas_info.columns

# %% Target atlas
target_atlas = 'aal3_bilateral'

if 'aal' in target_atlas:
    large_atlas = ['TD', 'whole', 'Yeo']
else:
    large_atlas = ['TD', 'whole']
# %%
# Reading the csv file containing the effects and tested

input_folder_path = os.path.abspath(r'../input')
# input_csv_path = list(filter(lambda x: target_atlas in x, os.listdir(input_folder_path)))[0]
input_csv_path = 'input_cole12_n428.csv'
input_file = pd.read_csv(os.path.join(input_folder_path, input_csv_path), sep=';', header=0)

# %% Input variables

critere_conv = 'max'
critere_conv_TD = 0.001
group_list = ['SZ-HC', 'BD-HC', 'MDD-HC']
metho_intra = ['ICA', 'VMHC', 'ALFF', 'ReHo']
metho_inter = ['SBVW', 'SBCW', 'STR']
metho = metho_intra + metho_inter
# Initialisation of table which contains origin label -> target label
table_conversion_recap = pd.DataFrame(columns=['origin', 'target'])

# %%
# Results
results_csv_folder_path = os.path.abspath(r'../results')
results_folder_created = '_'.join([target_atlas, input_csv_path.replace('.csv', ''), str(critere_conv)])

results_folder_path = os.path.join(results_csv_folder_path, results_folder_created)

if not os.path.exists(results_folder_path):
    os.makedirs(results_folder_path)

# %%

target_file_without_dup, table_conversion = lib.convert_in_target_atlas(input_file, table_conversion_recap, atlas_info,
                                                                        atlas_labels, target_atlas, large_atlas,
                                                                        metho_intra, critere_conv=critere_conv,
                                                                        critere_conv_TD=critere_conv_TD)
target_file_without_dup.to_csv(
    os.path.join(results_folder_path, 'target_file_atlas.csv'.replace('atlas', target_atlas)),
    index=False)

if 'aal' in target_atlas:
    networks_CAB_NP = pd.read_csv(r'../input/conversion_CAB-NP_DMN-CON-FPN_aal3.csv', header=0)
    target_atlas_labels = networks_CAB_NP['aal3'].to_list()
    target_file_filtered = target_file_without_dup[(target_file_without_dup['region_1'].isin(target_atlas_labels))
                                                   & target_file_without_dup['region_2'].isin(target_atlas_labels)]

    target_atlas_labels_comb = [','.join(i) for i in
                                list(itertools.combinations_with_replacement(target_atlas_labels, 2))]

    target_file_filtered.to_csv(
        os.path.join(results_folder_path, 'target_file_atlas_filtered.csv'.replace('atlas', target_atlas)),
        index=False)
    target_file = target_file_filtered
else:
    target_atlas_labels = atlas_labels[target_atlas].dropna()
    target_atlas_labels_comb = [','.join(i) for i in
                                list(itertools.combinations_with_replacement(target_atlas_labels, 2))]
    target_file = target_file_without_dup

# %%

tested_file = target_file[target_file['direction'] == 0]
pos_effect_file = target_file[target_file['direction'] == 1]
neg_effect_file = target_file[target_file['direction'] == -1]

# %%
for group in group_list:
    # Study selection from tested dataframe
    tested_input = lib.select_studies(tested_file, group, method=metho)
    list_study_number = tested_input['study'].unique()
    network_linked_results = pd.DataFrame(0, index=list_study_number, columns=target_atlas_labels_comb)

    # Effects selection
    pos_effect_input = lib.select_studies(pos_effect_file, group, method=metho)
    neg_effect_input = lib.select_studies(neg_effect_file, group, method=metho)

    # Baseline file treatment
    tested_target = lib.fill_binary_matrix(tested_input, network_linked_results)
    tested_target.to_csv(os.path.join(results_folder_path,
                                      'tested_target_group.csv'.replace('group', group).replace('target',
                                                                                                target_atlas)),
                         index_label=None)

    # Effects file treatment
    pos_effect_target = lib.fill_binary_matrix(pos_effect_input, network_linked_results)
    (pos_effect_target * tested_target).to_csv(os.path.join(results_folder_path,
                                                            'effects_target_group_hyper.csv'.replace('group',
                                                                                                     group).replace(
                                                                'target', target_atlas)), index_label=None)

    neg_effect_target = lib.fill_binary_matrix(neg_effect_input, network_linked_results)
    (neg_effect_target * tested_target).to_csv(os.path.join(results_folder_path,
                                                            'effects_target_group_hypo.csv'.replace('group',
                                                                                                    group).replace(
                                                                'target', target_atlas)), index_label=None)
