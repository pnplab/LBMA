# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:14:58 2019

@author: asus
"""

from nilearn import datasets
from nilearn import image
import pandas as pd
import numpy as np
import re
import functools


# %%
def get_atlaslabels(csv_path):
    """ This function read and process a CSV file from its file path
     Return numpy array containing Atlas ROI's labels """

    atlas_csv = pd.read_csv(csv_path, sep=';', header=0)
    return atlas_csv.iloc[:,1].str.strip()


def read_atlasinfo(csv_path):
    """ This function read and process a CSV file from its file path
     Return numpy array containing Atlas ROI's name and labels """

    atlas_csv = pd.read_csv(csv_path, sep=';', header=0)
    return atlas_csv


def select_studies(file_df, contrast, method=None):
    studiesSelected = file_df[file_df['contrast'] == contrast]

    if (method != None):
        studiesSelected = studiesSelected[studiesSelected['method'].isin(method)]

    return studiesSelected


# %%
def find_origin_atlas(atlas_labels, region1_origin):
    """ This function return region's atlas name """

    name = atlas_labels.isin(region1_origin).any()
    origin_atlas = [*filter(name.get, name.index)]

    return origin_atlas


# %%

def check_region_inAtlas(atlas_labels, region_origin):
    """ This function checks if the region exists in one of the atlas  """
    origin_atlas = find_origin_atlas(atlas_labels, region_origin)
    return not origin_atlas


def check_Atlas(region_origin):
    """ This function checks if region_origin is empty """
    return not region_origin


def process_inputFile(input_file_df, atlas_labels):
    for i in range(len(input_file_df)):
        interaction = input_file_df.iloc[i]
        region1_origin = [interaction['region_1'].strip()]
        region2_origin = [interaction['region_2'].strip()]

        if check_Atlas(region1_origin) or check_Atlas(region2_origin):
            print('La ligne @ oublie des régions'.replace('@', str(i)))
            continue

        if (check_region_inAtlas(atlas_labels, region1_origin) or
                check_region_inAtlas(atlas_labels, region2_origin)):
            print('La ligne @ utilise un ou plusieurs labels inexistants'.replace('@', str(i)))

    return


# %%
@functools.lru_cache(maxsize=64)
def atlas_conversion(origin_path='AtlasNiftiFile.nii', target_path='NetworkNiftiFile.nii',
                     origin_CSV_path='Atlas.csv', target_CSV_path='Network.csv'):
    """ This function calculate and store the overlap between two atlas in a pandas dataframe """

    # ____________________________________________________________________________________
    # Preprocessing Steps
    # 1. Loading Nii filesg
    origin_nifti = image.load_img(origin_path)
    target_nifti = image.load_img(target_path)

    # 2.Resampling NIFTI IMAGE to MNI152 template
    template = datasets.load_mni152_template()
    origin_nifti_resampled = image.resample_to_img(origin_nifti, template, interpolation='nearest', copy=True,
                                                   order='F', clip=True, fill_value=0)
    target_nifti_resampled = image.resample_to_img(target_nifti, template, interpolation='nearest', copy=True,
                                                   order='F', clip=True, fill_value=0)

    # 3. Create 3D array containing voxel intensity
    origin_roi_array = np.round(origin_nifti_resampled.get_data())
    target_roi_array = np.round(target_nifti_resampled.get_data())

    # 3. Loading CSV file and processing to generate INDICES numpy array
    origin_csv_to_array = read_atlasinfo(origin_CSV_path).to_numpy()
    target_csv_to_array = read_atlasinfo(target_CSV_path).to_numpy()
    origin_roi_index = origin_csv_to_array[:, 0]
    target_roi_index = target_csv_to_array[:, 0]
    origin_roi_labels = [x.strip() for x in origin_csv_to_array[:, 1]]
    target_roi_labels = [x.strip() for x in target_csv_to_array[:, 1]]

    # 4.Intersection calculation
    intersection_matrix = np.zeros((len(origin_roi_index), len(target_roi_index)))

    for index1 in range(0, len(origin_roi_index)):
        region_origin_mask = (origin_roi_array == int(origin_roi_index[index1]))
        region_origin_size = np.count_nonzero(region_origin_mask)

        for index2 in range(0, len(target_roi_index)):
            region_target_mask = (target_roi_array == int(target_roi_index[index2]))
            inter_mask = np.logical_and(region_origin_mask, region_target_mask)
            intersection_matrix[index1, index2] = (np.count_nonzero(inter_mask) / region_origin_size)

    # 5. Put results in Dataframe
    overlap_rate = pd.DataFrame(intersection_matrix, index=origin_roi_labels, columns=np.transpose(target_roi_labels))

    return overlap_rate

#%%
def convert_region_to_net(region_origin, atlas_info, atlas_labels, target_atlas, large_atlas,
                          critere_conv, critere_conv_TD):
    """ convert a single origin region into the target atlas """

    origin_atlas = find_origin_atlas(atlas_labels, region_origin)[0]

    origin_atlas_path = atlas_info.loc['nii_path', origin_atlas]
    origin_atlas_csv = atlas_info.loc['csv_path', origin_atlas]
    target_atlas_path = atlas_info.loc['nii_path', target_atlas]
    target_atlas_csv = atlas_info.loc['csv_path', target_atlas]

    table_conversion = atlas_conversion(origin_path=origin_atlas_path,
                                        target_path=target_atlas_path,
                                        origin_CSV_path=origin_atlas_csv,
                                        target_CSV_path=target_atlas_csv)

    if any(atlas in origin_atlas for atlas in large_atlas):
        region1_target = table_conversion.loc[region_origin]
        region_target_temp = region1_target[region1_target > critere_conv_TD].any()
        region_target = [*filter(region_target_temp.get, region_target_temp.index)]

    elif critere_conv == 'max':
        region_target = table_conversion.loc[region_origin].idxmax(axis=1).tolist()

    else:
        region1_target = table_conversion.loc[region_origin]
        region_target_temp = region1_target[region1_target > critere_conv].any()
        region_target = [*filter(region_target_temp.get, region_target_temp.index)]
        if not region_target:
            region_target = table_conversion.loc[region_origin].idxmax(axis=1).tolist()

    return region_target

# %%
def convert_in_target_atlas(origin_data, table_conversion_recap, atlas_info, atlas_labels, target_atlas, large_atlas,
                            metho_intra, critere_conv, critere_conv_TD):
    """ Return the origin data table with region converted into the target atlas """

    data_converted = pd.DataFrame(columns=origin_data.columns)

    for i in range(len(origin_data)):
        interaction = origin_data.iloc[i]
        metho_interaction = interaction['method']

        region1_origin = [interaction['region_1'].strip()]
        region2_origin = [interaction['region_2'].strip()]

        if (check_Atlas(region1_origin) or check_Atlas(region2_origin) or check_region_inAtlas(atlas_labels,
                                                                                               region1_origin) or check_region_inAtlas(
            atlas_labels, region2_origin)):
            continue
        # _____________________________________________________________________________________________________________

        if table_conversion_recap['origin'].isin(region1_origin).any():
            region1_target = table_conversion_recap['target'].loc[
                table_conversion_recap['origin'] == region1_origin[0]].tolist()

        else:
            region1_target = convert_region_to_net(region1_origin, atlas_info, atlas_labels, target_atlas, large_atlas,
                                                   critere_conv, critere_conv_TD)

            for l in range(len(region1_target)):
                table_conversion_recap = table_conversion_recap.append(
                    pd.DataFrame([[region1_origin[0], region1_target[l]]], columns=table_conversion_recap.columns),
                    ignore_index=True)
        # _____________________________________________________________________________________________________________

        if table_conversion_recap['origin'].isin(region2_origin).any():
            region2_target = table_conversion_recap['target'].loc[
                table_conversion_recap['origin'] == region2_origin[0]].tolist()

        else:
            region2_target = convert_region_to_net(region2_origin, atlas_info, atlas_labels, target_atlas, large_atlas,
                                                   critere_conv, critere_conv_TD)
            for l in range(len(region2_target)):
                table_conversion_recap = table_conversion_recap.append(
                    pd.DataFrame([[region2_origin[0], region2_target[l]]], columns=table_conversion_recap.columns),
                    ignore_index=True)
        # ____________________________________________________________________________________________________________

        if metho_interaction in metho_intra:
            for reg1, reg2 in zip(region1_target, region2_target):
                row = interaction.copy()
                row.loc[['region_1', 'region_2']] = [reg1, reg2]
                data_converted = data_converted.append(row)
        else:
            for j in range(len(region1_target)):
                for k in range(len(region2_target)):
                    row = interaction.copy()
                    row.loc[['region_1', 'region_2']] = [region1_target[j], region2_target[k]]
                    data_converted = data_converted.append(row)

    target_file = data_converted.copy().reset_index(drop=True)
    region_sorted = target_file[['region_1', 'region_2']].apply(np.sort, axis=1)
    target_file[['region_1', 'region_2']] = pd.DataFrame(region_sorted.tolist())
    target_file_without_dup = target_file.drop_duplicates(ignore_index=True)

    return target_file_without_dup, table_conversion_recap


# %%
def fill_binary_matrix(input_df, empty_binary_matrix):
    network_linked_results1 = empty_binary_matrix.copy()

    for i in range(len(input_df)):
        interaction = input_df.iloc[i]
        study_number = interaction['study']

        region1_target = interaction['region_1'].strip()
        region2_target = interaction['region_2'].strip()

        if region1_target == region2_target:
            col = [','.join([region1_target, region2_target])]
        else:
            col = [x for x in empty_binary_matrix.columns if
                   (re.search(region1_target, x) and re.search(region2_target, x))]

        network_linked_results1.loc[study_number, col] = 1

    return network_linked_results1
