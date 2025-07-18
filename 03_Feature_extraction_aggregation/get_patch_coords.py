# coding=utf-8
"""
提取每个patch中的细胞坐标，输出为json文件
添加并行化
"""

import h5py
import json
import os
import multiprocessing as mp

def getPatchJson(args):

    seg_path, h5_path, json_path, samp, patch_size = args

    path_patch = os.path.join(h5_path, samp + '.h5')
    path_cell = os.path.join(seg_path, samp + '.json')
    path_output = os.path.join(json_path, samp + '.json')

    # 读入patch坐标
    h5 = h5py.File(path_patch, 'r')
    coords = h5['coords'][:]

    # 读入细胞坐标
    bbox_list_wsi = []
    centroid_list_wsi = []
    contour_list_wsi = []
    type_list_wsi = []

    with open(path_cell) as json_file:
        data = json.load(json_file)
        mag_info = data['mag']
        nuc_info = data['nuc']
        for inst in nuc_info:
            inst_info = nuc_info[inst]
            inst_centroid = inst_info['centroid']
            centroid_list_wsi.append(inst_centroid)
            inst_contour = inst_info['contour']
            contour_list_wsi.append(inst_contour)
            inst_bbox = inst_info['bbox']
            bbox_list_wsi.append(inst_bbox)
            inst_type = inst_info['type']
            type_list_wsi.append(inst_type)

    json_file.close()

    print(samp)

    # 提取每个patch中的细胞坐标
    result = {}
    len_coords = len(coords[:])

    for i in range(len_coords):
        row = coords[:][i]
        index = []
        x0, y0 = row  # patch的起始坐标
        x1, y1 = x0 + patch_size, y0 + patch_size  # patch的终点坐标

        for cell in centroid_list_wsi:
            x, y = cell  # 细胞质心坐标
            if x0 < x < x1 and y0 < y < y1:
                index.append(True)
            else:
                index.append(False)

        ind = [x for x, value in enumerate(index) if value]  # 属于当前patch的细胞索引
        res = {
            "patch_coords": [int(x0), int(y0)],
            "patch_mag": mag_info,
            "patch_size": patch_size,
            "cell_num": len(ind),
            "cell_centroid": [centroid_list_wsi[i] for i in ind],
            "cell_type": [type_list_wsi[i] for i in ind]
        }
        result[str(i)] = res

    # 写出patch的json文件
    with open(path_output, "w") as json_file:
        json.dump(result, json_file)


if __name__ == "__main__":

    seg_path = '/home/ljc/0_project/0_ESCC/data/segment/xiaoshu_ren/json'
    h5_path = '/home/ljc/0_project/0_ESCC/data/patch/xiaoshu_ren_2048/patches'
    json_path = '/home/ljc/0_project/0_ESCC/data/patch/xiaoshu_ren_2048/json'
    patch_size = 2048

    # 需要处理的文件list
    json_list = [os.path.splitext(file_name)[0] for file_name in os.listdir(seg_path)]
    exist_list = [os.path.splitext(file_name)[0] for file_name in os.listdir(json_path)]
    todo_list = list(set(json_list) - set(exist_list))  # 取差集

    # # 并行处理
    # pool = mp.Pool(processes=4) # 开启并行核心
    # args = [[seg_path, h5_path, json_path, samp, patch_size] for samp in todo_list] # 打包参数
    # pool.map(getPatchJson, args)

    # 串行处理
    todo_list = list(set([os.path.splitext(file_name)[0] for file_name in os.listdir(seg_path)]) - set(
        [os.path.splitext(file_name)[0] for file_name in os.listdir(json_path)]))
    for samp in todo_list:
        getPatchJson([seg_path, h5_path, json_path, samp, patch_size])