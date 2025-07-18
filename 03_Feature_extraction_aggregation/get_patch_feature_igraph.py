# coding=utf-8
"""
igraph提取patch拓扑特征，输出为csv文件
"""

from collections import defaultdict
import multiprocessing as mp
from scipy import stats
from scipy.spatial import cKDTree
import numpy as np
import igraph as ig
import json
import time
import os
import multiprocessing as mp
import pandas as pd
import sys

try:
    mp.set_start_method('spawn')
except:
    pass

# noinspection PyTypeChecker
def getstats(input_list, name):

    result = defaultdict(list)
    input_list_valid = np.ma.masked_invalid(input_list)

    result['min' + name] = np.min(input_list_valid)
    result['max' + name] = np.max(input_list_valid)
    result['mean' + name] = np.mean(input_list_valid)
    result['std' + name] = np.std(input_list_valid)
    result['skewness' + name] = float(stats.skew(input_list_valid, nan_policy='omit'))
    result['kurtosis' + name] = float(stats.kurtosis(input_list_valid, nan_policy='omit'))

    return result

def getGraphDisKnnFeatures(name, disKnnList):

    result = defaultdict(list)
    result['edge'] = ['v1','v2', 'v3', 'v4', 'v5', 'all'] # 由短到长的5条边 + 全部边

    disKnnList[np.isinf(disKnnList)] = np.nan
    disKnnList_valid = np.ma.masked_invalid(disKnnList)

    # 边数目特征
    result['NEdge'] += np.sum(disKnnList_valid.mask, axis=0).tolist()
    result['NEdge'].append(np.sum(disKnnList_valid.mask))

    # 边长度特征
    result['minEdgeLength'] += np.min(disKnnList_valid, axis=0).tolist()
    result['minEdgeLength'].append(np.min(disKnnList_valid))

    result['maxEdgeLength'] += np.max(disKnnList_valid, axis=0).tolist()
    result['maxEdgeLength'].append(np.max(disKnnList_valid))

    result['meanEdgeLength'] += np.mean(disKnnList_valid, axis=0).tolist()
    result['meanEdgeLength'].append(np.mean(disKnnList_valid))

    result['stdEdgeLength'] += np.std(disKnnList_valid, axis=0).tolist()
    result['stdEdgeLength'].append(np.std(disKnnList_valid))

    result['skewnessEdgeLength'] += stats.skew(disKnnList, axis=0, nan_policy='omit').tolist()
    res = stats.skew(disKnnList, axis=None, nan_policy='omit')
    if isinstance(res, float):
        result['skewnessEdgeLength'].append(res)
    else:
        result['skewnessEdgeLength'].append(res.tolist())

    result['kurtosisEdgeLength'] += stats.kurtosis(disKnnList, axis=0, nan_policy='omit').tolist()
    res = stats.kurtosis(disKnnList, axis=None, nan_policy='omit')
    if isinstance(res, float):
        result['kurtosisEdgeLength'].append(res)
    else:
        result['kurtosisEdgeLength'].append(res.tolist())


    return result

# noinspection PyTypeChecker
def getSingleGraphFeatures(args):

    subgraph, cmd = args
    result = defaultdict(list)
    n = subgraph.vcount()

    # 图结构特征
    if cmd == 'Diameter':
        result['Diameter'] = subgraph.diameter()  # 直径
    elif cmd == 'Girth':
        result['Girth'] = subgraph.diameter()  # 周长
    elif cmd == 'Radius':
        result['Radius'] = subgraph.radius()  # 半径
    elif cmd == 'Assortativity':
        result['Assortativity'] = subgraph.assortativity_degree()  # 度相关性
    elif cmd == 'Density':
        result['Density'] = subgraph.density()  # 边密度
    elif cmd == 'Transitivity':
        result['Transitivity'] = subgraph.transitivity_undirected()  # 传递性
    elif cmd == 'Transitivity_avglocal':
        result['Transitivity_avglocal'] = subgraph.transitivity_undirected()  # 平均局部传递性

    # 图聚类特征
    elif cmd == 'Ncell':
        result['Ncell'] = n  # 细胞数目
    elif cmd == 'Nsubgraph':
        result['Nsubgraph'] = len(subgraph.decompose())  # 子图数目

    ## 太慢
    # elif cmd == 'Ncohesive':
    #     result['Ncohesive'] = len(subgraph.cohesive_blocks()) # 内聚块数目
    # elif cmd == 'Nclique':
    #     result['Nclique'] = subgraph.clique_number() # 集团数目
    # elif cmd == 'Nmotifs':
    #     result['Nmotifs'] = subgraph.motifs_randesu_no() # motif数目

    ## 容易报错
    # result['NclusterEdgeBetweenness'] = len(subgraph.community_edge_betweenness().as_clustering()) # 基于边中介度的聚类数目
    # result['NclusterFastgreedy'] = len(subgraph.community_fastgreedy().as_clustering()) # 基于最大贪婪算法的聚类数目
    # result['NclusterInfomap'] = len(subgraph.community_infomap()) # 基于infomap算法的聚类数目
    # result['NclusterLabelPropagation'] = len(subgraph.community_label_propagation())  # 基于label propagation算法的聚类数目
    # result['NclusterCommunityLeadingEigenvector'] = len(subgraph.community_leading_eigenvector()) # 基于leading eigenvector算法的聚类数目
    # result['NclusterCommunityLeiden'] = len(subgraph.community_leiden(objective_function='modularity')) # 基于leiden算法的聚类数目
    # result['NclusterCommunityWalktrap'] = len(subgraph.community_walktrap().as_clustering()) # 基于walktrap算法的聚类数目

    # 图节点特征
    elif cmd == 'Degree':
        result.update(getstats(subgraph.degree(), 'Degree'))
    elif cmd == 'Coreness':
        result.update(getstats(subgraph.coreness(), 'Coreness'))
    elif cmd == 'Eccentricity':
        result.update(getstats(subgraph.eccentricity(), 'Eccentricity'))
    elif cmd == 'HarmonicCentrality':
        result.update(getstats(subgraph.harmonic_centrality(), 'HarmonicCentrality'))
    elif cmd == 'Closeness':
        result.update(getstats(subgraph.closeness(), 'Closeness'))
    elif cmd == 'Betweenness':
        result.update(getstats(subgraph.betweenness(), 'Betweenness'))

    return result


def getGraphCenterFeatures(subgraph: ig.Graph):
    norm_cmds = ['Diameter', 'Girth', 'Radius', 'Assortativity', 'Density', 'Transitivity', 'Transitivity_avglocal', 'Ncell', 'Nsubgraph', 'Ncohesive', 'Nclique', 'Nmotifs',
                 'Degree', 'Coreness'] # 单线程计算特征
    multi_cmds = ['Eccentricity', 'HarmonicCentrality', 'Closeness', 'Betweenness']  # 多线程计算特征

    result = defaultdict(list)

    # 单线程特征
    for cmd in norm_cmds:
        args = [subgraph, cmd]
        ans = getSingleGraphFeatures(args)
        result.update(ans)

    # 多线程特征
    if subgraph.vcount() > 10000:  # Huge graph, use multiprocessing
        args = [[subgraph, cmd] for cmd in multi_cmds]
        with mp.Pool() as p:
            ans = p.map(getSingleGraphFeatures, args)
            for q_info in ans:
                result.update(q_info)
    else:  # Small graph, directly calculate
        for cmd in multi_cmds:
            args = [subgraph, cmd]
            ans = getSingleGraphFeatures(args)
            result.update(ans)

    return result


def constructGraphFromDict(
        nucleusInfo: dict, distanceThreshold: float,
        knn_n: int = 5, level: int = 0, offset=np.array([0, 0])
):
    r"""Construct graph from nucleus information dictionary

    Parameters
    ----------
    nucleusInfo : dict
        'mag': int
            magnification of the result
        'nuc': dict
            nucleus information
            'nuclei ID' : dict
                note that ID generated from HoverNet is not continuous
                'bbox' : list
                    [[left, top], [right, bottom]]
                'centroid' : list
                    [column, row]
                'contour' : list, from cv2.findContours
                    [[column1, row1], [column2, row2], ... ]
                'type_prob' : float
                    The probability of current nuclei belonging to type 'type'
                'type' : int

    distanceThreshold : maximum distance in magnification of 40x

    typeDict : dict
        "0" : "nolabe"
        "1" : "neopla"
        "2" : "inflam"
        "3" : "connec"
        "4" : "necros"
        "5" : "no-neo"

    cellSize : int, odd
        size of cell cropped for extracting GLCM features

    level : int
        level for reading WSI
        0 : 40x
        1 : 20x
        ...

    Returns
    -------
    graph :

    """
    typeDict2 = {
        'neolabe': 0, # 无
        'neopla': 1, # 肿瘤
        'inflame': 2, # 炎症
        'connect': 3, # 基质
        'necros': 4, # 坏死
        'normal': 5 # 正常
    }
    mag = nucleusInfo['patch_mag']
    distanceThreshold = distanceThreshold / (40.0 / mag)  # 细胞之间的最大像素距离

    t3 = time.time()

    # 1、从json中提取信息
    centroids, types = nucleusInfo['cell_centroid'], nucleusInfo['cell_type']
    vertex_len = nucleusInfo['cell_num']

    # 创建全局图
    globalGraph = ig.Graph()
    names = [str(i) for i in range(vertex_len)]
    globalGraph.add_vertices(vertex_len, attributes={
        'name': names, 'Centroid': centroids, 'CellType': types})  # 向图中添加细胞作为节点

    # 2、创建细胞之间的边
    nolabeIDs, neoplaIDs, inflamIDs, connecIDs, necrosIDs, normalIDs = \
        [np.where(np.array(types) == i)[0].tolist() for i in range(6)] # 6种细胞的索引

    edge_info = pd.DataFrame({'source': [], 'target': [], 'length': [], 'featype': []})  # 边信息表

    # Neopla->T, Inflam->I, Connec->S, Normal->N
    featype_dict = {'T-T': [neoplaIDs, neoplaIDs],
                    'I-I': [inflamIDs, inflamIDs],
                    'S-S': [connecIDs, connecIDs],
                    # 'N-N': [normalIDs, normalIDs],
                    'T-I': [neoplaIDs, inflamIDs],
                    'T-S': [neoplaIDs, connecIDs],
                    # 'T-N': [neoplaIDs, normalIDs],
                    'I-S': [inflamIDs, connecIDs],
                    # 'I-N': [inflamIDs, normalIDs],
                    # 'S-N': [connecIDs, normalIDs]
                    }

    feats_edge = {}
    feats_graph = {}
    feats_block = {}

    for featype, featype_index_list in zip(featype_dict.keys(),
                                           featype_dict.values()):  # 遍历边类型，以及对应的细胞索引

        print(f'Getting {featype} graph feature')
        print(f'---Creating edges')
        # Treat neopla and normal as the same cell type by making the same cellTypeMark,
        # and delete the edge between vertexs which have the same cellTypeMark

        pairs = np.array([]).reshape((0, 3))  # crz修改

        disKnnList = np.array([]).reshape((0, knn_n))
        subgraph_names = []
        featype_index = []
        for index in featype_index_list:
            featype_index += index

        if (len(featype_index) == 0): # patch中没有target and source细胞类型，不构建网络
            pairs = []
        else:
            for src, i_src in enumerate(featype_index_list):  # source细胞类型和对应的索引
                for tar, i_tar in enumerate(featype_index_list):  # target细胞类型和对应的索引
                    if src != tar:
                        centroid_tar = globalGraph.induced_subgraph(i_tar).vs['Centroid']  # 根据细胞索引，提取target的子图，并获得像素坐标
                        centroid_src = globalGraph.induced_subgraph(i_src).vs['Centroid']  # 提取包含source细胞的子图，并获得像素坐标
                        n_tar = len(i_tar)  # target细胞数目
                        n_src = len(i_src)  # source细胞数目
                        if n_tar == 0 or n_src == 0:  # patch中没有target/source细胞类型，不构建网络
                            pairs = []
                        else:
                            tree = cKDTree(centroid_tar)  # 根据target坐标构建CKD树
                            if i_src == i_tar:  # 如果source和target的坐标均相同（T-T/S-S/I-I）
                                disknn, vertex_index = tree.query(centroid_src, k=knn_n + 1,
                                                                  distance_upper_bound=distanceThreshold,
                                                                  p=2)  # 寻找距离每个source最近的6个target，使用欧几里得距离，距离上限为40个像素，超出设置为inf，相应的索引被设置为索引最大值
                                disknn = disknn[..., -knn_n:]  # k个最近邻的距离（删掉距离为0的第一列）
                                vertex_index = vertex_index[..., -knn_n:]  # k个最近邻的索引（删掉距离为0的第一列）
                            else:  # 如果source和target的坐标不同
                                disknn, vertex_index = tree.query(centroid_src, k=knn_n,
                                                                  distance_upper_bound=distanceThreshold,
                                                                  p=2)  # 寻找距离每个source最近的5个target

                            knn_mask = vertex_index != n_tar  # delete the vertex whose distance upper bound （删除距离超过上限的点）
                            v_src = np.tile(np.array(i_src, dtype='str').reshape((n_src, -1)), (1, knn_n))[
                                knn_mask]  # 创建source的一维数组
                            v_tar = np.array(i_tar, dtype='str')[vertex_index[knn_mask]]  # 将target压缩一维数组
                            v_len = disknn[knn_mask]  # 边距离

                            pairs = np.concatenate([pairs, np.stack((v_src, v_tar, v_len), axis=1)],
                                                   axis=0)  # 带有距离的边关系表、crz修改
                            disKnnList = np.concatenate([disKnnList, disknn], axis=0)  # 距离表
                            subgraph_names += i_src

        if len(pairs) == 0: # 没有边

            print(f'---Getting DisKnn features')
            feats = {
                'edge': ['v1', 'v2', 'v3', 'v4', 'v5', 'all'],
                'NEdge': [np.nan] * 6,
                'minEdgeLength': [np.nan] * 6,
                'maxEdgeLength': [np.nan] * 6,
                'meanEdgeLength': [np.nan] * 6,
                'stdEdgeLength': [np.nan] * 6,
                'skewnessEdgeLength': [np.nan] * 6,
                'kurtosisEdgeLength': [np.nan] * 6
            }

            for k, v in zip(feats.keys(), feats.values()):
                if k != 'edge':
                    for i in range(6):
                        feats_edge[featype + '_edge_' + feats['edge'][i] + '_' + k] = v[i]

            print(f'---Getting GraphCenter features')
            feats = {
             'Diameter': np.nan,
             'Girth': np.nan,
             'Radius': np.nan,
             'Assortativity': np.nan,
             'Density': np.nan,
             'Transitivity': np.nan,
             'Transitivity_avglocal': np.nan,
             'Ncell': np.nan,
             'Nsubgraph': np.nan,
             # 'Ncohesive': np.nan,
             # 'Nclique': np.nan,
             # 'Nmotifs': np.nan,
             # 'NclusterEdgeBetweenness': np.nan,
             # 'NclusterFastgreedy': np.nan,
             # 'NclusterInfomap': np.nan,
             # 'NclusterLabelPropagation': np.nan,
             # 'NclusterCommunityLeadingEigenvector': np.nan,
             # 'NclusterCommunityLeiden': np.nan,
             # 'NclusterCommunityWalktrap': np.nan,
             'minDegree': np.nan,
             'maxDegree': np.nan,
             'meanDegree': np.nan,
             'stdDegree': np.nan,
             'skewnessDegree': np.nan,
             'kurtosisDegree': np.nan,
             'minCloseness': np.nan,
             'maxCloseness': np.nan,
             'meanCloseness': np.nan,
             'stdCloseness': np.nan,
             'skewnessCloseness': np.nan,
             'kurtosisCloseness': np.nan,
             'minBetweenness': np.nan,
             'maxBetweenness': np.nan,
             'meanBetweenness': np.nan,
             'stdBetweenness': np.nan,
             'skewnessBetweenness': np.nan,
             'kurtosisBetweenness': np.nan,
             'minCoreness': np.nan,
             'maxCoreness': np.nan,
             'meanCoreness': np.nan,
             'stdCoreness': np.nan,
             'skewnessCoreness': np.nan,
             'kurtosisCoreness': np.nan,
             'minEccentricity': np.nan,
             'maxEccentricity': np.nan,
             'meanEccentricity': np.nan,
             'stdEccentricity': np.nan,
             'skewnessEccentricity': np.nan,
             'kurtosisEccentricity': np.nan,
             'minHarmonicCentrality': np.nan,
             'maxHarmonicCentrality': np.nan,
             'meanHarmonicCentrality': np.nan,
             'stdHarmonicCentrality': np.nan,
             'skewnessHarmonicCentrality': np.nan,
             'kurtosisHarmonicCentrality': np.nan
            }

            for k, v in zip(feats.keys(), feats.values()):
                feats_graph[featype + '_graph_' + k] = v

            continue

        else: # 有边

            subgraph = globalGraph.induced_subgraph(set(featype_index))  # 提取包含source和target的子图
            subgraph.add_edges(pairs[:, 0:2])  # 向子图中添加边
            multiple_edge_index = np.where(np.array(subgraph.is_multiple()))[0].tolist()
            multiple_edge = subgraph.es[multiple_edge_index]
            subgraph.delete_edges(multiple_edge)  # delete the multiple edge 去掉重复的边

            subgraph_edge = subgraph.get_edge_dataframe()  # 提取边
            subgraph_vname = subgraph.get_vertex_dataframe()['name']  # 提取点
            subgraph_edge['source'] = [subgraph_vname[a] for a in subgraph_edge['source'].values]  # 修改边的source名
            subgraph_edge['target'] = [subgraph_vname[a] for a in subgraph_edge['target'].values]  # 修改边的target名
            subgraph_edge.insert(subgraph_edge.shape[1], 'featype', featype)  # 加入边类型

            pairs = np.delete(pairs, multiple_edge_index, axis=0)  # 去掉重复边 crz修改
            subgraph_edge.insert(subgraph_edge.shape[1], 'length', pairs[:, 2])
            edge_info = pd.concat([edge_info, subgraph_edge], sort=True)  # 边信息的汇总表

            print(f'---Getting DisKnn features')
            feats = getGraphDisKnnFeatures(subgraph_names, disKnnList)  # 边特征
            for k, v in zip(feats.keys(), feats.values()):
                if k != 'edge':
                    for i in range(6):
                        feats_edge[featype + '_edge_' + feats['edge'][i] + '_' + k] = v[i]

            print(f'---Getting GraphCenter features')
            feats = getGraphCenterFeatures(subgraph) # 网络特征
            for k, v in zip(feats.keys(), feats.values()):
                feats_graph[featype + '_graph_' + k] = v

    print(f"{'Graph features cost':#^40s}, {time.time() - t3:*^10.2f}")

    # Stroma blocker
    # For each inflam node, add it to neoplaAddConnecGraph, compute blocker and delete
    # ! Why select the maximum subgraph, if distanceThreshold wasn't set appropriately, shortestPathsLymCancer would be empty
    t4 = time.time()

    if (globalGraph.vcount() != 0):

        centroid_T = globalGraph.induced_subgraph(neoplaIDs).vs['Centroid']
        centroid_I = globalGraph.induced_subgraph(inflamIDs).vs['Centroid']
        centroid_S = globalGraph.induced_subgraph(connecIDs).vs['Centroid']

        if len(centroid_T) != 0 and len(centroid_I) != 0 and len(centroid_S) != 0:

            Ttree = cKDTree(centroid_T)
            STree = cKDTree(centroid_S)
            dis, pairindex_T = Ttree.query(centroid_I, k=1)
            paircentroid_T = np.array(centroid_T)[pairindex_T]
            blocker = []

            for Tcoor, Icoor, r in zip(centroid_I, paircentroid_T, dis):
                set1 = set(STree.query_ball_point(Tcoor, r))
                set2 = set(STree.query_ball_point(Icoor, r))
                blocker.append(len(set1 & set2))

            feats_block['sumBlock'] = np.sum(blocker)
            feats_block.update(getstats(blocker, name='Block'))

        else:

            feats_block = {
                'sumBlock': np.nan,
                'minBlock': np.nan,
                'maxBlock': np.nan,
                'meanBlock': np.nan,
                'stdBlock': np.nan,
                'skewnessBlock': np.nan,
                'kurtosisBlock': np.nan
            }

    else:

        feats_block = {
            'sumBlock': np.nan,
            'minBlock': np.nan,
            'maxBlock': np.nan,
            'meanBlock': np.nan,
            'stdBlock': np.nan,
            'skewnessBlock': np.nan,
            'kurtosisBlock': np.nan
        }

    print(f"{'stroma blocker cost':#^40s}, {time.time() - t4:*^10.2f}")

    feats_all = {}
    feats_all.update(feats_edge)
    feats_all.update(feats_graph)
    feats_all.update(feats_block)

    return feats_all


def getFeatureCSV(args):

    distanceThreshold, k, level, json_path, output_path, sample_name = args
    json_file = json_path + sample_name + '.json'

    # 读入patch json
    with open(json_file) as fp:
        patchInfo = json.load(fp)

    # 特征文件路径
    output_path_patch = os.path.join(output_path, sample_name)
    csvfile = os.path.join(output_path_patch + '.csv')

    print(sample_name)

    # patch特征提取和写出
    ## 逐行写出
    # with open(csvfile, mode='w', newline='') as file:
    #     len_patch = len(patchInfo.keys())
    #     for patch in range(len_patch):
    #         # print(patch)
    #         nucleusInfo = patchInfo[str(patch)]
    #         sys.stdout = open(os.devnull, 'w')  # 关闭print输出
    #         feats = constructGraphFromDict(nucleusInfo, distanceThreshold, k, level)
    #         sys.stdout = sys.__stdout__ # 打开print输出
    #         if patch == 0:
    #             writer = csv.writer(file)
    #             writer.writerow(['Patch_index'] + list(feats.keys()))
    #         writer.writerow([str(patch)] + list(feats.values()))

    ## 合并写出 - 方便排错
    len_patch = len(patchInfo.keys())
    for patch in range(len_patch):
        nucleusInfo = patchInfo[str(patch)]
        # print(patch)
        sys.stdout = open(os.devnull, 'w')  # 关闭print输出
        feats = constructGraphFromDict(nucleusInfo, distanceThreshold, k, level)
        sys.stdout = sys.__stdout__ # 打开print输出
        if (patch == 0):
            result = pd.DataFrame.from_dict(feats, orient='index')
        else:
            result[str(patch)] = pd.Series(feats)
    result.to_csv(csvfile, index=True, header=True)


if __name__ == "__main__":

    distanceThreshold = 100 # 像素距离阈值
    level = 1 # 20倍放大
    k = 5 # 5个最近邻

    json_path = '/home/ljc/0_project/0_ESCC/data/patch/xiaoshu_ren_2048/json/'
    output_path = '/home/ljc/0_project/0_ESCC/data/patch/xiaoshu_ren_2048/feature/'

    # 需要处理的文件list
    json_list = [os.path.splitext(file_name)[0] for file_name in os.listdir(json_path)]
    exist_list = [os.path.splitext(file_name)[0] for file_name in os.listdir(output_path)]
    todo_list = list(set(json_list) - set(exist_list)) # 取差集

    # 并行处理
    args = [[distanceThreshold, k, level, json_path, output_path, sample_name] for sample_name in todo_list]  # 打包参数
    pool = mp.Pool(processes=4)  # 开启并行核心
    pool.map(getFeatureCSV, args)

    # # 串行处理
    # for samp in todo_list:
    #     args = [distanceThreshold, k, level, json_path, output_path, samp]
    #     getFeatureCSV(args)