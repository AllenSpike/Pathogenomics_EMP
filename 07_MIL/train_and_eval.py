import numpy as np
import torch
from PIL import Image
import slideflow as sf
import slideflow.tfrecord
from slideflow.slide import qc
from slideflow.mil import mil_config
from mytools import img_utils
import os
import random
from sklearn.metrics import classification_report, roc_auc_score, roc_curve
from sklearn.preprocessing import label_binarize
import re
import pandas as pd
import matplotlib.pyplot as plt
from slideflow.mil.eval import predict_mil, generate_mil_features
import json

os.environ['SF_SLIDE_BACKEND'] = 'libvips'
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2,50))
qc = [qc.Gaussian(), qc.Otsu()]

### 创建项目 ####
# project = sf.create_project(
#     root='xj_tcga_gdph',
#     slides='/home/ljc/DataDisk/0_ESCC/XJ_TCGA_GDPH',
#     annotations='/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/annotations.csv'
# ) # 创建项目
# project.extract_tiles(qc = qc,
#                       tile_px = 299,
#                       tile_um = 302,
#                       normalizer ='reinhard_fast',
#                       num_threads = 20
# ) # qc, 提取patch，使用Reinhard-fast归一化
project = sf.load_project(root='xj_tcga_gdph') # 导入项目
dataset = project.dataset(tile_px=299, tile_um=302) # 创建数据集

#### 特征提取器 ####
resnet = sf.build_feature_extractor("resnet50_imagenet", tile_px=299, resize=True)
retccl = sf.build_feature_extractor('retccl', ckpt = '/home/ljc/0_project/0_ESCC/0_slideflow/model_weight/retccl.pth', resize=True)
ctranspath = sf.build_feature_extractor('ctranspath', resize=True)
phikon = sf.build_feature_extractor('phikon', weights='/home/ljc/0_project/0_ESCC/0_slideflow/model_weight/ibot_vit_base_pancan.pth', resize=True)
# uni = sf.build_feature_extractor('uni', weights='/home/ljc/0_project/0_ESCC/0_slideflow/model_weight/uni/pytorch_model.bin', resize=True) # uni训练有NAN
# giga = sf.build_feature_extractor('gigapath') # gigapath导入失败
# plip = sf.build_feature_extractor('plip') # plip导入失败

#### 特征提取 ####
project.generate_feature_bags(resnet, dataset, outdir='/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_resnet')
project.generate_feature_bags(retccl, dataset, outdir='/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_retccl')
project.generate_feature_bags(ctranspath, dataset, outdir='/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_ctranspath')
project.generate_feature_bags(phikon, dataset, outdir='/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_phikon')

#### 数据集划分 ####
random.seed(42)
splits = dataset.kfold_split(k=5, labels='category', preserved_site=False)
np.save('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/splits.npy', splits)

splits = np.load('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/splits.npy', allow_pickle=True)

#### 模型训练 ####
config = mil_config('clam_sb', lr=1e-4)
config = mil_config('clam_mb', lr=1e-4)
config = mil_config('attention_mil', lr=1e-4)
config = mil_config('transmil', lr=1e-4)

for train, val in splits:
    project.train_mil(
        config=config,
        attention_heatmaps=True,
        cmap='jet',
        interpolation='bicubic',
        bags='/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_ctranspath',
        outcomes='category',
        train_dataset=train,
        val_dataset=val
    )

#### 模型指标 ####
path = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil'
files = [file for file in os.listdir(path)]

# 提取指标
files = sorted([file for file in files if re.search('transmil', file)])
metrics = []
metrics_df = pd.DataFrame(columns=['accuracy', 'auc', 'macro avg f1-score', 'weighted avg f1-score'])

for i in range(len(files)):
    file = files[i]
    model_path = os.path.join(path, file)
    df_clam = project.evaluate_mil(
        model=model_path,
        bags='/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_ctranspath',
        outcomes='category',
        dataset=splits[i][1])

    y_true = df_clam.iloc[:, 1].values
    y_ture_bin = label_binarize(y_true, classes=[0, 1, 2])
    y_score = df_clam.iloc[:, 2:].values
    y_pred = df_clam.iloc[:, 2:].values.argmax(axis=1)

    fpr = dict()
    tpr = dict()
    for x in range(3):
        fpr[x], tpr[x], _ = roc_curve(y_ture_bin[:,x], y_score[:,x])
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(3)]))
    mean_tpr = np.zeros_like(all_fpr)
    for x in range(3):
        mean_tpr += np.interp(all_fpr, fpr[x], tpr[x])
    mean_tpr /= 3
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    auc = roc_auc_score(y_true, y_score, multi_class='ovr')

    metric = classification_report(y_true, y_pred, output_dict=True)
    metric.update({'fpr': fpr})
    metric.update({'tpr': tpr})
    metric.update({'auc': auc})

    metrics.append(metric)
    metrics_df.loc[i] = [
        metrics[i]['accuracy'],
        metrics[i]['auc'],
        metrics[i]['macro avg']['f1-score'],
        metrics[i]['weighted avg']['f1-score']
    ]
metrics_df.mean()

# 绘制ROC曲线
plt.figure()
for i in range(len(files)):
    fpr = metrics[i]['fpr']
    tpr = metrics[i]['tpr']
    auc = metrics[i]['auc']
    plt.plot(fpr["macro"], tpr["macro"],
             label='ROC of fold {0} (area = {1:0.2f})'
             ''.format(i, auc),
             linestyle = ':',
             linewidth=3)
plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.title('Multi-class ROC of trans_mil')
# plt.show()
plt.savefig(fname = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Trans_mil.pdf')
plt.close()

####
## Attention_mil: 0.648
## Clam_sb:0.658
## Clam_mb: 0.668
## Trans_mil: 0.626
####

#### Patch特征和注意力分值 ####
path = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil'
files = [file for file in os.listdir(path)]

random.seed(42)
def sample_dict_with_indices(data, n_samples=200, seed=None):
    if seed is not None:
        random.seed(seed)

    sampled_data = {}
    sampled_indices = {}

    for key, matrix in data.items():
        if len(matrix) < n_samples:
            raise ValueError(f"Matrix '{key}' has fewer than {n_samples} rows.")

        all_indices = list(range(len(matrix)))
        selected_indices = random.sample(all_indices, n_samples)
        sampled_rows = [matrix[i] for i in selected_indices]

        sampled_data[key] = sampled_rows
        sampled_indices[key] = selected_indices

    return sampled_data, sampled_indices

files = sorted([file for file in files if re.search('clam_mb', file)])
config = mil_config('clam_mb', lr=1e-4)

preds = []
atts = []
features = []

for i in range(len(files)):

    file = files[i]
    model_path = os.path.join(path, file)
    pred, att = predict_mil(
        model = model_path,
        config = config,
        bags = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_ctranspath',
        outcomes = 'category',
        dataset = dataset,
        attention = True
    )
    feature = generate_mil_features(
        weights = model_path,
        dataset = dataset,
        bags = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/pt_files_ctranspath',
        config = config
    )
    # 修改文件以返回patch特征：
    # /home/ljc/anaconda3/envs/slideflow/lib/python3.9/site-packages/slideflow/mil/feature.py、
    # /home/ljc/anaconda3/envs/slideflow/lib/python3.9/site-packages/slideflow_gpl/clam/model.py

    att = [pd.DataFrame(x).T.to_dict(orient="records") for x in att]
    pred = pred.to_dict(orient="records")
    feature = {key: value.tolist() for key, value in feature.patch_activations.items()}
    feature, indice = sample_dict_with_indices(feature, n_samples=200) # 每个样本随机抽取200个patch
    att = [[att[index][:][i] for i in value] for index, (key, value) in enumerate(indice.items())] # 根据patch特征的索引抽取注意力分值
    feature = [pd.DataFrame(value).to_dict(orient="records") for key, value in feature.items()]

    preds.append(pred)
    atts.append(att)
    features.append(feature)

with open('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Clam_mb_prediction.json', 'w') as f:
    json.dump(preds, f)
with open('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Clam_mb_attention.json', 'w') as f:
    json.dump(atts, f)
with open('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Clam_mb_feature.json', 'w') as f:
    json.dump(features, f)


#### 空转样本推断 ####
from slideflow.mil import predict_slide

slide = '/home/ljc/DataDisk/0_ESCC/GDPH/1909933-6.svs'
wsi = sf.WSI(slide, tile_px=299, tile_um=302)
wsi.qc(qc)
model = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/00006-clam_mb-category'
y_pred, y_att = predict_slide(model, wsi, attention=True)


#### 获得空转样本特征 ####
# project_val = sf.create_project(
#     root='gdph',
#     slides='/home/ljc/DataDisk/0_ESCC/GDPH'
# ) # 创建项目
# project_val.extract_tiles(qc = qc,
#                       tile_px = 299,
#                       tile_um = 302,
#                       normalizer ='reinhard_fast',
#                       num_threads = 20
# ) # qc, 提取patch，使用Reinhard-fast归一化
project_val = sf.load_project(root='gdph') # 导入项目
dataset_val = project_val.dataset(tile_px=299, tile_um=302) # 创建数据集

path = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil'
files = [file for file in os.listdir(path)]
files = sorted([file for file in files if re.search('clam_mb', file)])
config = mil_config('clam_mb', lr=1e-4)

features_val = []

for i in range(len(files)):

    file = files[i]
    model_path = os.path.join(path, file)

    feature = generate_mil_features(
        weights = model_path,
        dataset = dataset_val,
        bags = '/home/ljc/0_project/0_ESCC/0_slideflow/gdph/pt_files_ctranspath',
        config = config
    )

    feature = {key: value.tolist() for key, value in feature.patch_activations.items()}
    feature = [pd.DataFrame(value).to_dict(orient="records") for key, value in feature.items()]

    features_val.append(feature)

with open('/home/ljc/0_project/0_ESCC/0_slideflow/gdph/Clam_mb_feature.json', 'w') as f:
    json.dump(features_val, f)

#### 获得空转样本的热图 ####
from matplotlib.colors import ListedColormap
import openslide

# 自己改写热图函数
def get_wsi(slide_path, level):

    image = openslide.open_slide(slide_path + ".svs")
    dimensions = list(image.dimensions)
    img = np.array(image.read_region(location=(0, 0), level=level, size=image.level_dimensions[level]))

    return dimensions, img

def get_heatmap(dimensions, coords, feature_value, rate, colors, size):

    # 降采样
    heatmap_data = np.full((int(dimensions[1] / rate), int(dimensions[0] / rate)), np.nan)
    for i in range(len(coords)):
        x, y = [int(t / rate) for t in coords[i]]
        heatmap_data[y:(y + int(size / rate)), x:(x + int(size / rate))] = feature_value.iloc[i]
    cmap = ListedColormap(colors)

    # Plot heatmap
    fig, ax = plt.subplots()
    im = ax.imshow(heatmap_data, cmap=cmap)
    cbar = ax.figure.colorbar(im, ax=ax)
    plt.xticks([])
    plt.yticks([])
    plt.show()

locations = np.load('/home/ljc/0_project/0_ESCC/0_slideflow/gdph/pt_files_ctranspath/1930377-13.index.npz')['arr_0']
values_df = pd.read_csv('/home/ljc/0_project/0_ESCC/0_slideflow/gdph/Clam_mb_1930377-13.csv')
values = values_df['predicted.d'].to_numpy()

dimensions, img = get_wsi(slide_path = '/home/ljc/DataDisk/0_ESCC/GDPH/1930377-13', level = 1)
get_heatmap(dimensions = dimensions,
            coords = locations,
            feature_value = values_df['predicted.b'],
            rate = 4,
            colors = ['dodgerblue', 'tomato', 'mediumseagreen', 'darkorange', 'rebeccapurple', 'gold'],
            size = 1144)