# %% [markdown]
# 更改阈值为4

# %%
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, roc_curve, auc, precision_recall_curve
from itertools import combinations
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from imblearn.under_sampling import RandomUnderSampler
from collections import Counter
from xgboost import XGBClassifier
df = pd.read_csv("train10.csv", index_col=None)
df = df.drop(columns=["depmap_id", "drugname1", "drugname2", "drug1", "drug2", "cell line", "score", "seneitive10",
                      "drug_combination", "cancer"])

# 特征集列表
feature_sets = [
    df.filter(regex='^DD_maccs_'),
    df.filter(regex='^DD_tox_'),
    df.filter(regex='^DDP_'),
    df.filter(regex='^exp_'),
    df.filter(regex='^cnv_'),
    df.filter(regex='^methy_'),
    df.filter(regex='^mut_'),
    df.filter(regex='^RNAi_')
]

# 归一化
def normalize_features(features):
    transfer = MinMaxScaler()
    return transfer.fit_transform(features)

# 获取所有特征集的组合
all_combinations = []
for r in range(1, len(feature_sets) + 1):
    all_combinations.extend(combinations(feature_sets, r))

# 创建一个空的DataFrame来存储结果
result_df = pd.DataFrame(columns=["Feature Combination (Regex)", "Accuracy", "Precision", "Recall", "F1 Score", "AUC", "AUPR"])

# 遍历不同的特征组合
for feature_combination in all_combinations:
    # 合并特征
    X = pd.concat(list(feature_combination), axis=1)
    y = df["cut4"]

    # 归一化
    X = normalize_features(X)

    # 数据分割与建模
    print('随机分类情况：{}'.format(Counter(y)))
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=22)

    # 随机森林超参数调整
    param_grid = {
    'learning_rate': [0.1],
    'max_depth': [7],
    'n_estimators': [700]
}

    scale_pos_weight = len(y_train[y_train==0]) / len(y_train[y_train==1])
    xgb = XGBClassifier(scale_pos_weight=scale_pos_weight)

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    grid_search = GridSearchCV(xgb, param_grid=param_grid, cv=skf)
    grid_search.fit(X_train, y_train)
    
    best_params = grid_search.best_params_
    xgb = XGBClassifier(**best_params)
    xgb.fit(X_train, y_train.astype('int'))
    XBpredictions = xgb.predict(X_test)

    # 获取 regex 组合
    regex_combination = [feature_set.columns[0] for feature_set in feature_combination if feature_set.shape[1] > 0]
    accuracy = accuracy_score(y_test, xgb.predict(X_test))
    precision = precision_score(y_test, xgb.predict(X_test))
    recall = recall_score(y_test, xgb.predict(X_test))
    f1 = f1_score(y_test, xgb.predict(X_test))

    y_score = xgb.fit(X_train,y_train).predict_proba(X_test)
    xgb_fpr,xgb_tpr,xgb_thresholds=roc_curve(y_test, y_score[:,1])
    xgb_auc =auc(xgb_fpr,xgb_tpr)

    xgb_precision, xgb_recall, xgb_thresholds = precision_recall_curve(y_test, y_score[:,1])
    xgb_aupr = auc(xgb_recall, xgb_precision)

    # 输出完成消息
    print(f"Feature Combination (Regex) completed: {regex_combination}")
    print(f"Accuracy: {accuracy}")
    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"F1 Score: {f1}")
    print(f"AUC: {xgb_auc}")
    print(f"AUPR: {xgb_aupr}")
    print("\n")
    # 添加到DataFrame
    result_df = pd.concat([result_df, pd.DataFrame({
        "Feature Combination (Regex)": [regex_combination],
        "Accuracy": [accuracy],
        "Precision": [precision],
        "Recall": [recall],
        "F1 Score": [f1],
        "AUC": [xgb_auc],
        "AUPR": [xgb_aupr]
    })], ignore_index=True)

# 将结果保存到CSV文件
result_df.to_csv('model_results.csv', index=False)



