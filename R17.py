# %%
######不做其余处理，cut=4#########
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler


df = pd.read_csv("train10.csv",index_col=None)
df = df.drop(columns=["Unnamed: 0","depmap_id","drugname1","drugname2","drug1","drug2","cell line","score","seneitive10","drug_combination","cancer"])
#R1之后做特征集筛选删除mut数据#
df = df.loc[:, ~df.columns.str.startswith('mut_')]
X = df.drop(columns=["cut4"])
y = df["cut4"]

#归一化#
from sklearn.preprocessing import MinMaxScaler
transfer=MinMaxScaler()
X=transfer.fit_transform(X)


#数据分割与建模#
from imblearn.under_sampling import RandomUnderSampler
from collections import Counter
from sklearn.model_selection import train_test_split 
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve 
from sklearn.metrics import precision_recall_curve, auc
import matplotlib.pyplot as plt 
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import Lasso
from xgboost import XGBClassifier
from sklearn.svm import SVC
from joblib import parallel_backend

print('随机分类情况：{}'.format(Counter(y)))
X_train,X_test,y_train,y_test = train_test_split(X, y, test_size=0.2, random_state=22)

# %% [markdown]
# logistic

# %%

from joblib import parallel_backend
class_weight = {0: len(y_train) / len(y_train[y_train==0]), 1: 1.0}
#Logistic#
param_grid = {
    'C': [0.1, 1, 10],
    'penalty': ['l1', 'l2'],
    'solver': ['liblinear']
}

logreg = LogisticRegression(class_weight=class_weight)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
grid_search = GridSearchCV(logreg, param_grid=param_grid, cv=skf)
grid_search.fit(X_train, y_train)
print(f"Best accuracy: {grid_search.best_score_}")
print(f"Best parameters: {grid_search.best_params_}")

best_params = grid_search.best_params_
logreg = LogisticRegression(**best_params)
logreg.fit(X_train, y_train.astype('int'))

LRpredictions = logreg.predict(X_test)
print("Accuracy:",accuracy_score(y_test,LRpredictions))
print("Precision:",precision_score(y_test,LRpredictions))
print("Recall:",recall_score(y_test,LRpredictions))
print("F1 Score:",f1_score(y_test,LRpredictions))
print("混淆矩阵:",confusion_matrix(y_test,LRpredictions))

y_score = logreg.fit(X_train,y_train).predict_proba(X_test)
LR_fpr,LR_tpr,thresholds=roc_curve(y_test, y_score[:,1])
LR_auc =auc(LR_fpr,LR_tpr)
print("Logistic Regression AUC:",LR_auc)
plt.plot(LR_fpr,LR_tpr)
plt.title("Logistic Regression: ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate(recall)")
plt.savefig("logistic_auc17.png")  # Save the plot as PNG image
plt.close()

LR_precision, LR_recall, thresholds = precision_recall_curve(y_test, y_score[:,1])
LR_aupr = auc(LR_recall, LR_precision)
print("Logistic Regression AUPR:", LR_aupr)
plt.plot(LR_recall, LR_precision)
plt.title('Precision-Recall curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.savefig("logistic_aupr17.png")  # Save the plot as PNG image
plt.close()

# %%
import joblib
joblib.dump(logreg, "models/logistic17.pkl")
# 调用模型
# logreg = joblib.load("/models/logistic1.pkl")
# print(logreg)


# %% [markdown]
# lasso

# %%
from joblib import parallel_backend
#lasso#
param_grid = {
    'C': [0.1, 1, 10],
    'penalty': ['l1'],
    'solver': ['liblinear']
}

lasso = LogisticRegression(class_weight=class_weight)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
grid_search = GridSearchCV(lasso, param_grid=param_grid, cv=skf)
grid_search.fit(X_train, y_train)
print(f"Best accuracy: {grid_search.best_score_}")
print(f"Best parameters: {grid_search.best_params_}")

best_params = grid_search.best_params_
lasso = LogisticRegression(**best_params)
lasso.fit(X_train, y_train.astype('int'))

LApredictions = lasso.predict(X_test)
print("Accuracy:",accuracy_score(y_test,LApredictions))
print("Precision:",precision_score(y_test,LApredictions))
print("Recall:",recall_score(y_test,LApredictions))
print("F1 Score:",f1_score(y_test,LApredictions))
print("混淆矩阵:",confusion_matrix(y_test,LApredictions))

y_score = lasso.fit(X_train,y_train).predict_proba(X_test)
LA_fpr,LA_tpr,LA_thresholds=roc_curve(y_test, y_score[:,1])
LA_auc =auc(LA_fpr,LA_tpr)
print("Lasso Regression AUC:",LA_auc)
plt.plot(LA_fpr,LA_tpr)
plt.title("Lasso Regression: ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate(recall)")
plt.savefig("lasso_auc17.png")  # Save the plot as PNG image
plt.close()

LA_precision, LA_recall, LA_thresholds = precision_recall_curve(y_test, y_score[:,1])
LA_aupr = auc(LA_recall, LA_precision)
print("Lasso AUPR:", LA_aupr)
plt.plot(LA_recall, LA_precision)
plt.title('Precision-Recall curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.savefig("lasso_aupr17.png")  # Save the plot as PNG image
plt.close()

# %%
import joblib
joblib.dump(lasso, "models/lasso17.pkl")


# %% [markdown]
# RF

# %%
param_grid = {
    'n_estimators': [10,20,30,50,70,100],
    'max_depth': [5,7,9],
}

rfc = RandomForestClassifier(class_weight=class_weight)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42) #不平衡处理
grid_search = GridSearchCV(rfc, param_grid=param_grid, cv=skf)
grid_search.fit(X_train, y_train)
print(f"Best accuracy: {grid_search.best_score_}")
print(f"Best parameters: {grid_search.best_params_}")

best_params = grid_search.best_params_
rfc = RandomForestClassifier(**best_params)
rfc.fit(X_train, y_train.astype('int'))

RFpredictions = rfc.predict(X_test)
print("Accuracy:",accuracy_score(y_test,RFpredictions))
print("Precision:",precision_score(y_test,RFpredictions))
print("Recall:",recall_score(y_test,RFpredictions))
print("F1 Score:",f1_score(y_test,RFpredictions))
print("混淆矩阵:",confusion_matrix(y_test,RFpredictions))

y_score = rfc.fit(X_train,y_train).predict_proba(X_test)
RF_fpr,RF_tpr,RF_thresholds=roc_curve(y_test, y_score[:,1])
RF_auc =auc(RF_fpr,RF_tpr)
print("Random Forest AUC:",RF_auc)
plt.plot(RF_fpr,RF_tpr)
plt.title("Random Forest: ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate(recall)")
plt.savefig("RF_auc17.png")  # Save the plot as PNG image
plt.close()
RF_precision, RF_recall, RF_thresholds = precision_recall_curve(y_test, y_score[:,1])
RF_aupr = auc(RF_recall, RF_precision)
print("Random Forest AUPR:", RF_aupr)
plt.plot(RF_recall, RF_precision)
plt.title('Precision-Recall curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.savefig("rf_aupr17.png")  # Save the plot as PNG image
plt.close()

# %%
import joblib
joblib.dump(rfc, "models/rfc17.pkl")



# %% [markdown]
# svm

# %%
param_grid = {

     'C': [5],
     'kernel': ['poly'],
     'gamma': ['scale'],

 }

class_weight = {0: len(y_train) / len(y_train[y_train==0]), 1: 1.0}
svm = SVC(class_weight=class_weight)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
grid_search = GridSearchCV(svm, param_grid=param_grid, cv=skf)
grid_search.fit(X_train, y_train)
print(f"Best accuracy: {grid_search.best_score_}")
print(f"Best parameters: {grid_search.best_params_}")

best_params = grid_search.best_params_
svm = SVC(**best_params)
svm.fit(X_train, y_train.astype('int'))

SVMpredictions = svm.predict(X_test)
print("Accuracy:",accuracy_score(y_test,SVMpredictions))
print("Precision:",precision_score(y_test,SVMpredictions))
print("Recall:",recall_score(y_test,SVMpredictions))
print("F1 Score:",f1_score(y_test,SVMpredictions))
print("混淆矩阵:",confusion_matrix(y_test,SVMpredictions))

y_score = svm.fit(X_train,y_train).decision_function(X_test)
plt.figure(figsize=(5,5))
SVM_fpr,SVM_tpr,SVM_thresholds=roc_curve(y_test, y_score)
SVM_auc =auc(SVM_fpr,SVM_tpr)
print("SVM AUC:",SVM_auc)
plt.plot(SVM_fpr,SVM_tpr)
plt.title("SVM: ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate(recall)")
plt.savefig("svm_auc17.png")  # Save the plot as PNG image
plt.close()

SVM_precision, SVM_recall, SVM_thresholds = precision_recall_curve(y_test, y_score)
SVM_aupr = auc(SVM_recall, SVM_precision)
print("SVM AUPR:", SVM_aupr)
plt.plot(SVM_recall, SVM_precision)
plt.title('Precision-Recall curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.savefig("svm_aupr17.png")  # Save the plot as PNG image
plt.close()

# %%
import joblib
joblib.dump(svm, "models/svm17.pkl")

# %% [markdown]
# xgboost

# %%

param_grid = {
    'learning_rate': [0.05],
    'max_depth': [11],
    'n_estimators': [300,500,700]
}
scale_pos_weight = len(y_train[y_train==0]) / len(y_train[y_train==1])
xgb = XGBClassifier(scale_pos_weight=scale_pos_weight)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
grid_search = GridSearchCV(xgb, param_grid=param_grid, cv=skf)
grid_search.fit(X_train, y_train)

print(f"Best accuracy: {grid_search.best_score_}")
print(f"Best parameters: {grid_search.best_params_}")

best_params = grid_search.best_params_
xgb = XGBClassifier(**best_params)
xgb.fit(X_train, y_train.astype('int'))
XBpredictions = xgb.predict(X_test)

print("Accuracy:",accuracy_score(y_test,XBpredictions))
print("Precision:",precision_score(y_test,XBpredictions))
print("Recall:",recall_score(y_test,XBpredictions))
print("F1 Score:",f1_score(y_test,XBpredictions))
print("混淆矩阵:",confusion_matrix(y_test,XBpredictions))

y_score = xgb.fit(X_train,y_train).predict_proba(X_test)
plt.figure(figsize=(5,5))
XB_fpr,XB_tpr,thresholds=roc_curve(y_test, y_score[:,1])
XB_auc =auc(XB_fpr,XB_tpr)
print("XGbooost_auc:",XB_auc)
plt.plot(XB_fpr,XB_tpr)
plt.title("XGbooost:Roc curve")
plt.xlabel("FPR")
plt.ylabel("TPR(recall)")
plt.savefig("xgboost_auc17.png")  # Save the plot as PNG image
plt.close()

XB_precision, XB_recall, XB_thresholds = precision_recall_curve(y_test, y_score[:,1])
XB_aupr = auc(XB_recall, XB_precision)
print("XGboost AUPR:", XB_aupr)
plt.plot(XB_recall, XB_precision)
plt.title('Precision-Recall curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.savefig("xgboost_aupr17.png")  # Save the plot as PNG image
plt.close()

# %%
import joblib
joblib.dump(xgb, "models/xgboost17.pkl")

# %%
plt.plot(LR_fpr, LR_tpr, color='red', lw=2, label='Logistic Regression (AUC=%0.4f)' % LR_auc)
plt.plot(LA_fpr, LA_tpr, color='green', lw=2, label='Lasso Regression (AUC=%0.4f)' % LA_auc)
plt.plot(RF_fpr, RF_tpr, color='purple', lw=2, label='Random Forest (AUC=%0.4f)' % RF_auc)
plt.plot(SVM_fpr, SVM_tpr, color='orange', lw=2, label='SVM (AUC=%0.4f)' % SVM_auc)
plt.plot(XB_fpr, XB_tpr, color='yellow', lw=2, label='XGBOOST (AUC=%0.4f)' % XB_auc)

plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig("auc17.png")  # Save the plot as PNG image
plt.close()

# %%
plt.plot(LR_recall, LR_precision, color='red', lw=2, label='Logistic Regression (AUPR=%0.4f)' % LR_aupr)
plt.plot(LA_recall, LA_precision, color='green', lw=2, label='Lasso Regression (AUPR=%0.4f)' % LA_aupr)
plt.plot(RF_recall, RF_precision, color='purple', lw=2, label='Random Forest (AUPR=%0.4f)' % RF_aupr)
plt.plot(SVM_recall, SVM_precision, color='orange', lw=2, label='SVM (AUPR=%0.4f)' % SVM_aupr)
plt.plot(XB_recall, XB_precision, color='yellow', lw=2, label='XGBOOST (AUPR=%0.4f)' % XB_aupr)

#plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall curve')
plt.legend(loc="lower right")
plt.savefig("aupr17.png")  # Save the plot as PNG image
plt.close()


