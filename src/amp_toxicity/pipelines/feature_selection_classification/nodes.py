"""
This is a boilerplate pipeline 'feature_selection_classification'
generated using Kedro 0.18.7

Created on Tue Apr 18 14:01:18 2023

@author: JJvRijn
"""
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import  matthews_corrcoef

def Classification (Data_complete):
    #happenn formating
    Data_complete["toxic"] = "~"
    #temp fix
    for i in range(len(Data_complete)):
        if Data_complete["hem_activity"][i]==None:
            Data_complete =Data_complete.drop(index = [i] , axis=0)
    Data_complete =Data_complete.reset_index()
    for v in range(len(Data_complete)):
        if float(Data_complete["hem_activity"][v]) <= 5:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 30:
                Data_complete.loc[v, "toxic"] = 0        
        if 5 <= float(Data_complete["hem_activity"][v]) <= 10:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 60:
                Data_complete.loc[v, "toxic"] = 0   
        if 10 <= float(Data_complete["hem_activity"][v]) <= 15:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 90:
                Data_complete.loc[v, "toxic"] = 0
        if 15 <= float(Data_complete["hem_activity"][v]) <= 20:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 120:
                Data_complete.loc[v, "toxic"] = 0
        if 20 <= float(Data_complete["hem_activity"][v]) <= 25:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 180:
                Data_complete.loc[v, "toxic"] = 0  
        if 25 <= float(Data_complete["hem_activity"][v]) <= 30:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 180:
                Data_complete.loc[v, "toxic"] = 0  
        if 30 <= float(Data_complete["hem_activity"][v]) <= 35:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 210:
                Data_complete.loc[v, "toxic"] = 0         
        if 35 <= float(Data_complete["hem_activity"][v]) <= 40:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 240:
                Data_complete.loc[v, "toxic"] = 0  
        if 40 <= float(Data_complete["hem_activity"][v]) <= 45:
            if (Data_complete["CONCENTRATION_µM"][v]) >= 270:
                Data_complete.loc[v, "toxic"] = 0  
        if float(Data_complete["hem_activity"][v]) >= 50:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 300:
                Data_complete.loc[v, "toxic"] = 1
        if float(Data_complete["hem_activity"][v]) >= 55:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 330:
                Data_complete.loc[v, "toxic"] = 1
        if float(Data_complete["hem_activity"][v]) >= 60:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 360:
                Data_complete.loc[v, "toxic"] = 1
        if float(Data_complete["hem_activity"][v]) >= 65:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 390:
                Data_complete.loc[v, "toxic"] = 1 
        if float(Data_complete["hem_activity"][v]) >= 70:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 420:
                Data_complete.loc[v, "toxic"] = 1 
        if float(Data_complete["hem_activity"][v]) >= 75:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 450:
                Data_complete.loc[v, "toxic"] = 1 
        if float(Data_complete["hem_activity"][v]) >= 80:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 480:
                Data_complete.loc[v, "toxic"] = 1 
        if float(Data_complete["hem_activity"][v]) >= 85:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 510:
                Data_complete.loc[v, "toxic"] = 1 
        if float(Data_complete["hem_activity"][v]) >= 90:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 540:
                Data_complete.loc[v, "toxic"] = 1
        if float(Data_complete["hem_activity"][v]) >= 95:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 570:
                Data_complete.loc[v, "toxic"] = 1
        if float(Data_complete["hem_activity"][v]) >= 100:
            if (Data_complete["CONCENTRATION_µM"][v]) <= 600:
                Data_complete.loc[v, "toxic"] = 1 
    return Data_complete[["seq", "toxic"]]

def Remove_unclear_tox(Tox_label, remove_draws = True):
    #remove rows were toxicity is unclear
    Tox_label = Tox_label[Tox_label["toxic"] != '~']
    Tox_label = Tox_label.reset_index(drop=True)
    #merge same peptides
    #where majority vote decides toxicity and a draw gives "~"
    peptides = Tox_label["seq"].unique()
    Tox_label_final = pd.DataFrame(peptides, columns=['seq'])
    for a in range(len(Tox_label)):
        for n in range(len(peptides)):
            toxic = []
            toxic_count = 0
            if Tox_label["seq"][a] == peptides[n]:
                toxic.append(Tox_label["toxic"][a])
                toxic_count = sum(toxic)/len(toxic)
                if toxic_count == len(toxic)/2:
                    Tox_label_final.loc[n, "toxic"] = "~"
                if toxic_count < len(toxic)/2:
                    Tox_label_final.loc[n, "toxic"] = 0
                if toxic_count > len(toxic)/2:
                    Tox_label_final.loc[n, "toxic"] = 1
    if remove_draws is True:
        Tox_label_final = Tox_label_final[Tox_label_final["toxic"] != '~']
    return Tox_label_final[["seq", "toxic"]]

def sort_best_features(all_features, seq_and_label):
    all_features = all_features.loc[all_features['seq'].isin(seq_and_label.iloc[:,0])].reset_index(drop = True)
    all_features = all_features.sort_values(by=["seq"])
    seq_and_label = seq_and_label.sort_values(by=["seq"])
    no_seq = all_features.drop(columns = ["seq"])
    if list(all_features["seq"]) == list(seq_and_label.iloc[:,0]):
        fs = SelectKBest(score_func=f_classif, k="all")
        fs.fit(all_features[no_seq.columns], seq_and_label.iloc[:,1],)
        ordered_feat = fs.get_feature_names_out(input_features=no_seq.columns)
    else:
        print("features and labels dont match")
    return ordered_feat
    

def split_data_Kmeans (all_features, seq_and_label, percent_test=0.1):
    all_features = all_features.loc[all_features['seq'].isin(seq_and_label.iloc[:,0])].reset_index(drop = True)
    all_features = all_features.sort_values(by=["seq"])
    seq_and_label = seq_and_label.sort_values(by=["seq"])
    no_seq = all_features.drop(columns = ["seq"])
    if list(all_features["seq"]) == list(seq_and_label.iloc[:,0]):
        pos_samples = all_features[seq_and_label.iloc[:,1]==1]
        pos_samples = pos_samples.reset_index(drop=True)
        neg_samples = all_features[seq_and_label.iloc[:,1]==0]
        neg_samples = neg_samples.reset_index(drop=True)
        #clustering the for both toxic and non-toxic and saving the labels
        number_test = round(percent_test * len(all_features) *0.5)
        
        Kmeans_cluster_pos = KMeans(n_clusters = number_test).fit(pos_samples[no_seq.columns])
        pos_labels = Kmeans_cluster_pos.labels_
        Kmeans_cluster_neg = KMeans(n_clusters = number_test).fit(neg_samples[no_seq.columns])
        neg_labels = Kmeans_cluster_neg.labels_
        #grabbing the first peptide of a cluster for the test set, left overs in to training set
        pos_test = pos_samples[pos_labels == 1].iloc[:1]
        neg_test = neg_samples[neg_labels == 1].iloc[:1]
        pos_train = pos_samples[pos_labels == 1].iloc[1:]
        neg_train = neg_samples[neg_labels == 1].iloc[1:]
        #merging the toxic and non-toxic parts of the sets
        KMean_train = pd.concat([pos_train, neg_train], axis = 0)
        KMean_test = pd.concat([pos_test, neg_test], axis = 0)
        for j in range(1, number_test):
            pos_test = pos_samples[pos_labels == j].iloc[:1]
            neg_test = neg_samples[neg_labels == j].iloc[:1]
            pos_train = pos_samples[pos_labels == j].iloc[1:]
            neg_train = neg_samples[neg_labels == j].iloc[1:]
            KMean_train = pd.concat([KMean_train, pos_train, neg_train], axis = 0)
            KMean_test = pd.concat([KMean_test, pos_test, neg_test], axis = 0)
   
        KMean_test_feat = KMean_test.drop_duplicates(subset=["seq"]).reset_index(drop=True).sort_values(by=["seq"])
        KMean_train_feat = KMean_train.drop_duplicates(subset=["seq"]).reset_index(drop=True).sort_values(by=["seq"])
        Kmean_test_labels = seq_and_label[seq_and_label['seq'].isin(KMean_test_feat['seq'])].sort_values(by=["seq"])
        Kmean_train_labels = seq_and_label[seq_and_label['seq'].isin(KMean_train_feat['seq'])].sort_values(by=["seq"])
        return KMean_train_feat, KMean_test_feat, Kmean_train_labels, Kmean_test_labels
    else:
        return print("features and labels dont match")
    

def test_n_features (KMean_train_feat, Kmean_train_labels, top_labels, n_feat=[10, 20, 50, 100, 110, 120, 130, 140, 150 , 160, 170, 180, 190 ,200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 5000]):
    mcc = []    
    for i in n_feat:
        used_feat_labels = top_labels[:i]
        used_feat = KMean_train_feat[used_feat_labels]
        train_X, val_X, train_y, val_y = train_test_split(used_feat ,Kmean_train_labels.iloc[:,1], random_state = 445)
        RFC = RandomForestClassifier()
        RFC.fit(train_X, train_y)
        #results before optimization of parameters
        RF_pred = RFC.predict(val_X)
        mcc.append(matthews_corrcoef(val_y, RF_pred))
    plt.plot(n_feat, mcc, marker='.')
    plt.xlabel('Number of features')
    plt.ylabel('MCC')
    return plt