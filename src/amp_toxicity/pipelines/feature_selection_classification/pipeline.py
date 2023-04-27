"""
This is a boilerplate pipeline 'feature_selection_classification'
generated using Kedro 0.18.7

Created on Tue Apr 18 14:01:18 2023

@author: JJvRijn
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import Classification, Remove_unclear_tox, sort_best_features, split_data_Kmeans, test_n_features
import pandas as pd

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
        node(
            func=Classification,
            inputs="data_complete",
            outputs="Labeled_data",
            name="Classify_the_data",
        ),node(
            func=Remove_unclear_tox,
            inputs="Labeled_data",
            outputs="Labeled_data_v2",
            name="Remove_unclear_labels",
        ),
        node(
            func=lambda *dfs: pd.concat(dfs, axis=1).T.drop_duplicates().T,
            inputs=["AAC_features", "DAAC_features", 
                    "MAAC_features" ,"Compact_fingerprint_features"],
            outputs="combined_features",
            name="combine_feat",
        ),
        node(
            func=sort_best_features,
            inputs=["combined_features", "Labeled_data_v2"],
            outputs="ranked_features",
            name="Rank_the_features",
        ),
        node(
            func=split_data_Kmeans,
            inputs=["combined_features","Labeled_data_v2"],
            outputs=["KMean_train_feat", "KMean_test_feat",
                     "Kmean_train_labels", "Kmean_test_labels"],
            name="Split_data_using_Kmeans",
        ),
        node(
            func=test_n_features,
            inputs=["KMean_train_feat", "Kmean_train_labels", "ranked_features"],
            outputs="n_feat_graph",
            name="plot_graph_n_feat",
        )
    ]
)
