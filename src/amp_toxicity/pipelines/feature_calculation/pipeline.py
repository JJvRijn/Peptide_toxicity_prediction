"""
This is a boilerplate pipeline 'feature_calculation'
generated using Kedro 0.18.7
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import AAC_df, DAAC_df, MAAC_df, compact_fp_df, BLOSUM62_df, distrubution_df, unique_peptides

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=unique_peptides,
                inputs="data_complete",
                outputs="unique_peptides",
                name="unique_peptides_from_data",
            ),node(
                func=AAC_df,
                inputs="unique_peptides",
                outputs="AAC_features",
                name="AminoAcidComposition",
            ),
            node(
                func=DAAC_df,
                inputs="unique_peptides",
                outputs="DAAC_features",
                name="DiAminoAcidComposition",
            ),
            node(
                func=MAAC_df,
                inputs="unique_peptides",
                outputs="MAAC_features",
                name="MultiAminoAcidComposition",
            ),
            node(
                func=compact_fp_df,
                inputs="unique_peptides",
                outputs="Compact_fingerprint_features",
                name="Compact_fingerprint",
            ),
            node(
                func=BLOSUM62_df,
                inputs="unique_peptides",
                outputs="BLOSUM_matrix",
                name="BLOSUM_matrix",
            ),
            node(func=distrubution_df,
                 inputs="unique_peptides",
                 outputs="AA_distribution_features",
                 name="AA_distribution",                 
                )
        ]
    )