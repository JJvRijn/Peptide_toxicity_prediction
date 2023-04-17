# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 17:41:30 2023

@author: JJvRijn
"""

from kedro.pipeline import Pipeline, node, pipeline

from .nodes import Scraping_DBAASP, Scraping_DRAMP, string_to_act_conc, concat_w_source


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
            node(
                func= Scraping_DBAASP,
                inputs="DBAASP_complete",
                outputs="DBAASP_processed",
                name="Scraping_DBAASP_node",
            ),
            node(
                func=Scraping_DRAMP,
                inputs="DRAMP_raw",
                outputs="DRAMP_processed",
                name="Scraping_DRAMP_node",
            ),
             node(
                func=string_to_act_conc,
                inputs="hemolytik_raw_data",
                outputs="hemolytik_raw_processed",
                name="Hemolytik_string_to_conc_node",
            ),
            node(
                func=concat_w_source,
                inputs=["DBAASP_processed", "DRAMP_processed", "hemolytik_raw_processed", "DADP_data"],
                outputs="data_complete",
                name="concat_with_source",
            )
        ])
