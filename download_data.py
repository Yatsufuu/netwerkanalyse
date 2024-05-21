"""This is the python code to perform the entire check for
unsewered plots (ongerioleerde percelen).
"""

import os
import logging
import cProfile
import pstats
import geopandas as gpd
import pandas as pd
from pathlib import Path

import sys

# Internal packages
from dataportaal.calc.select_plots import clip_buffer_select_part
from dataportaal.calc.select_plots import select_items_with_spatial_index
import dataportaal.calc.utilities as ut
import dataportaal.io.pdok_wfs as pdok_wfs
import dataportaal.io.pdok_api_kadastralekaart as api_kad
import dataportaal.io.gwsw_wfs as gwsw_wfs
import dataportaal.io.bgt_api as bgt_api
from settings import MUNICIPALITIES


for municipality in MUNICIPALITIES:
    print(f"Downloading and preparing data for '{municipality}'")
    data_folder = f"data/{municipality}"
    # create the data folder if it is not already done
    Path(f"./{data_folder}").mkdir(parents=True, exist_ok=True)

    # Start logging
    logging.basicConfig(
        filename=f"./{data_folder}/{municipality.lower()}_get_data.log",
        level=logging.INFO,
        filemode="w",
    )
    logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)

    profiler = cProfile.Profile()
    profiler.enable()

    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats("cumulative")
    stats.print_stats(20)

    gdf_municipality, polygon_municipality = pdok_wfs.get_contour_municipalities(
        municipality, data_folder, False
    )

    logging.info("Downloading new data for municipality %s", municipality)
    logging.info("Downloading plots")
    gdf_plots = api_kad.get_gdf_plots(
        polygon_municipality, data_folder, cache_data=True
    )
    logging.info("Downloading sewer lines")
    gdf_sewer = gwsw_wfs.get_gwsw_sewer_lines(
        municipality, data_folder, cache_data=True
    )
    logging.info("Downloading pump points")
    gdf_pump = gwsw_wfs.get_gwsw_pump_points(municipality, data_folder, cache_data=True)
    logging.info("Downloading BGT roads and water")
    bgt_dict = bgt_api.get_bgt_features_poly(
        polygon_municipality.wkt, data_folder, cache_data=True
    )
    gdf_water = bgt_dict["bgt_waterdeel"]

    logging.info("Extracting railroads from BGT roads.")
    gdf_roads = bgt_dict["bgt_wegdeel"]
    gdf_railroads = gdf_roads.loc[gdf_roads["function"] == "spoorbaan"]

    # ===== Some extra logging ============
    if gdf_plots.empty:
        logging.warning("Empty Plots GeoDataFrame found")
    if gdf_sewer.empty:
        logging.warning("Empty Sewer GeoDataFrame found")
    if gdf_pump.empty:
        logging.warning("Empty Pump GeoDataFrame found")

    gdf_plots.to_file(str(Path(data_folder) / "01_plots.shp"))
    gdf_sewer.to_file(str(Path(data_folder) / "02_sewerlines.shp"))
    gdf_pump.to_file(str(Path(data_folder) / "03_pumps.shp"))
    gdf_water.to_file(str(Path(data_folder) / "04_water.shp"))
    gdf_roads.to_file(str(Path(data_folder) / "05_roads.shp"))
    gdf_railroads.to_file(str(Path(data_folder) / "06_railroads.shp"))
