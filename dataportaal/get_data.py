"""This is the python code to perform the entire check for
unsewered plots (ongerioleerde percelen).
"""

import os
import logging
import cProfile
import pstats
import geopandas as gpd
import pandas as pd

import sys

cur_dir = os.getcwd()
if not cur_dir in sys.path:
    sys.path.insert(0, cur_dir)

# Internal packages
from dataportaal.calc.select_plots import clip_buffer_select_part
from dataportaal.calc.select_plots import select_items_with_spatial_index
import dataportaal.calc.utilities as ut
import dataportaal.io.pdok_wfs as pdok_wfs
import dataportaal.io.pdok_api_kadastralekaart as api_kad
import dataportaal.io.gwsw_wfs as gwsw_wfs
import dataportaal.io.bgt_api as bgt_api

# Mogelijke gemeentes
# "Leiden",
# "Gouda",
# "Lisse",
# "Teylingen",
# "Hillegom",
# "Katwijk",
# "Wassenaar",
# "Leiderdorp",
# "Noordwijk",
# "Zoeterwoude",
# "Bodegraven-Reeuwijk",
# "Haarlemmermeer",
# "Alphen aan den Rijn",
MUNICIPALITY = "Noordwijk"

# Start logging
logging.basicConfig(
    filename=f"./{MUNICIPALITY.lower()}_get_data.log",
    level=logging.DEBUG,
    filemode="w",
)
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)


def run_opc_per_municipality(municipalities, data_folder, use_cache=True, test=True):
    """
    Runs the check for a given list of municipalities,
    reads and stores per municipality data in
    <workingdirectory>\\data\\<municipality name>
    Calls the function 'ongerioleerde_percelen_check' in which the
    calculations are performed.

    Parameters
    ----------
    municipalities: array of str
        The names of the municipalities to do the calculation for.
    data_folder: string
        Base path where the output should be written to.
    usecache = bool, optional
        If true checks if data already exist as a file in the data directory.
        Default = True
    test = bool, optional
        If True, only performs calculation for a small area and not for the
        given municipalities. Default = False

    Returns
    -------
    None
    """
    for municipality in municipalities:
        logging.info(
            "==============> Processing municipality %s <===========", municipality
        )
        ongerioleerde_percelen_check(
            municipality,
            os.path.join(data_folder, municipality),
            use_cache=use_cache,
            test=test,
        )
        logging.info(
            "==============> Finished %s successfully<===========", municipality
        )


def read_data(municipality, data_folder, use_cache, test):
    """
    Reads all data that are needed for the calculation from several web
    services or -if available- from local files.
    If use_cache is false, reads all data from the native source, otherwise
    the data will be read from the given data_folder.
    If downloaded, the data will be saved to '<data_folder>/<municipality>'.

    Parameters
    ----------
    municipality: str
        The name of the municipality
    data_folder: str
        The base output directory
    use_cache: bool
        If true, it will be tried to read data from disk. If in that
        case no data is found, it will be downloaded.
    test: bool
        If true a small area will be used to perform the calculation on.

    Returns
    -------
    gdf_plots: geodataframe
        Contains all plots for the municipality. (from BRK)
    gdf_sewer: geodataframe
        Contains the sewer lines. (from GWSW)
    gdf_pump: geodataframe
        Contains sewer pumps (from GWSW)
    gdf_water: geodataframe
        Contains water polygons (from BGT)
    gdf_railroads:
        Contains railroad polygons (from BGT)
    """
    # get the municipality contour
    gdf_municipality, polygon_municipality = pdok_wfs.get_contour_municipalities(
        municipality, data_folder, use_cache
    )

    if test:
        # Choose a small area within the polygon of the municipality to limit
        # the number of features for test purposes.
        (
            xmin_test,
            ymin_test,
            xmax_test,
            ymax_test,
            gdf_municipality,
            polygon_municipality,
        ) = ut.create_test_bbox(gdf_municipality)

    if use_cache:
        logging.info("Cache mode")
        try:
            logging.info("Loading plots GeoDataFrame from cache")
            gdf_plots = gpd.read_parquet(os.path.join(data_folder, "1_plots.parquet"))

        except IOError as error:
            logging.info("Loading from cache failed, file not found. %s", error)
            gdf_plots = api_kad.get_gdf_plots(
                polygon_municipality, data_folder, cache_data=True
            )

        try:
            logging.info("Loading sewer lines GeoDataFrame from cache")
            gdf_sewer = gpd.read_parquet(
                os.path.join(data_folder, "2a_sewer_lines.parquet")
            )

        except IOError as error:
            logging.info("Loading from cache failed, file not found, %s", error)
            gdf_sewer = gwsw_wfs.get_gwsw_sewer_lines(
                str.replace(municipality, " ", ""), data_folder, cache_data=True
            )
        try:
            logging.info("Loading pump points GeoDataFrame from cache")
            gdf_pump = gpd.read_parquet(
                os.path.join(data_folder, "2b_pump_points.parquet")
            )

        except IOError as error:
            logging.info("Loading from cache failed, file not found, %s", error)
            gdf_pump = gwsw_wfs.get_gwsw_pump_points(
                str.replace(municipality, " ", ""), data_folder, cache_data=True
            )

        try:
            logging.info("Loading BGT GeoDataFrame from cache")
            gdf_water = gpd.read_parquet(
                os.path.join(data_folder, "bgt_waterdeel.parquet")
            )
            gdf_roads = gpd.read_parquet(
                os.path.join(data_folder, "bgt_wegdeel.parquet")
            )
            gdf_railroads = gdf_roads.loc[gdf_roads["function"] == "spoorbaan"]
        except IOError as error:
            logging.info("Loading BGT from cache failed, file not found")
            logging.info("Downloading BGT roads and water, %s.", error)
            bgt_dict = bgt_api.get_bgt_features_poly(
                polygon_municipality.wkt, data_folder, cache_data=True
            )
            gdf_water = bgt_dict["bgt_waterdeel"]
            gdf_roads = bgt_dict["bgt_wegdeel"]
            # Attention: there are municipalities without railroads!
            gdf_railroads = gdf_roads.loc[gdf_roads["function"] == "spoorbaan"]
    else:
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
        gdf_pump = gwsw_wfs.get_gwsw_pump_points(
            municipality, data_folder, cache_data=True
        )
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
    # ===== Limit the data to the test bounding box if test is True =============
    if test:
        logging.info("Make test selection for development")
        gdf_sewer = gdf_sewer.to_crs(gdf_plots.crs)
        gdf_pump = gdf_pump.to_crs(gdf_plots.crs)

        # gdf_plots = gdf_plots.cx[xmin_test:xmax_test,
        #                          ymin_test:ymax_test]
        gdf_sewer = gdf_sewer.cx[xmin_test:xmax_test, ymin_test:ymax_test]
        gdf_pump = gdf_pump.cx[xmin_test:xmax_test, ymin_test:ymax_test]
        gdf_water = gdf_water.cx[xmin_test:xmax_test, ymin_test:ymax_test]
        gdf_railroads = gdf_railroads.cx[xmin_test:xmax_test, ymin_test:ymax_test]
    return gdf_plots, gdf_sewer, gdf_pump, gdf_water, gdf_railroads


def ongerioleerde_percelen_check(municipality, data_folder, use_cache=False, test=True):
    """
    Performs the check for a given municipality name. Downloaded data will be saved
    to a directory: <data_folder>/<municipality>

    Parameters
    ----------
    municipality: str
        The name of the municipality
    data_folder: str
        The base output directory
    use_cache: bool
        If true, it will be tried to read data from disk. If in that
        case no data is found, it will be downloaded.
    test: bool
        If true a small area will be set to perform the test on.

    """
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)

    # Step 1: ====== Read the data ====================
    gdf_plots, gdf_sewer, gdf_pump, gdf_water, gdf_railroads = read_data(
        municipality, data_folder, use_cache, test
    )

    # This step is added because otherwise the selection with
    # spatial indices goes wrong.
    # gdf_plots = gdf_plots.reset_index()
    # # create a folder to save the shape files
    result_folder = os.path.join(data_folder, "shapes")
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)

    # Save the input files as shape.
    gdf_plots.to_file(os.path.join(result_folder, "kadaster_plots.shp"))
    gdf_sewer.to_file(os.path.join(result_folder, "sewers.shp"))
    gdf_pump.to_file(os.path.join(result_folder, "pumps.shp"))
    gdf_water.to_file(os.path.join(result_folder, "bgt_water.shp"))
    if not gdf_railroads.empty:
        gdf_railroads.to_file(os.path.join(result_folder, "bgt_rail.shp"))

    # Step 2: =========== Create buffer for sewer and pumps ==========
    logging.info(
        "Create polygon of the 40m buffer %s", " around sewer lines and pump points"
    )
    gdf_pump_sewer = gpd.GeoDataFrame(pd.concat([gdf_pump, gdf_sewer])).reset_index()

    buffer_total = gdf_pump_sewer.buffer(40)

    # Step 3: ============ Create mask from the railroads and water ========
    logging.info("Create mask from bgt to clip the buffer.")
    # Add water and railroads in one dataframe.
    bgt_gdf = gpd.GeoDataFrame(
        pd.concat(
            [gdf_water["geometry"], gdf_railroads["geometry"]],
            ignore_index=True,
        )
    ).reset_index()
    bgt_gdf = gpd.GeoDataFrame(
        pd.concat(
            [gdf_water["geometry"], gdf_railroads["geometry"]],
            ignore_index=True,
        )
    ).reset_index()

    logging.info(
        "limit the number of bgt records using a spatial index to increase speed."
    )
    bgt_gdf_limited, bgt_excluded = select_items_with_spatial_index(
        bgt_gdf, buffer_total
    )

    bgt_mask = bgt_gdf_limited.geometry.unary_union

    # Step 4: ============== Clip the buffer with the created mask =============
    logging.info("Clip the buffers (pump and sewer).")
    clipped_buffer = clip_buffer_select_part(buffer_total, bgt_mask, gdf_pump_sewer)

    # export clipped buffer to shape.
    gs_buffer_union = clipped_buffer.geometry.unary_union

    # df_buffer_union = gpd.GeoDataFrame()
    # df_buffer_union = df_buffer_union.set_geometry([gs_buffer_union])
    # df_buffer_union = df_buffer_union.set_crs("EPSG:28992")

    try:
        clipped_buffer.to_file(os.path.join(result_folder, "clipped_buffer.shp"))
    except IOError as error:
        logging.error("saving unary union of buffer failed. %s", error)

    # Step 5: =========== find all plots that are outside the clipped buffer =========
    logging.info("Calculate plots outside the 40m buffer")
    logging.info("Limit plot candidates with spatial index")
    # Here the plots of with the spatial index is outside the buffer are
    # excluded. The others are included.
    # Using the clipped buffer here will result in selecting all plots!!
    # Using the simple unclipped buffer gives the desired result.
    gdf_plots_included, gdf_plots_excluded = select_items_with_spatial_index(
        gdf_plots, buffer_total
    )
    # Here the plots of the included dataframe that are outside the buffer
    # are selected.
    logging.info("Calculating plots outside buffer...")
    gdf_outside_buf = gdf_plots_included[
        gdf_plots_included.geometry.disjoint(gs_buffer_union)
    ]

    # Here all plots outside the buffer are joined in one dataframe
    gdf_outside_buf = gpd.GeoDataFrame(
        pd.concat([gdf_outside_buf, gdf_plots_excluded], ignore_index=True)
    ).reset_index()

    # Export plots outside buffer to shape.
    gdf_outside_buf.to_file(os.path.join(result_folder, "plots_outside_buf.shp"))

    # TODO: Step 6: ============ Select plots that do not have an IBA =================
    # Not implemented yet.
    logging.info(
        "Select plots that do not have an IBA %s",
        "(Individuele Behandeling Afvalwater)",
    )

    # Step 7: ================== Only keep plots with a building which is in use ===============
    logging.info("Add bag information to gdf")
    gdf_outside_buf = gdf_outside_buf.copy(deep=True)

    # This is the check that uses the webservice kadaster - BAG webservice
    # (use either this method, or the original check below):
    # gdf_bag_outside_buf = api_kad.filter_plots_with_bagobject(gdf_outside_buf)
    gdf_bag_outside_buf = pdok_wfs.check_allplots_on_bag(gdf_outside_buf)

    if gdf_bag_outside_buf is None:
        raise ValueError(
            "The building plots shape seems to be empty. Cannot proceed script."
        )

    # export result with bag information to shape
    gdf_bag_outside_buf.to_file(
        os.path.join(result_folder, "plots_bag_outside_buf.shp")
    )


def main():
    """
    Here the work is done. Calculations for municipalities Bodegraven Reeuwijk,
    Haarlemmermeer, Alphen aan den Rijn take a lot of time (huge area). Normally
    it is not possible to run all municipalities at once on a laptop because it
    takes more than one day in total.
    Also downloading plots for big municipalities sometimes results in errors
    like 'too many requests'. Workarounds have been programmed for this, for
    Haarlemmermeer it appeared to work.  But it is not an ideal situation of
    course.

    Variables listmunicipalities and existing municipalities were created to
    split the calculations over serveral days. In ut.merge_shapes all results are
    merged to one shape.
    I always adapt the lists manually if more municipalities are calculated.
    """
    # Initialize the logger. If no logging is required set level lower.
    logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)
    profiler = cProfile.Profile()
    profiler.enable()

    # The municipalities for which the check should be done.
    listmunicipalities = [
        "Leiden",
        # "Gouda",
        # "Lisse",
        # "Teylingen",
        # "Hillegom",
        # "Katwijk",
        # "Wassenaar",
        # "Leiderdorp",
        # "Noordwijk",
        # "Zoeterwoude",
        # "Bodegraven-Reeuwijk",
        # "Haarlemmermeer",
        # "Alphen aan den Rijn",
    ]

    # The municipalities that should be merged to the new output
    # shape with all municipalities. (Kind of reminder.)
    existingmunicipalities = []

    output_folder = os.path.join(os.getcwd(), "data")

    ##########################################################
    # Sometimes the bag service  does not work. This function
    # does this bag analysis separately in case something went
    # wrong.
    # Saves a lot of calculation time because all previous
    # analysis do not need to be done anew then.
    # Exmple here is for Bodegraven-Reeuwijk.
    # In bag service too an try except block is made to prevent
    # a total crash after one error.
    #########################################################
    # ut.get_bag_objects('Bodegraven-Reeuwijk', output_folder)

    #########################################################
    # This is the real analysis, here all work is done.
    #########################################################
    run_opc_per_municipality(
        listmunicipalities, output_folder, use_cache=True, test=False
    )

    ##########################################################
    # In utilities is a function merge_shapes that is meant to
    # merge final results to one shape.
    # The function expects the municipality name as submap and
    # looks for a given shape name in that map.
    # (Final result for every municipality:
    # plots_bag_outside_buf.shp)
    ##########################################################
    # ut.merge_shapes(existingmunicipalities,
    #                 output_folder,
    #                 'plots_bag_outside_buf.shp')

    # ut.merge_shapes(existingmunicipalities,
    #                 output_folder,
    #                 'clipped_buffer.shp')

    # Show some statistics like calculation time etc.
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats("cumulative")
    stats.print_stats(20)


if __name__ == "__main__":
    main()
