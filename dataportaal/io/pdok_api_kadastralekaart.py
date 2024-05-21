"""
Functions to download kadaster plots.
See explanation here:
https://api.pdok.nl/kadaster/kadastralekaart/download/v4_0/ui/

Every download is a full custom download, delta's are not used here.
The downloaded parcels have a geometry of type:
shapely.geometry.multilinestring.multilinestring

These are converted to polygons, but this proces causes a lot of warnings.
TODO: This part of the process needs a thorough review
"""

import geopandas as gpd
import pandas as pd
import sys
import requests
import wget
import os
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from time import sleep
import logging
import numpy as np
import shapely.ops
from shapely.geometry import MultiLineString
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon


def get_gdf_plots(polygon, data_folder, cache_data=False):
    """
    TODO: Discuss with Yanick to use the new functions ('create plots')
    Returns a dataframe with kadaster plots that lie within a
    given polygon.

    Parameters
    ----------
    polygon: shapely polygon
    data_folder: str
    cache_data: bool, optional

    Returns
    -------
    gdf: geodataframe

    """
    api_url = "https://api.pdok.nl/kadaster/kadastralekaart/download/v5_0/full/custom"

    # TODO api is not functioning anymore..
    # TODO replace by the new WFS server

    gdfs = get_gdf_pdok_api(
        str(polygon),
        api_url,
        layer_names=[
            "kadastralekaart_kadastralegrens.gml",
            "kadastralekaart_perceel.gml",
        ],
        format="gml",
        projection_crs="epsg:28992",
    )
    gdf_grens = gdfs[0]  # Dataframe with all borders as Linestring or Multilinestring
    gdf_perc = gdfs[1]  # Dataframe with administrative info about plots

    # Limit the plots to the desired polygon of the municipality
    polygon = polygon.buffer(
        10
    )  # make sure to include plots that are just on the border
    gdf_perc = gdf_perc.loc[gdf_perc.within(polygon)]
    # Drop the column geometry to prevent two geometry columns in the final result
    gdf_perc = gdf_perc.drop(columns=["geometry"])

    logging.info("Parsing %s percelen, this might take some time", len(gdf_perc.gml_id))
    # gdf = kadastralegrens_to_perceel(gdf_perc, gdf_grens)
    gdf = create_plots(gdf_perc, gdf_grens)

    if cache_data:
        path = os.path.join(data_folder, "1_plots.parquet")
        gdf.to_parquet(path)
    return gdf


def create_polygon(collection):
    """
    Creates a polygon from a collection of Multilinestrings and/or Linstrings
    using the function shapely.ops.polygonize.
    In case of polygons with holes: returns the right outer polygon.

    Parameters
    ----------
    collection: geodataframe
        The dataframe that contains the geometries with Linestrings and/or
        Multilinestrings

    Returns
    -------
    final_result: polygon
        The polygon representation of the input.

    """
    final_result = None
    polygons = list(shapely.ops.polygonize(collection.geometry))
    num_polygons = len(polygons)
    if num_polygons == 1:
        final_result = polygons[0]
    else:
        n = 0
        while (n < num_polygons - 1) | (final_result is None):
            if len(polygons[n].interiors) == (num_polygons - 1):
                final_result = polygons[n]
            n = n + 1
    return final_result


def create_plots(gdf_plots, gdf_borders):
    """
    Receives the downloads from the kadastral website with plots and borders of the plots.
    All borders are shapely.geometry.LineString or shapely.geometry.MultilineString.
    All linestrings are merged and polygonized with shapely function polygonize.

    Parameters
    ----------
    gdf_plots: geodataframe
        The geodataframe containing all information about the plots like id's etc.
    gdf_borders: geodataframe
        The geodataframe with the borders for all plots as downloaded from the kadastral website.

    Returns
    -------
    gdf_result: geodataframe
        A dataframe with columns:
        - gml_id: the Id of the plot
        - geometry: the geometry of the plot as Polygon
        - plus all columns of gdf_plots
    """
    # # create geodataframe for the final result

    gdf_result = gpd.GeoDataFrame(
        columns=["gml_id", "geometry"], geometry="geometry", crs="epsg:28992"
    )

    result = []

    # loop through alle unique plot id's. Some id's only occur in one of these columns.
    # note that since v5.0 the gml_id is defined as NL.IMKAD.KadastraalObject.23970659070000
    # the perceelLinks and perceelRechts only contain the last 14 digits
    for gml_id in gdf_plots["gml_id"]:
        # select all lineparts for the left and the right, merge these to one geometry
        identificatie = int(gml_id.split(".")[-1])  # only get the 14 digits
        collection = gdf_borders.loc[
            (gdf_borders["perceelLinks"] == identificatie)
            | (gdf_borders["perceelRechts"] == identificatie)
        ][["geometry"]]
        try:
            plot = create_polygon(collection)
            row = gpd.GeoDataFrame(
                [{"gml_id": gml_id, "geometry": plot}], crs="epsg:28992"
            )
            gdf_result = gpd.GeoDataFrame(
                pd.concat([gdf_result, row], ignore_index=True)
            )
        except Exception as e:
            # Many plots are partly outside the border, will result in error.
            # The check is left out now probably this increases the speed.
            logging.debug(f"Perceel {gml_id} cannot be polygonized, got error '{e}'.")

    gdf_result = gdf_result.merge(gdf_plots, on="gml_id", suffixes=["", "_perc"])
    return gdf_result.reset_index()


def get_gdf_pdok_api(
    polygon_str,
    api_url,
    layer_names=["kadastralekaart_kadastralegrens.gml"],
    format="gml",
    projection_crs="epsg:28992",
):
    """Download a gml file from PDOK API and read as GeoDataFrame.
    Use polygon to filter the area you want to receive.

    Parameters
    ----------
    polygon_str : str
        polygon shape used as filter
            kadastralekaart accepts format: POLYGON((x1 y1, x2 y2, x3 y3, x1 y1))
            bgt accepts a list of bounding boxes: POLYGON(([xmin, ymin, xmax, ymax]))
    api_url : str
        PDOK API url.For kadastralekaart:
        old; https://downloads.pdok.nl/kadastralekaart/api/v4_0/full/custom
        new; https://api.pdok.nl/kadaster/kadastralekaart/download/v5_0/full/custom
    layer_names: array of str
        The names of the output files for the layers in the request. In the code
        the featuretypes are 'perceel' and 'kadastralegrens'.
    format: str
        The format in which the data is downloaded. Default = 'gml'
    projection_crs: str
        The epsg code of the desired projection. Default = 'epsg:28992'

    Returns
    -------
    gdfs: array of geodataframe
        Containing the requested PDOK data in two dataframes
        (one containing plots, another containing the border of the frames)
    """
    logging.info(f"Send download request to PDOK kadastralekaart API server")
    # build request to API server
    r = requests.post(
        api_url,
        json={
            "featuretypes": ["perceel", "kadastralegrens"],
            "format": format,
            "geofilter": polygon_str,
        },
    )
    logging.info(f"response status code {r.status_code}")

    # check response
    result = r.json()
    downloadID = result["downloadRequestId"]

    sleep(5)

    # Get status response
    url_download = api_url + f"/{downloadID}/status"
    response = requests.get(url_download)
    status = response.json()

    # Check progress status
    while status["status"] != "COMPLETED":
        response = requests.get(url_download)
        status = response.json()
        ## This try except is build in by peter venema because suddenly errors
        ## occured during download Haarlemmermeer (> 60.000 plots).
        ## The error message was 'Too many requests'
        try:
            if divmod(status["progress"], 1)[1] == 0:
                job_progress("status of job: {} %".format(status["progress"]))
            logging.info(f"status response is {status['status']}")
        except Exception as e:
            ## This is build in as workaround if server says too many requests.
            ## Danger is that we will be blacklisted.
            logging.error("Exception occured while downloading plots, %s", status)
            status["status"] = "RUNNING"

    # Perform zipfile download
    with TemporaryDirectory() as temp_dir:
        if status["status"] == "COMPLETED":
            logging.info(f"Download zipfile")
            pdok_bgt_theme_zip = wget.download(
                f"https://api.pdok.nl{status['_links']['download']['href']}",
                os.path.join(temp_dir, "pdok_download.zip"),
                bar=download_progress,
            )

        # Unpack downloaded zip file containing wanted gml's
        with ZipFile(pdok_bgt_theme_zip, "r") as zipObj:
            logging.info(f"Extract zipfile")
            zipObj.extractall()

        gdfs = []
        logging.info(f"Merge layers")
        for layer in layer_names:
            gdf = gpd.read_file(layer)
            gdf.crs = projection_crs
            gdfs.append(gdf)
    return gdfs


def job_progress(progress_message):
    """Flushes progress message to stdout, i.e. a dynamic
    changing message appears as notebook cell output.

    Parameters
    ----------
    progress_message : str
        Progress message text
    """
    sys.stdout.write("\r" + progress_message)
    sys.stdout.flush()


def download_progress(current, total, width=80):
    """Flush download progress message to stdout, i.e.
    a dynamic changing one line message appears as
    notebook cell output.

    Parameters
    ----------
    current : float
        downloaded part of file so far
    total : float
        total size of file to download
    width : int, optional
        maximum allowed line width, by default 80
    """
    progress_message = "Downloading: %d%% [%d / %d] bytes" % (
        current / total * 100,
        current,
        total,
    )
    sys.stdout.write("\r" + progress_message)
    sys.stdout.flush()


###################################################################################################
# OLD FUNCTIONS, TROW AWAY IN NEXT UPDATE PLEASE.
###################################################################################################


def old_kadastralegrens_to_perceel(gdf_perc, gdf_grens):
    """
    NOT USED ANYMORE
    TODO: @Yanick Mampaey: please explain better what happens here.

    Creates polygons from two geodataframes from BRK.

    Parameters
    ----------
    gdf_perc: geodataframe
        The dataframe that contains the id of a plot. Geometry is
        a point.
    gdf_grens: geodataframe
        The dataframe that contains the actual coordinates of the
        border of the plot. Geometry is of type
        shapely.geometry.multilinestring.multilinestring

    Returns
    -------
    gdf: geodataframe
        A geodataframe with all plots as polygons.

    """
    gdf = gdf_perc.copy(deep=True)

    # gdf = gdf.explode() # to get rid of multipolygons --> no effect :-(
    def get_cadasteral_boundaries(row):
        """
        Internal function, is called as lambda function  in
        apply on dataframe gdf_perc, a row from that
        dataframe is passed.

        Parameters
        ----------
        row: dataframe row
            The row that contains a shapely geometry.

        Returns
        -------
        gdf: geodataframe
        """
        geometry = gdf_grens[
            (gdf_grens["perceelLinks|gml_id"].values == row.gml_id)
            | (gdf_grens["perceelRechts|gml_id"].values == row.gml_id)
        ].geometry.unary_union
        # geometry is either a linestring or a multilinestring and need to be handeled seperately

        try:  # for a simple linestring.
            # this will always go wrong so seems to make no sense.
            #     poly_coordinates = list(geometry.coords)
            # except Exception as e:
            # for the multilinestring (pve: the geometry is a multilinestring!)
            poly_coordinates = list()
            for line in geometry.geoms:
                if not poly_coordinates:
                    poly_coordinates = [*list(line.coords)]
                    remaining_linestrings = list(geometry[1:])
                else:
                    poly_coordinates, remaining_lines = find_next_linestring(
                        poly_coordinates, remaining_linestrings
                    )
        except Exception as e:
            logging.error("Exception occured : %s", e)
            # poly_coordinates = list()
            # for line in geometry.geoms:
            # poly_coordinates = [*poly_coordinates, *list(line.coords)]
        try:
            polygon_plot = Polygon(poly_coordinates)
        except Exception as e:
            # logging.warning(f'Failed creating polygon of geometry {geometry} of gml_id {row.gml_id}. Error: {e}')
            logging.warning("tst multipolygon")
            polygon_plot = MultiPolygon(poly_coordinates)
            # polygon_plot = np.NaN
            # TODO: code for multilines will do the trick

        # import matplotlib.pyplot as plt; plt.plot(*polygon_plot.exterior.xy); plt.show()
        return polygon_plot

    gdf["geometry"] = gdf.apply(lambda row: get_cadasteral_boundaries(row), axis=1)
    return gdf


def old_find_next_linestring(poly_coordinates, remaining_linestrings):
    """
    NOT USED ANYMORE
    The geometry of get_cadasteral_boundaries gives a multilinestring that
    is not ordered correctly for a polygon. This functions finds the next
    linestring in the sequence and returns the polygon coordinates and
    the remaining linestrings

    Parameters
    ----------
    poly_coordinates : list
        containing coordinate points as tuples
    remaining_linestrings : list
        containing shapely lines

    Returns
    -------
    poly_coordinates : list
        containing coordinate points as tuples, with points of the next line added
    remaining_linestrings : list
        containing shapely lines, with one line substracted
    """
    counter = 0
    for line in remaining_linestrings:
        if poly_coordinates[-1] == tuple(line.coords[0]):
            poly_coordinates = [*poly_coordinates, *list(line.coords)[1:]]
            remaining_linestrings.pop(counter)
            break
        if poly_coordinates[-1] == tuple(line.coords[-1]):
            reverse_linestring = list(line.coords)
            reverse_linestring.reverse()
            poly_coordinates = [*poly_coordinates, *reverse_linestring[1:]]
            remaining_linestrings.pop(counter)
            break
        counter += 1
    return poly_coordinates, remaining_linestrings


def old_merge_lines(gdf):
    """
    NOT USED
    Receives a collection of geometries of the left part and the right
    part of a plot border. In de kadastral plot export these geometries consist
    of shapely.geometry.LineString and shapely.geometry.MultiLineString.
    All separate linestrings are merged to one linestring if possible.
    Otherwise a multilinestring will be created.

    Parameters
    ----------
    gdf: geodataframe
        The geodataframe with the linestring parts for the left and right part of the borders

    Returns
    -------
    final_linestring: shapely linestring or multilinestring
        The final result is mostly a linestring but in some situations a multiline
        string is returned.
    """
    # print(shapely.ops.linemerge(gdf_left))
    result = []
    for geom in gdf["geometry"]:
        result.append(geom)
    final_linestring = shapely.ops.linemerge(result)
    return final_linestring


def old_create_plots(gdf_plots, gdf_borders, polygon):
    """
    NOT USED ANYMORE, THERE WERE MANY PROBLEMS WITH THE PLOTS. THEREFORE I WANTED TO
    KEEP THIS CODE AT HAND. KAN BE REMOVED LATER IF THE CODE APPEARS TO BE STABLE.

    Comment PVE: this code is more simple and straightforward and causes no warnings. Is
    a bit slower, but made quicker to select only linestring within the polygon.
    Lisse takes 5 minutes calculation time.


    TODO: I this method is used the gdf_grens should be used to obtain the administrative
    kadastral info. In this function this is not implemented. The function is meant to get
    better understanding of the downloaded gml.

    Receives the download from the kadastral website with borders of the plots. All borders are
    shapely.geometry.LineString.
    All linestrings are merged, for every merged linestring a check is done if it is a closed
    linestring. If so, a polygon is created. The polygon is stored with the plot id in a row of
    a geodataframe.
    In case of a multilinestring, a check is done on each individual linestring if it is closed.
    If so, a Polygon is created. All resulting polygons are then merged into a multipolygon and
    added to the final result.

    Parameters
    ----------
    gdf_borders: geodataframe
        The geodataframe with the borders for all plots as downloaded from the kadastral website.

    Returns
    -------
    gdf_plots: geodataframe
        A dataframe with two columns:
        - perceelID: the Id of the plot
        - geometry: the geometry of the plot as Polygon
    """
    # # create geodataframe for the final result
    gdf_result = gpd.GeoDataFrame(
        columns=["gml_id", "geometry"], geometry="geometry", crs="epsg:28992"
    )
    polygon = polygon.buffer(10)
    # gdf_plots = gdf_plots.loc[gdf_plots.within(polygon)].copy()
    # gdf_plots = gdf_plots.rename(columns={'geometry': 'point_geom'})
    gdf_plots = gdf_plots.drop(columns=["geometry"])
    # loop through alle unique plot id's. Some id's only occur in one of these columns.
    for id in gdf_plots["gml_id"]:
        # select all lineparts for the left and the right, merge these to one geometry
        collection = gdf_borders.loc[
            (gdf_borders["perceelLinks|gml_id"] == id)
            | (gdf_borders["perceelRechts|gml_id"] == id)
        ][["geometry"]]
        ls = merge_lines(collection)
        if ls.within(polygon):
            geom = None
            if ls.is_closed:
                # if the geometry is closed, create a polygon
                geom = Polygon(ls.coords)
                # gdf_plots.loc[gdf_plots['gml_id'] == id][['geometry']] = geom
            elif type(ls) is MultiLineString:
                # In case of a multilinestring check each individual linestring
                # if it is closed and if yes, create a polygon with holes.
                # Note: first a Multipolygon was created but that's not right. It will
                # lead to errors.
                # It appears that linestring occur that are not closed. This makes the
                # code a bit more complicated but a check is_closed is needed!
                outerbound = None
                holes = []
                for hole in ls:
                    # check if the hole is closed
                    if id == 23150428970000:
                        print(f"\n hole: {hole} \n {hole.is_closed}\n\n")
                    if hole.is_closed:
                        if outerbound is None:
                            outerbound = hole
                        else:
                            if id == 23150428970000:
                                print(f"outerbound = {outerbound}\n hole = {hole}")

                            if Polygon(outerbound).contains(hole):
                                holes.append(hole.coords)
                            else:
                                # if not then the hole must be the outerbound
                                holes.append(outerbound.coords)
                                outerbound = hole

                geom = Polygon(outerbound.coords, holes)
                # gdf_plots.loc[gdf_plots['gml_id'] == id][['geometry']] = gpd.GeoSeries([geom])
            row = {"gml_id": id, "geometry": geom}
            gdf_result = gdf_result.append(row, ignore_index=True)
    gdf_result = gdf_result.merge(gdf_plots, on="gml_id", suffixes=["", "_perc"])
    test = Polygon
    test.interiors
    return gdf_result
