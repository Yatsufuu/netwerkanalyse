"""Contains functions to download data from the BGT download service.
"""
from http.client import HTTPException
import time
import os
import zipfile
import tempfile
import logging
import geopandas as gpd
from osgeo import gdal
import wget
import requests

import dataportaal.io.pdok_wfs as pdok_wfs


def get_bgt_features_mun(
    municipality, outputpath="", cache_data=True, features='["wegdeel", "waterdeel"]'
):
    """
    Gets feature layers from bgt download for a
    given municipality. Saves the data as parquet to a given output path.
    Uses function 'get_wkt_bbox_municipality(mun)' from pdok_wfs to get the
    bounding box of the municipality. After that calls function
    'get_bgt_features_poly' to actually retrieve the data.

    Parameters
    ----------
    municipality: str
        The name of the municipality.

    outputpath: str
        The path to which the data should be written as parquet

    cache_data: bool, optional
        If True the data are saved to the outputpath as parquet
        Default value = True

    features: str, optional
        A string containing the names of the features to get.
        Default = '["wegdeel","waterdeel"]'

    Returns
    -------
    data: dict of geodataframes
        Contains for each feature layer a geodataframe
    """
    wkt_bbox = pdok_wfs.get_wkt_bbox_municipality(municipality)
    data = get_bgt_features_poly(wkt_bbox, outputpath, cache_data, features)
    return data


def get_bgt_features_poly(
    selectpolygon, outputpath="", cache_data=True, features='["wegdeel", "waterdeel"]'
):
    """
    Gets feature layers from bgt download for a
    given wkt polygon. Saves the data as parquet to a given output path if
    cache_data=True.
    The raw data is saved to a internal temp directory, so not directly
    available.

    Parameters
    ----------
    selectpolygon: str
        The polygon in which the features are selected.
        In wkt format.

    outputpath: str
        The path to which the data should be written as parquet

    cache_data: bool, optional
        If True, the data will be saved to the outputpath as parquet.
        Default value = True

    features: str, optional
        A string containing the names of the features to get.
        Default = '["wegdeel","waterdeel"]'

    Returns
    -------
    data: dict of geodataframes
        Contains for each feature layer a geodataframe
    """
    temp_dir = tempfile.TemporaryDirectory()
    temp_path = temp_dir.name
    zippath = os.path.join(temp_path, "bgt_extract.zip")

    file_path = download_bgt_data(selectpolygon, zippath, features)
    returnvalue = {}

    if file_path is not None:
        # Retrieve all gml data from the downloaded zip file.
        # Convert the citygml to normal gml and read it as gdf.
        with zipfile.ZipFile(zippath) as z_file:
            z_file.extractall(temp_path)
            for file in z_file.infolist():
                input_file = file.filename
                output_file = "conv_" + input_file
                gdal.VectorTranslate(
                    os.path.join(temp_path, output_file),
                    os.path.join(temp_path, input_file),
                    options='-f "gml" -nlt CONVERT_TO_LINEAR',  # converting curved polygons to lineair to avoid fiona read error
                )
                # raises error, unsupported type 10 -> curvepolygon
                # https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry#Well-known_binary
                # solution: https://gis.stackexchange.com/questions/433216/read-xml-with-curvepolygon-with-geopandas-fiona
                gdf = gpd.read_file(os.path.join(temp_path, output_file))
                if cache_data:
                    gdf.to_parquet(
                        os.path.join(outputpath, input_file.replace("gml", "parquet"))
                    )
                returnvalue[input_file.replace(".gml", "")] = gdf
    return returnvalue


def download_bgt_data(selectpolygon, outputpath, features):
    """
    Downloads bgt features of a given selection polygon (wkt) as
    citygml (imgeo) in a zip file to a given outputpath.

    Parameters
    ----------
    selectpolygon: str
        The polygon of the selection (wkt format)

    outputpath: str
        The path to save the zip file to

    features: str
        The features to download. This is a string containing an
        array with featurenames (e.g.: '["wegdeel", "waterdeel"]')

    Returns
    -------
    file_url: str
        The path to the donwloaded file, in case of error None is returned.
    """
    file_url = None
    # These are the eindpoints for the BGT API
    baseurl = "https://api.pdok.nl"
    requesturi = "/lv/bgt/download/v1_0/full/custom"
    statusuri = "/lv/bgt/download/v1_0/full/custom/-id-/status"

    # Here the parameters for the request are prepared.
    data = (
        '{"featuretypes": -features-,'
        + '"format":"citygml",'
        + '"geofilter": "-polygon-"}'
    )
    data = data.replace("-features-", features)
    data = data.replace("-polygon-", selectpolygon)

    # Step 1 is send the request to the BGT api
    resp_download_id = requests.post(
        baseurl + requesturi,
        data=data,
        headers={"accept": "application/json", "Content-Type": "application/json"},
    )
    # Step 2 is check if the request is accepted and wait until
    #        the data are prepared.
    if resp_download_id.status_code == 202:
        download_id = resp_download_id.json()["downloadRequestId"]
        statusurl = baseurl + statusuri.replace("-id-", download_id)
        try:
            resp_download_url = requests.get(statusurl)
            while resp_download_url.json()["status"] != "COMPLETED":
                logging.info(
                    "Request in progress (%s)", resp_download_url.json()["progress"]
                )
                time.sleep(1)
                resp_download_url = requests.get(statusurl)
        except requests.exceptions.RequestException as ex:
            logging.error("Error during preparing BGT data. {%s}", ex)

        try:
            downloaduri = resp_download_url.json()["_links"]["download"]["href"]
            downloadurl = baseurl + downloaduri
            file_url = wget.download(downloadurl, outputpath)

        except HTTPException as ex:
            logging.error("Error requesting BGT data {%s}.", ex)
    else:
        logging.info("Fout bij versturen aanvraag.")

    return file_url


if __name__ == "__main__":
    # test = get_bgt_features_mun("Leiden", os.path.join(os.getcwd(),'research', 'data'))
    gdf_result = gpd.read_parquet(
        os.path.join(os.getcwd(), "research", "data", "bgt_wegdeel.parquet")
    )
    print(gdf_result["function"].unique())
