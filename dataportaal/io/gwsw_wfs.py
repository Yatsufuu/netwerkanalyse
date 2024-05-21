"""
Functions to read data from GWSW.
"""
import os
import shapely.wkt
from owslib.wfs import WebFeatureService
from requests import Request
import pandas as pd
import geopandas as gpd


# GWSW_NETWERK_WFS_URL = 'https://geodata.gwsw.nl/{}/netwerk'
GWSW_NETWERK_WFS_URL = "https://geodata.gwsw.nl/geoserver/{}-netwerk/wfs"
GWSW_BEHEER_WFS_URL = "https://geodata.gwsw.nl/geoserver/{}-beheer/wfs"
NAME_CONVERSION = {
    "Bodegraven-Reeuwijk": "Bodegraven",
    "Haarlemmermeer": "Haarlemmermeer",
    "Alphen aan den Rijn": "AlphenAanDenRijn",
    "Zoeterwoude": "Zoeterwoude",
    "Gouda": "Gouda",
    "Wassenaar": "Wassenaar",
    "Leiderdorp": "Leiderdorp",
    "Leiden": "Leiden",
    "Noordwijk": "Noordwijk",
    "Lisse": "Lisse",
    "Teylingen": "Teylingen",
    "Hillegom": "Hillegom",
    "Katwijk": "Katwijk",
    "Voorschoten": "Voorschoten",
}


def get_gwsw_sewer_lines(city, data_folder, cache_data=False, filter_value=None):
    """
    Gets data from gwsw WFS and filters only lines that are relevant for
    'ongerioleerde percelen check.' Filters on these two types:
    - 'http://data.gwsw.nl/1.5/totaal/GemengdRiool'
    - 'http://data.gwsw.nl/1.5/totaal/Vuilwaterriool'
    (Type 'http://data.gwsw.nl/1.5/totaal/Drukleiding' is removed because sinks were added.)
    Uses url 'https://geodata.gwsw.nl/<city>/beheer' and feature name
    'gwsw:Beheer_Leiding'.

    Parameters
    ----------
    city: str
        The city name for which the data should be downloaded.
    data_folder: str
        The output folder where the output should be written to.
    cache_data: bool
        If True, the data will be saved as parquet file (2a_sewer_lines.parquet)
    filter_value: array of str, optional
        An array with the specific typenames to download.
        Default value is ['http://data.gwsw.nl/1.5/totaal/GemengdRiool',
                          'http://data.gwsw.nl/1.5/totaal/Vuilwaterriool']

    Returns
    -------
    data: geodataframe
        Contains the requested data in crs EPSG:28992.
    """
    if filter_value is None:
        filter_value = (
            r"http://data.gwsw.nl/\d.\d/totaal/GemengdRiool"
            + r"|http://data.gwsw.nl/\d.\d/totaal/Vuilwaterriool"
        )
        # ['/totaal/GemengdRiool'|
        # '/totaal/Vuilwaterriool'] # ,
        # 'http://data.gwsw.nl/1.5/totaal/Drukleiding']
    layer = "gwsw:beheer_leiding"
    gwsw_city = NAME_CONVERSION[city]
    wfs_url = GWSW_BEHEER_WFS_URL.format(gwsw_city)
    data = get_wfs_data(wfs_url, "2.0.0", layer)
    # data = data[filter_value.isin(data['type'])]
    data = data.loc[data["type"].str.match(filter_value)]
    if cache_data:
        path = os.path.join(data_folder, "2a_sewer_lines.parquet")
        data.to_parquet(path)
    return data


def get_gwsw_pump_points(city, data_folder, cache_data=False, filter_value=None):
    """
    Gets data from gwsw WFS and filters only lines that are relevant for
    'ongerioleerde percelen check.' Filters daefault on these two types:
    ['http://data.gwsw.nl/1.5/totaal/Rioolgemaal',
     'http://data.gwsw.nl/1.5/totaal/Pompunit']

    Uses url 'https://geodata.gwsw.nl/<city>/beheer' and feature name
    'gwsw:Beheer_Pomp'.

    Parameters
    ----------
    city: str
        The city name for which the data should be downloaded.
    data_folder: str
        The output folder where the output should be written to.
    cache_data: bool
        If True, the data will be saved as parquet file (2b_punp_points.parquet)
    filter_value: array of str, optional
        An array with the specific typenames to download.
        Default value is ['http://data.gwsw.nl/1.5/totaal/Rioolgemaal',
                          'http://data.gwsw.nl/1.5/totaal/Pompunit']

    Returns
    -------
    data: geodataframe
        Contains the requested data in crs EPSG:28992.
    """
    if filter_value is None:
        filter_value = (
            r"http://data.gwsw.nl/\d.\d/totaal/Rioolgemaal"
            + r"|http://data.gwsw.nl/\d.\d/totaal/Pompunit"
            + r"|http://data.gwsw.nl/\d.\d/totaal/Gemaal"
            + r"|http://data.gwsw.nl/\d.\d/totaal/Pompput"
        )
    layer = "gwsw:beheer_pomp"
    gwsw_city = NAME_CONVERSION[city]
    wfs_url = GWSW_BEHEER_WFS_URL.format(gwsw_city)
    data = get_wfs_data(wfs_url, "2.0.0", layer)
    # data = data[data['type'].isin(filter_value)]
    data = data.loc[data["type"].str.match(filter_value)]

    if cache_data:
        path = os.path.join(data_folder, "2b_pump_points.parquet")
        data.to_parquet(path)
    return data


def get_wfs_data(url, version, layer, srs="EPSG:28992", outputformat="json"):
    """
    Gets featuredata from a WFS service and returns it a a geodataframe.

    Parameters
    ----------
    url: str
        The base url of the WFS service.
    version: str
        The version of WFS to use.
    layer: str
        The name of the layer to get the data from.
    srs: str, optional
        The projection to use.
        Default is EPSG:28992
    outputformat: str, optional
        Output format the service should provide.
        Default is json.

    Returns
    -------
    data: geodataframe
        A geodataframe with the requested data.
    """
    params = dict(
        service="WFS",
        version=version,
        request="GetFeature",
        typeName=layer,
        srsname=srs,
        outputFormat=outputformat,
    )
    # Parse the URL with parameters
    req_url = Request("GET", url, params=params).prepare().url

    # Read data from URL, gives ERROR: Could not resolve host: geo,
    # however it does return the correct data
    data = gpd.read_file(req_url)
    return data


def get_wfs_connection(wfs_url, version):
    """
    Connects to a WFS service and returns an WFS object and it's metadata.

    Parameters
    ----------
    wfs_url: str
        The url of the WFS service
    version: str
        The WFS version to use.

    Returns
    -------
    wfs: WebFeatureService
        An instance of owslib.wfs.WebFeatureService
    contents: dict
        A dictionary with the metadata of the WFS service.
    """
    wfs = WebFeatureService(url=wfs_url, version=version)
    # Get metadata like layernames etc.
    return wfs, wfs.contents


def save_as_type(path, gdf, output_format="GML"):
    """
    Saves a geodataframe to a given path in given format.

    Parameters
    ----------
    path: str
        String with path for the output file
    gdf: geodataframe
        The geodataframe to export
    output_format: str, optional
        The type of file to export
        Default is GML
        Possible values are: 'AeronavFAA', 'ARCGEN',' BNA', 'DXF', 'CSV',
        'OpenFileGDB','ESRIJSON', 'ESRI Shapefile', 'FlatGeobuf', 'GeoJSON',
        'GeoJSONSeq', 'GPKG', 'GML', 'OGR_GMT', 'GPX', 'GPSTrackMaker',
        'Idrisi', 'MapInfo File', 'DGN', 'PCIDSK', 'OGR_PDS', 'S57', 'SEGY',
        'SUA', 'TopoJSON'

    Returns
    -------
    None
    """
    gdf.to_file(path, driver=output_format)


def read_geo_csv(filepath):
    """
    Reads a geodataframe from a csv file. It is expected that the csv file has
    a column geometry and crs EPSG:28992.

    Parameters
    ----------
    filepath: str
        The path to the csv file

    Returns
    -------
    data: geodataframe
        The data as geodataframe with crs 28992.
    """
    data = pd.read_csv(filepath)
    geometry = data["geometry"].map(shapely.wkt.loads)
    data = gpd.GeoDataFrame(data, crs="EPSG:28992", geometry=geometry)
    return data
