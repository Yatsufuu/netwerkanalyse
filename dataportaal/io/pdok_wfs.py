"""This file contains functions to read data from
pdok wfs for administrative areas (bestuurlijke gebieden)
and the Administration Adresses and Buildings (BAG -
Basisregistratie Adressen en Gebouwen).


Ter info:
Per 1 juli 2021 is de opvolger van Bestuurlijke grenzen live gezet namelijk de
Bestuurlijke gebieden.
De dataset Bestuurlijke grenzen wordt niet meer bijgewerkt;
de WMS en WFS services blijven tot eind 2021 live, daarna worden deze
definitief uit productie genomen.
De downloads blijven wel beschikbaar maar zullen per jaar verdwijnen.
Nieuwe url:
https://service.pdok.nl/kadaster/bestuurlijkegebieden/wfs
                                       /v1_0?request=GetCapabilities&service=WFS
"""
import logging
import xml.etree.ElementTree as ET
import geopandas as gpd
import requests
from shapely.ops import unary_union
import numpy as np
import shapely
import owslib.fes as fes
from pathlib import Path
import pandas as pd

BEST_GR_WFS = "https://service.pdok.nl/kadaster/bestuurlijkegebieden/wfs/v1_0"
GEMEENTEN_LAYER = "bestuurlijkegebieden:Gemeentegebied"

BAG_WFS = "https://service.pdok.nl/lv/bag/wfs/v2_0"
# pdok_bag_url_wfs = "https://geodata.nationaalgeoregister.nl/bag/wfs/v1_1"


def get_wktpolygon_municipality(name):
    """
    Returns wkt geometry for a given municipality from the pdok
    WFS service 'bestuurlijkegrenzen'.
    The filter for the name is applied directly to the WFS request.

    Function is obsolete, bhbox is more usefull.

    Parameters
    ----------
    name: str
        The name of the municipality

    Returns: array of str
        The strings contain the wkt Polygon(s) for the given municipality.
    """
    version = "2.0.0"
    srs = "EPSG:28992"
    outputformat = "json"

    # Apply the filter to the wfs request
    filter1 = fes.PropertyIsEqualTo("naam", name)
    f_r = fes.FilterRequest()
    filter_fes = f_r.setConstraint(filter1, tostring=True)

    # Set all parameters for the wfs request
    params = dict(
        service="WFS",
        version=version,
        request="GetFeature",
        typeName=GEMEENTEN_LAYER,
        srsname=srs,
        filter=filter_fes,
        outputFormat=outputformat,
    )
    # Parse the URL with parameters
    req_url = requests.Request("GET", BEST_GR_WFS, params=params).prepare().url

    # Read data from URL, gives ERROR: Could not resolve host: geo,
    # however it does return the correct data
    data = gpd.read_file(req_url)
    polygon = unary_union(data["geometry"])

    # polygon_list = []
    # polygons = json.\
    #     loads(data.to_json())['features'][0]['geometry']['coordinates']
    # for polygon in polygons:
    #     wktpolygon = create_wkt_polygon(polygon)
    #     polygon_list.append(wktpolygon)
    return polygon


def get_wkt_bbox_municipality(mun):
    """
    Returns wkt geometry of the bbox for a given municipality from
    the pdok WFS service 'bestuurlijkegrenzen'.
    The filter for the name is applied directly to the WFS request.

    Parameters
    ----------
    name: str
        The name of the municipality

    Returns: array of str
        The string contain the wkt polygon of the bbox of the
        polygon(s) for the given municipality.
    """
    version = "2.0.0"
    srs = "EPSG:28992"
    outputformat = "json"

    # Apply the filter to the WFS request
    filter1 = fes.PropertyIsEqualTo("naam", mun)
    f_r = fes.FilterRequest()
    filter_fes = f_r.setConstraint(filter1, tostring=True)

    # Fill the parameters of the WFS request
    params = dict(
        service="WFS",
        version=version,
        request="GetFeature",
        typeName=GEMEENTEN_LAYER,
        srsname=srs,
        filter=filter_fes,
        outputFormat=outputformat,
    )
    # Parse the URL with parameters
    req_url = requests.Request("GET", BEST_GR_WFS, params=params).prepare().url
    # print(req_url)

    # Read data from URL, gives ERROR: Could not resolve host: geo,
    # however it does return the correct data
    data = gpd.read_file(req_url)
    bbox = data.geometry.bounds
    xmin = bbox.iloc[0]["minx"]
    ymin = bbox.iloc[0]["miny"]
    xmax = bbox.iloc[0]["maxx"]
    ymax = bbox.iloc[0]["maxy"]

    bbox_municipality = (
        f"POLYGON(({xmin} {ymin},{xmin} {ymax},{xmax} {ymax},"
        + f"{xmax} {ymin},{xmin} {ymin}))"
    )
    return bbox_municipality


def create_wkt_polygon(polygon):
    """
    Receives polygon coordinates in json format and returns a string
    containing a wkt polygon.

    Parameters
    ----------
    polygon: json
        A list of coordinates of a polygon.

    Returns
    -------
    wktpolygon: str
        A str containing a wkt polygon with all coordinates.
    """
    wktpolygon = "POLYGON(("
    for coordinate in polygon[0]:
        wktpolygon += str(coordinate[0]) + " " + str(coordinate[1]) + ","
    wktpolygon = wktpolygon[:-1] + "))"
    return wktpolygon


def get_contour_municipalities(municipality, data_folder, use_cache):
    """Get a GeoDataFrame with municipality contours as polygons
    from a list of municipality strings

    Retrieves information via the PDOK WFS service. Source:
    https://www.pdok.nl/introductie/-/article/bestuurlijke-grenzen

    Parameters
    ----------
    municipalities : list
        Containing municipalities as strings, eg: ['Leiden', 'Leiderdorp']
    data_folder : str
        Folder to write cache data to
    use_cache: bool
        If true, it will be tried to read data from disk. If in that
        case no data is found, it will be downloaded.

    Returns
    -------
    GeoDataFrame
        Containing the geometries of municipalities
    """

    contour_filename = Path(data_folder) / f"contour_municipality_{municipality}.json"

    # if we already have the file we will use that one
    # unless we ask for an update
    if not contour_filename.exists() or not use_cache:
        logging.debug("Start downloading municipality: %s", municipality)

        # Apply the filter to the wfs request
        filter1 = fes.PropertyIsEqualTo("naam", municipality)
        f_r = fes.FilterRequest()
        filter_fes = f_r.setConstraint(filter1, tostring=True)

        params = dict(
            service="WFS",
            request="GetFeature",
            typeName=GEMEENTEN_LAYER,
            srsname="EPSG:28992",
            filter=filter_fes,
            version="2.0.0",
            outputFormat="json",
        )
        # Parse the URL with parameters
        req_url = requests.Request("GET", BEST_GR_WFS, params=params).prepare().url

        # Download as a file
        response = requests.get(req_url)
        with open(contour_filename, "w") as f:
            f.write(response.text)

    # Read data from URL
    gdf = gpd.read_file(contour_filename)
    # gdf_filter = gdf[gdf['naam'] == municipality]
    # alternative method for selecting multiple municipalities:
    # gdf_filter = gpd.GeoDataFrame()
    # for municipality in municipalities:
    #     gdf_filter = gdf_filter.\
    #           append(gdf[gdf['gemeentenaam'] == municipality])
    polygon_municipalities = unary_union(gdf["geometry"])
    return gdf, polygon_municipalities


def check_allplots_on_bag(gdf):
    """
    Receives a geodataframe with plots, and checks for each plot on bag wfs
    if it contains address information. All found addresses are stored in a
    dataframe that is returned.

    TODO: Find out why sometimes the service crashes on a timeout. Maybe
    because requests are to quick after each other??

    Parameters
    ----------
    gdf: Geodataframe
        A geodataframe that contains plot data, read from the kadaster webservice.

    Returns
    -------
    eindresultaat: geodataframe
        A geodataframe with the plot polygons and all addresses for each plot.
    """
    tot_rows = len(gdf)
    eindresultaat = None
    for count in range(0, tot_rows - 1):
        row = gdf.iloc[count]
        try:
            check_result = check_plot_on_bag(row)
        except:
            logging.error(
                "Error in plot with localID: %s. No check performed", row["lokaalID"]
            )

        if check_result is not None:
            if eindresultaat is None:
                eindresultaat = check_result.copy()
            else:
                eindresultaat = gpd.GeoDataFrame(
                    pd.concat([eindresultaat, check_result], ignore_index=True)
                )
    return eindresultaat.reset_index()


def check_plot_on_bag(plot_row):
    """
    Receives a row from the plots dataframe and checks if the geometry contains
    bag locatieons with function 'get_bag_locations'.
    - First checks if an addresslocation is found, if not:
    - Checks if an 'ligplaats' is found, if not:
    - Checks if an 'standplaats' is found.

    Parameters
    ----------
    plot_row: dataframe.row
        The row of the dataframe with the plot that should be checked.

    Returns
    -------
    gdf_bag_result: dataframe or None if nothing found
        A dataframe that contains all addressess found. Each row in the dataframe contains
        the geometry of the plot.
    """
    # In kadaster data some very long prefixes occur, this is one of them:
    c_gemcode = (
        "kadastraleAanduiding|TypeKadastraleAanduiding|"
        + "aKRKadastraleGemeenteCode|AKRKadastraleGemeenteCode|waarde"
    )
    c_perceelnr = "perceelnummer"
    c_letter = "sectie"

    kadasterid = (
        plot_row[c_gemcode]
        + "-"
        + plot_row[c_letter]
        + "-"
        + str(plot_row[c_perceelnr])
    )
    kadaster_locaalid = plot_row["lokaalID"]

    geometry = plot_row.geometry
    gdf_bag_result = get_bag_locations(geometry)

    if gdf_bag_result is None:
        gdf_bag_result = get_bag_locations(geometry, "ligplaats")
        if gdf_bag_result is None:
            gdf_bag_result = get_bag_locations(geometry, "standplaats")
            if gdf_bag_result is not None:
                gdf_bag_result["gebruiksdoel"] = "Standplaats"
        else:
            gdf_bag_result["gebruiksdoel"] = "Ligplaats"

    if gdf_bag_result is not None:
        gdf_bag_result = gdf_bag_result.set_crs("EPSG:28992")
        gdf_bag_result = gdf_bag_result[
            [
                "openbare_ruimte",
                "huisnummer",
                "huisletter",
                "toevoeging",
                "postcode",
                "woonplaats",
                "gebruiksdoel",
                "geometry",
            ]
        ]
        gdf_bag_result["geometry"] = geometry
        gdf_bag_result["kadaster"] = kadasterid
        gdf_bag_result["kadaster_locaalid"] = kadaster_locaalid

    return gdf_bag_result


def create_gml_multipolygon(geometry):
    """
    Receives a shapely MultiPolygon geometry and converts it to
    a gml Multipolygon that can serve as a filter in a WFS request.

    Parameters
    ----------
    geometry: shapely.geometry.Multipolygon
        The multipolygon to convert

    Returns
    -------
    gml_multipolygon: str
        A gml representation of the multipolygon
    """
    polygon_members = ""

    for polygon in list(geometry.geoms):
        gml_polygon = create_gml_polygon(polygon, False)
        polygon_members += f"""
        <gml:polygonMember>
            {gml_polygon}
        </gml:polygonMember>
        """

    gml_multipolygon = f"""<gml:MultiPolygon gml:id="filter" srsName="urn:ogc:def:crs:EPSG::28992">
                {polygon_members}
	        </gml:MultiPolygon>"""
    return gml_multipolygon


def create_gml_linearing(coords):
    """
    Receives coordinates from a shapely geometry and create an geml linear ring from them.

    Parameters
    ----------
    coords: array
        The coordinates of the geometry that should be converted to a gml linearRing.

    Returns
    -------
    linearRing: str
        The gml string with the coordinates as linearRing
    """

    coordinates = ""
    for coordinate in coords:
        coordinates += str(coordinate[0]) + " " + str(coordinate[1]) + " "

    coordinates = coordinates[:-1]
    linearRing = f"""<gml:LinearRing><gml:posList srsDimension="2">
    {coordinates}</gml:posList></gml:LinearRing>"""
    return linearRing


def create_gml_interiors(interiors):
    """
    Receives a shapely geometry.interiors and converts it to
    a gml interiors that can serve as a filter in a WFS request.

    Parameters
    ----------
    interiors: geometry.interiors
        The interiors

    Returns
    -------
    interiors: str
        A gml representation of the interiors (holes)
    """
    interiors_gml = ""

    for interior in interiors:
        linRing = create_gml_linearing(interior.coords)
        interiors_gml += f"""
        <gml:interior>{linRing}</gml:interior>
        """
    return interiors_gml


def create_gml_polygon(geometry, filter):
    """
    Receives a shapely Polygon geometry and converts this to
    a gml polygon. It can serve as part of a MultiPolygon or as
    a filter itself.

    Parameters
    ----------
    geometry: shapely.geometry.Polygon
        The polygon to convert
    filter: bool
        If filter is True, an attribute filter is added.
        This is needed if the polygon itself is the filter.

    Returns
    -------
    gml_polygon: str
        A gml presentation of the polygon

    """
    extRing = create_gml_linearing(geometry.exterior.coords)
    interiors = create_gml_interiors(geometry.interiors)

    filter_attribute = ""
    if filter:
        filter_attribute = 'gml:id="filter"'

    gml_polygon = f"""
        <gml:Polygon {filter_attribute} srsName="urn:ogc:def:crs:EPSG::28992">
        <gml:exterior>
        {extRing}
        </gml:exterior>
        {interiors}
        </gml:Polygon>
        """
    return gml_polygon


def get_bag_locations(geometry, typename="verblijfsobject"):
    """
    Retrieves from the pdok WFS service for bag a geodataframe with bag
    locations that meet two requirements:
    - Lie within a given polygon (geometry)
    - AND are in use (status = 'Verblijfsobject in gebruik' or
                      status = 'Plaats aangewezen')
    Uses a post request.
    See: https://www.kadaster.nl/-/beschrijving-bag-wfs

    Parameters:
    geometry: Polygon
        The polygon that filters the bag locations.
    typename: str
        The typename(s) of the layer for which objects should be returned.
        If more names are used, separate them with '&'
    """
    status = "Verblijfsobject in gebruik"
    if typename != "verblijfsobject":
        status = "Plaats aangewezen"

    filter_geometry = create_gml_polygon(geometry, True)

    # Intersects vervangen door Within
    postheader = f"""<?xml version="1.0" encoding="utf-8"?>
        <GetFeature xmlns="http://www.opengis.net/wfs/2.0" xmlns:gml="http://www.opengis.ne
        t/gml/3.2" service="WFS" version="2.0.0" outputformat="json"
        xmlns:xsi="http://www.w3.org/2001/XMLSchem
        ainstance" xsi:schemaLocation="http://schemas.opengis.net/wfs/2.0/wfs.xsd http://sch
        emas.opengis.net/wfs/2.0.0/WFS-transaction.xsd">
        <Query typeNames="{typename}" xmlns:bag="http://bag.geonovum.nl">
        <fes:Filter xmlns:fes="http://www.opengis.net/fes/2.0">
        <fes:And>
        <fes:Within>
        <fes:ValueReference>geometrie</fes:ValueReference>
        {filter_geometry}
        </fes:Within>
        <fes:PropertyIsEqualTo>
        <fes:PropertyName>status</fes:PropertyName>
        <fes:Literal>{status}</fes:Literal>
        </fes:PropertyIsEqualTo>
        </fes:And>
        </fes:Filter>
        </Query>
        </GetFeature>
        """
    # Step 1 is send the request to the BGT api
    req_url = requests.post(
        BAG_WFS,
        data=postheader,
        headers={"accept": "application/xml", "Content-Type": "application/json"},
    )
    gdf = None
    if req_url.status_code == 200:  # Status = O.k.
        # print(req_url.json())
        gdf = gpd.GeoDataFrame.from_features(req_url.json())
        if gdf.empty:
            gdf = None

    return gdf


def get_gdf_wfs(url, gdf, layers, sortby, n_cells=30):
    """
    Gets all data from BAG WFS using function 'get_gdf_wfs_loop'

    Available layers:
    https://www.nationaalgeoregister.nl/geonetwork
    /srv/dut/catalog.search#/metadata/1c0dcc64-91aa-4d44-a9e3-54355556f5e7  # nopep8

    Parameters
    ----------
    gdf : GeoDataFrame
        containing the shape of the area to load through the WFS service
    data_folder : str
        folder for cache data
    n_cells : int, optional
        amount of cells to divide the dataframe in, by default 30
    cache_data : bool, optional
        Caching gdf to data_folder if True, by default False
    """
    gdf_total = gpd.GeoDataFrame()
    for layer in layers:
        gdf_singlelayer = get_gdf_wfs_loop(url, gdf, layer, sortby, n_cells=n_cells)
        gdf_singlelayer = gdf_singlelayer.drop_duplicates()
        gdf_singlelayer["layer"] = layer
        gdf_total = gdf_total.append(gdf_singlelayer)
        gdf_total = gdf_total.drop_duplicates(subset=sortby)
    gdf_total = gdf_total.reset_index(drop=True)
    return gdf_total


def get_gdf_wfs_loop(url, gdf, layer, sortby, n_cells=30):
    """
    Read all BAG data for a specified layer from the pdok BAG WFS service
    that can be found in the bounds of a given geodataframe.

    Parameters
    ----------
    url: str
        The base url of the bag wfs service.
    gdf: geodataframe
        The geodataframe of which the bounds are used to limit the request.
    sortby:

    """
    # pdok_bag_url_wfs = "https://geodata.nationaalgeoregister.nl/bag/wfs/v1_1"
    type_name = layer
    coordinate_sys = "EPSG:28992"  # 'EPSG:4326'
    count = "1000"

    grid_cells = get_grid_cells_gdf(gdf, n_cells)
    # select only the cells in the municipality geometry
    grid_cells = grid_cells[~grid_cells.geometry.disjoint(gdf.geometry.unary_union)]
    # Uncomment to check the area to download:
    # fig, ax = plt.subplots(figsize=(15, 15))
    # grid_cells.plot(ax=ax)
    # gdf.plot(ax=ax, alpha=0.9, color="pink")

    # this used to be lambda x: shapely.wkt.dumps(x)
    grid_geometry_strings = grid_cells.geometry.apply(shapely.wkt.dumps)
    # drop the 'POLYGON ((...))' part of the geometry string to comply with the WFS filter
    grid_geometry_strings = grid_geometry_strings.apply(lambda x: x[10:-2])
    logging.info(
        "-> Start loop for layer %s and merging %i cells",
        layer,
        len(grid_geometry_strings),
    )

    gdf_bag_total = gpd.GeoDataFrame()
    counter = 0
    for cell in grid_geometry_strings:
        counter += 1
        logging.info(
            "--> Processing cell number %i/%i", counter, len(grid_geometry_strings)
        )
        # polygon_coordinates = '89786 460729,89786 461307,98560 461307,98560 460729,89786 460729'
        filter_str = create_ogc_polygon_filter(cell)
        params = dict(
            service="WFS",
            request="GetFeature",
            typeName=type_name,
            srsname=coordinate_sys,
            version="2.0.0",
            filter=filter_str,
        )
        gdf = gdf_wfs_appended(url, params, sortby, count)
        gdf_bag_total = gdf_bag_total.append(gdf)
    gdf_bag_total = gdf_bag_total.reset_index(drop=True)
    return gdf_bag_total


def get_grid_cells_gdf(gdf, n_cells):
    """
    Receives a geodataframe and desired number of cells that should fit (in x
    direction) in the bounds of the geodataframe.
    Calculates the xmin, xmax, ymin and ymax of each cell so that all cells
    will fit in the bounds of the given geodataframe.

    Parameters
    ----------
    gdf: geodataframe
        The given geodataframe of which the bounds are used to calculate the
        cells.
    n_cells: int
        The desired number of cells that should fit in the dataframe bounds.

    Returns
    -------
    cells: geodataframe
        A geodataframe that contains a square cell in each row with xmin, xmax, ymin, ymax.
    """
    # total area for the grid
    xmin, ymin, xmax, ymax = gdf.total_bounds
    # how many cells across and down
    cell_size = (xmax - xmin) / n_cells
    # projection of the grid
    crs = gdf.crs
    # "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181
    # +b=6371007.181 +units=m +no_defs"
    # create the cells in a loop
    grid_cells = []
    for x_0 in np.arange(xmin, xmax + cell_size, cell_size):
        for y_0 in np.arange(ymin, ymax + cell_size, cell_size):
            # bounds
            x_1 = x_0 - cell_size
            y_1 = y_0 + cell_size
            grid_cells.append(shapely.geometry.box(x_0, y_0, x_1, y_1))
    cells = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=crs)
    return cells


def gdf_wfs_appended(url, params, sortby, count):
    """Load all WFS features into one GeoDataFrame

    Parameters
    ----------
    total_features : int
        Total features in a WFS request

    Returns
    -------
    GeoDataFrame
        Containing records equal to total_features
    """
    startindex = 0
    total_features = request_count(url, params)
    page_indexes = range(startindex, total_features, int(count))
    gdf_final = gpd.GeoDataFrame()
    logging.info(
        "---> Downloading and appending %i page(s)", len(page_indexes)
    )  # nopep8
    for page_nr in page_indexes:
        logging.info("---> Downloading startindex %s", page_nr)
        gdf = request_gdf_wfs(url, params, page_nr, sortby, count)
        gdf_final = gdf_final.append(gdf)
    return gdf_final.reset_index(drop=True)


def request_count(url, params):
    """Counts the amount of features in a WFS request

    Parameters
    ----------
    params : dict
        query parameter

    Returns
    -------
    int
        amount of features returned by query
    """
    parameters = merge_dicts(params, dict(resulttype="hits"))
    req_url = requests.Request("GET", url, params=parameters).prepare().url
    resp = requests.get(req_url)
    xml_str = resp.content.decode("utf-8")
    xml = ET.fromstring(xml_str)
    numbers_matched = int(xml.attrib["numberMatched"])
    logging.info("--> request_count found %i feature in cell", numbers_matched)
    return numbers_matched


def request_gdf_wfs(url, params, startindex, sortby, count):
    """Build a WFS url and read into GeoDataFrame

    For more information on the count and startindex parameters, see:
    https://geoforum.nl/t/bag-wfs-request-geeft-maar-1000-objecten-maximaal-terug/2405/2

    Parameters
    ----------
    params : dict
        dictionary of parameters required for the WFS protocol
        Example for the bag {'service': 'WFS',
                             'request': 'GetFeature',
                             'typeName': 'verblijfsobject',
                             'srsname': 'EPSG:28992',
                             'version': '2.0.0',
                             'filter': <a geospatial query in ogc format
                                       see function create_ogc_polygon_filter>,
                             'sortby': 'identificatie',
                             'count': '1000',
                             'startindex': '2500'}
    count: str or int
        Number of features that one WFS request returns
        The WFS service has this usually capped at 1000
    startindex : str or int
        Number of the first returned feature in the request
    sortby:
        Column to sort, must be in the requested layer

    Returns
    -------
    GeoDataFrame
        containing features equal to count
    """
    if int(startindex) > 50000:
        logging.warning(
            "Most WFS services have startindexlimit of not more than 50000,\
 make your geometry filter smaller"
        )
    startindex = str(startindex)
    parameters = merge_dicts(
        params,
        dict(outputFormat="json", sortby=sortby, count=count, startindex=startindex),
    )
    req_url = requests.Request("GET", url, params=parameters).prepare().url
    gdf = gpd.read_file(req_url)
    return gdf


def create_ogc_polygon_filter(polygon_coordinates):
    """create a polygon filter for a WFS request, according to the ogc standard

    filter reference:
    https://docs.geoserver.org/latest/en/user/filter/filter_reference.html

    Parameters
    ----------
    polygon_coordinates : str
        Str with enclosed points in the coordinate system requested, eg:
        '89800 460800,89800 460900,89900 460900,89900 460800,89900 460900'

    Returns
    -------
    str
        Ogc filter as a str, in xml format
    """
    el_filter = ET.Element("ogc:Filter")
    el_not = ET.SubElement(el_filter, "Not")
    el_disjoint = ET.SubElement(el_not, "Disjoint")
    el_pn = ET.SubElement(el_disjoint, "PropertyName")
    el_pn.text = "Geometry"
    el_p = ET.SubElement(el_disjoint, "gml:Polygon")
    el_ob = ET.SubElement(el_p, "gml:outerBoundaryIs")
    el_ls = ET.SubElement(el_ob, "gml:LinearRing")
    el_pl = ET.SubElement(el_ls, "gml:posList")
    el_pl.text = polygon_coordinates
    filter_str = ET.tostring(el_filter, encoding="utf8", method="xml").decode()[38:]
    logging.debug("Ogc filter string is %s", filter_str)
    return filter_str


def merge_dicts(dict1, dict2):
    """Merge two dictionaries together"""
    dict_merged = dict1.copy()  # start with keys and values of x
    dict_merged.update(dict2)  # modifies z with keys and values of y
    return dict_merged
