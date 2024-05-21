"""Contains help functions for the check.
Placed in a separate file because otherwise the main file becomes very
long and messy.
"""
import os
import geopandas as gpd
import ogr
import shapely
import logging

import dataportaal.io.pdok_wfs as pdok_wfs


def create_test_bbox(gdf):
    """
    Receives a geodataframe and determines the bounding box of it.
    Next calculates a limited bounding box that is within the bounding
    box of the given geodataframe. This limited bounding box is meant for
    test purposes to limit the calculation time. 
    Returns the borders of this limited bounding box.
    
    Parameters
    ----------
    gdf: geodataframe
        The geodataframe of which a test bounding box will be 
        calculated.

    Returns
    -------
    xmin_test: int
        The x of the lower left corner of the test bbox
    ymin_test: int
        The y of the lower left corner of the test bbox
    xmax_test: int
        The x of the upper right corner of the test bbox
    ymax_test: int
        The y of the upper right corner of the test bbox
    gdf_test: geodataframe
        The input geodataframe limited to the test bbox
    polygon_test: shapely polygon
        The test bbox as shapely polygon
    """
    xmin = gdf.total_bounds[0]
    xmax = gdf.total_bounds[2]

    ymin = gdf.total_bounds[1]
    ymax = gdf.total_bounds[3]

    deltax = xmax - xmin
    xmin_test = int(xmin + deltax / 2 - 500)
    xmax_test = int(xmin_test + 1000)

    deltay = ymax - ymin
    ymin_test = int(ymin + deltay / 2 - 500)
    ymax_test = int(ymin_test + 1000)

    ## These test coordinates were used to solve a specific bug in some area.
    # xmin_test = 93000
    # ymin_test = 462500
    # xmax_test = 94000
    # ymax_test = 463500

    polygon_test = f"""POLYGON(({xmin_test} {ymin_test},
                                {xmin_test} {ymax_test},
                                {xmax_test} {ymax_test},
                                {xmax_test} {ymin_test},
                                {xmin_test} {ymin_test}))"""

    polygon_test= shapely.wkt.loads(polygon_test)
    gdf_test = gpd.GeoDataFrame(index=[0],
                                        crs='epsg:28992',
                                        geometry=[polygon_test])
    return xmin_test, ymin_test, xmax_test, ymax_test, \
           gdf_test, polygon_test


def merge_shapes(municipalities, data_folder, filename):
    """Merges shapefiles that are created during the proces
    (function ongerioleerde_percelen_check()) and writes the merged result to
    a subdirectory 'common_shapes' in the given output folder.

    Parameters
    ----------
    municipalities: array of str
        An array with the names of the municipalitiess for which the calculation
        should be performed.
    data_folder: str
        The directory in which all results are written, the check created for
        each municipality a subdirectory. In this folder a new folder 
        'common_shapes' is created.
    filename: str
        The name of the shape that has to be merged. The check uses standard names
        like 'buffer_total.shp' or 'plots_outside_buf.shp'.

    Returns
    -------
    None
    """
    result_path =  os.path.join(data_folder,
                                    'common_shapes')
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    merge_result = None
    for municipality in municipalities:
        filepath = os.path.join(data_folder, municipality, 'shapes', filename)
        if os.path.exists(filepath):
            gdf = gpd.read_file(filepath)
            gdf['gemeente'] = municipality
            if merge_result is None:
                merge_result = gdf
            else:
                merge_result = merge_result.append(gdf)

    # add_fields_and_save_as_shape(merge_result,
    #                              os.path.join(result_path,
    #                                           "test_merged_" + filename))
    merge_result.to_file(os.path.join(result_path, "merged_" + filename))

def add_fields_to_shape(path_to_shape):
    # driver = ogr.GetDriverByName('ESRI Shapefile')
    datasource = ogr.Open(path_to_shape, 1)  # 1 means read / write

    fldDef = ogr.FieldDefn('Omschrijvi', ogr.OFTString)
    fldDef.SetWidth(250)
    fldDef.SetDefault(' ')
    fldDef2 = ogr.FieldDefn('Maatregel', ogr.OFTString)
    fldDef2.SetWidth(25)
    fldDef.SetDefault(' ')
    fldDef3 = ogr.FieldDefn('I.E.', ogr.OFTInteger)
    fldDef3.SetWidth(8)
    fldDef.SetDefault('0')

    layer = datasource.GetLayer()
    layer.CreateField(fldDef)
    layer.CreateField(fldDef2)
    layer.CreateField(fldDef3)


def add_fields_and_save_as_shape(gdf, output_path):
    """
    Add fields to the result shape file. Until now not successfull.

    If fields are added the exported shape gives an error in ArcMap
    This function is created to explore several options to avoid the
    error. Not successfull until now.

    It seems like using shape is not a good option.
    
    """
    # gdf["Omschrijvi"] = '' 
    # gdf["Maatregel"] = ''
    # gdf["I.E."] = 0

    # Get geopandas to populate a schema dict so we don't have to build it from scratch
    schema = gpd.io.file.infer_schema(gdf)
    # print(schema)
    schema['properties']['Omschrijvi'] = 'str:250'  # 'Short integer' format
    schema['properties']['Maatregel'] = 'str:25'  # 'Long integer' format
    schema['properties']['I.E.'] = 'int32:10'  # 'Long integer' format
    print(schema)
    with open(output_path, 'w') as f:
        # f.write(gdf.to_json())
        f.write(gdf.to_file())
    # gdf.to_file(output_path, schema=schema)
    # print(fiona.drivers)
    # add_fields_to_shape(output_path)

def create_export(munlist, datapath):
    """
    Help function to create export of all results of a given list 
    of municipalities.

    Parameters
    ----------
    munlist: array of str
        The list of municipalities
    datapath: str
        The path to write all resultfiles to
    
    Returns
    -------
    None
    """
    
    outputf = os.path.join(datapath,
                           'tmp_export')
    if not os.path.exists(outputf):
        os.makedirs(outputf)
    for municipality in munlist:
        input_path = os.path.join(datapath, municipality, "shapes")
        gdf = gpd.read_file(os.path.join(input_path,
                            "plots_bag_outside_buf.shp"))
        gdf.to_file(os.path.join(outputf, f'{municipality}_result.shp'))

def get_bag_objects(municipality, datapath):
    """
    This function can be used in case the BAG service crashed (did happen few
    times). If a municipality name is given, the exported shape with plots that
    are outside the 40m buffer is loaded. After that the function
    'check_all_plots_on_bag' is called for this dataframe.
    Thus the opportunity is offered to get the final result for a 
    community in case of a crash in the end.

    The final result is stored in the output directory for the given 
    municipality.

    Parameters
    ----------
    municipality: str
        The name of the municipality of which the calculation should be done.
    datapath: str
        The basepath where the results of all municipalities are stored.

    Returns
    -------
    None
    """
    datapath = os.path.join(datapath,
                            municipality,
                            "shapes")
    gdf = gpd.read_file(os.path.join(datapath,
                        "plots_outside_buf.shp"))
    logging.info('Checking for %s, %s plots on bag, this might take some time', 
                 municipality, len(gdf.lokaalID))

    # Long names are abbreviated in a shape. Make sure to restore the
    # names that are used in the function.
    c_gemcode = 'kadastraleAanduiding|TypeKadastraleAanduiding|' + \
        'aKRKadastraleGemeenteCode|AKRKadastraleGemeenteCode|waarde'
    c_perceelnr = 'perceelnummer'

    gdf = gdf.rename(columns = {"kadastra_3": c_gemcode,"perceelnum": c_perceelnr})
    result = pdok_wfs.check_allplots_on_bag(gdf)

    result.to_file(os.path.join(
                   datapath,
                   "plots_bag_outside_buf.shp"))


def filter_plots_with_bagobject(gdf):
    """
    Receives an geodataframe with kadaster plot data. Uses function
    find_address_plot for each row in the dataframe to check if that
    plot has an address. If so, it will be selected otherwise it is removed.

    Parameters
    ----------
    gdf: geodataframe
        Contains kadaster plot data. Has at least these three fields:
        'gemeente', 'sectie', 'perceelnummer'

    Returns
    -------
    geodataframe with plots that have an BAG address.
    ATTENTION:
    In this version the addresses are not returned yet.
    """
    gdf['BAG'] = False
    for index, row in gdf.iterrows():
        test = find_address_plot(row['gemeente'],
                                 row['sectie'],
                                 row['perceelnummer'])
        if test is not None:
            gdf.loc[gdf.index==index, 'BAG'] = True
    return gdf.loc[gdf['BAG']]


def find_address_plot(gemcode, sectie, number):
    """
    Receives a plot id of a kadaster plot and uses the pdok webservice to find
    the addresses for that plot if any. Returns a dataframe with the addresses
    for the plot or None if no addresses are available

    Parameters
    ----------
    gemcode: str
        The community code (e.g. 'LDN01' for Leiden)
    sectie: str
        The letter for the kadastral section (e.g. 'R')
    number: str
        The number of the plot

    Returns
    -------
    A dataframe with the addresses that are found for the given parameters or
    None if no addresses were found.
    """

    url  = f"https://geodata.nationaalgeoregister.nl/locatieserver/v3/suggest?" +\
           f"q=gekoppeld_perceel:{gemcode}-{sectie}-{number}&" +\
           f"fl=type,weergavenaam,id,score,gekoppeld_perceel"

    retValue = None
    q = requests.get(url)
    if (q.status_code == 200):
        response = q.json()['response']
        if response['numFound'] > 0:
            retValue = pd.DataFrame(data=response['docs'])
    return retValue