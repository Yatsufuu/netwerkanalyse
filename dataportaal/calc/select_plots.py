"""This file contains functions that are used to do
the 'calculations' which are in fact selections of plots.
"""
import logging
import geopandas as gpd


def create_buffer(gdf, crs="EPSG:28992", bufsize=40):
    """
    Receives a geodataframe, sets the projection and then
    creates a buffer around its geometry. All buffered geometries are
    unified through a union and returned as one polygon.

    Parameters
    ----------
    gdf: geodataframe
        The geodataframe for which a buffer is created.
    crs: crs, optional
        The desired crs for the input dataframe.
        Default is 'EPSG:28992'
    bufsize: float, optional
        The buffersize. Default is 40

    Returns
    -------
    overlay: Polygon
        A polygon of the created buffer

    """
    gdf_buf = gdf.copy(deep=True)
    gdf_buf = gdf_buf.to_crs(crs)

    gdf_buf['geometry'] = gdf_buf.buffer(bufsize)
    return gdf_buf.unary_union


def select_items_with_spatial_index(source_gdf, mask_gdf):
    """
    Used to limit the number of rows in a geodataframe for an 'expensive'
    operations using the spatial index. In this way unary_union etc. will be
    limited to less rows so that the calculation will be much quicker.

    Parameters
    ----------
    source_gdf: geodataframe
        The geodataframe of which the records should be limited.
    mask_gdf: geodataframe
        The geodataframe that will be used as the mask to limit the
        number of records.

    Returns
    -------
    candidate_gdf: geodataframe
        The geodataframe with the candidates of which the spatial index
        intersects with the mask geometry.
    exclude_gdf: geodataframe
        The geodataframe with the excluded rows that do not intersect
        with the mask geometry.
    """
    matches = []
    source_sindex = source_gdf.sindex
    for values in mask_gdf.bounds.values:
        bounds = list(values)
        matches += list(source_sindex.intersection(bounds))

    unique_matches = list(set(matches))

    candidate_gdf = source_gdf.loc[unique_matches]
    exclude_gdf = source_gdf[~source_gdf.index.isin(candidate_gdf.index)]
    return candidate_gdf, exclude_gdf


def clip_buffer_select_part(gdf_buffer, mask, gdf_lines):
    """
    Clips a (buffer) polygon dataframe with a given (mask) polygon dataframe
    (e.g. water or railway).
    Next uses a line dataframe (e.g. sewers) to remove the (buffer) polygons
    that do not intersect with a line.

    Parameters
    ----------
    gdf_buffer: geodataframe
        The dataframe with the buffer polygons that should be cut.
    mask: geodatframe
        The dataframe with the mask that will be used to clip the
        buffer polygons.
    gdf_lines: geodataframe
        The dataframe with the (sewer) lines to select the polygons
        that intersect with a line.

    Returns
    -------
    gdf_cb: geodataframe
        Dataframe with clipped buffer polygons that intersect with
        a sewer line.
    """
    logging.info("Clipping buffergeometry.")
    # mask_poly = gdf_mask.geometry.unary_union
    # clipped_buffer = gdf_buffer.cx(mask_poly)
    clipped_buffer = gdf_buffer.difference(mask, align=False)

    # result = result.loc[~result.contains(gdf_leidingen)]
    gdf_cb = gpd.GeoDataFrame()
    gdf_cb = gdf_cb.set_geometry(clipped_buffer)
    gdf_cb = gdf_cb.set_crs("EPSG:28992")

    try:
        logging.info("Explode clipped buffergeometry.")
        gdf_cb = gdf_cb.explode(column="geometry", ignore_index=True)
        # gdf2 = gdf_cb.reset_index().rename(columns={0: 'geometry'})
        # gdf_cb = gdf2.merge(gdf_cb.drop('geometry',
        #                     axis=1),
        #                     left_on='level_0',
        #                     right_index=True)
        # gdf_cb = gdf_cb.set_index(['level_0',
        #                            'level_1']).set_geometry('geometry')
    except Exception as e:
        logging.error("Explode in clipped buffer error: {%s}", e)

    logging.info("Skip buffers without sewer/pump.")
    # 'within' 'intersects' 'contain'
    gdf_cb = gpd.sjoin(gdf_cb, gdf_lines, 'left', 'intersects')
    gdf_cb = gdf_cb.loc[~gdf_cb['index_right'].isna()]
    return gdf_cb


def select_plots(gdf, overlay):
    """
    Selects plots of a geodataframe with plots from kadaster that are outside a
    given overlay-geometry.
    The selection is performed through the 'intersect' method.

    Parameters
    ----------
    gdf: geodataframe
        The dataframe with the plot geometries.
    overlay: Polygon
        The geometry used to perform the selection.

    Returns
    -------
    gdf: geodataframe
        The geodataframe with the final result
    """
    gdf['intersect'] = gdf.geometry.apply(lambda x: x.intersects(overlay))
    gdf_outside = gdf.loc[~gdf['intersect']]
    # gdf_outside = gdf[gdf.disjoint(overlay)] # Dit doet precies hetzelfde
    return gdf_outside


def select_polygons(gdf_poly, gdf_object):
    """
    Returns array of index values for polygons that contain an
    other (e.g. bag) object. Returns a dataframe that contains an
    object and a dataframe dat does not contain an object.
    The second dataframe can be used for further selection.

    Parameters
    ----------
    gdf_poly: geodataframe
        Dataframe with the polygons to select

    gdf_object: geodataframe
        Dataframe with the (bag) objects used for the selection

    Returns
    -------
    gdf_selected: geodataframe
        Geodataframe with plots that do contain one or more object(s).
    gdf_not_selected: geodataframe
        Geodataframe with rest of the plots. This can be used
        for further selection with other objects.
    """
    select_index = []
    for i, row in gdf_poly.iterrows():
        if (gdf_object.within(row.geometry.convex_hull).any()):
            #          print ("plot met locaties erin: {}".format(i))
            select_index.append(i)

    gdf_selected = gdf_poly.loc[gdf_poly.index.isin(select_index)]
    gdf_not_selected = gdf_poly.loc[~gdf_poly.index.isin(select_index)]
    return gdf_selected, gdf_not_selected
