import folium
from folium.plugins import MarkerCluster
from folium.raster_layers import WmsTileLayer
import geopandas as gpd


class FoliumMap():
    """
    Implentation of folium map.

    Parameters
    ----------
    center_map: Array(2) of float, optional
        Contains the latitude and longitude to center the map on.
        (Northing, Easting) in WGS84
        Default is [52.16, 4.5]

    """
    wms_layers = []
    # wfs_layers = []
    center_map = None

    def __init__(self, center_map=[52.16, 4.5]):
        self.map = folium.Map(location=center_map, zoom_start=13)

    def add_wms_layer(self, url, name, visible=False):
        """
        Creates an instance of WmsTileLayer with given parameters and adds it to
        the internal list of wms services. They will be added to the map after 
        calling init_map.
        This has been implemented this way to make sure that the wms services are
        added after the data layers.

        Parameters
        ----------
        url: str
            The base url of the WMS service
        name: str
            The name of the layer to add

        Returns
        -------
        None

        """
        wms_layer = WmsTileLayer(url,
                                 name,
                                 styles='',
                                 fmt='image/png',
                                 transparent=True,
                                 version='1.1.1',
                                 name="{} - WMS".format(name),
                                 overlay=True,
                                 control=True,
                                 show=visible).add_to(self.map)
        # Add the layer to the internal array.
        # Will be added to the map at init_map.
        # self.wms_layers.append(wms_layer)

    def add_data_layer(self, geodata, name, fill_color=None,
                       line_color='black', geom="geometry"):
        """
        Adds a layer to the map with data of a given dataframe. Automatically
        detects if the layer contains points or other geometries. In case of
        points a Markercluster is added, for other geometries a
        folium.Chloropleth is added.

        Parameters
        ----------
        geodata: geodataframe
            The geodataframe with the geometries to show.
        name: str
            The name of the layer to add.
        fill_color: str, optional
            The fill color (only for polygons)
            Default is None
        line_color: str, optional
            The line color (used for markers too).
            Default is black
        geom: str, optional
            The name of the column with the geometry to show
            Default is 'black'

        Returns
        -------
        None

        """
        geometry = geodata[geom]
        # Make sure alle geodata are in WGS84 projection.

        # Onderstaande code werkt niet, snap niet waarom niet.
        # geometry = None
        # if isinstance(geodata, gpd.GeoDataFrame):
        #     geometry = geodata[geom]
        # else:
        #     geometry = geodata
        
        # Make sure the crs is wgs84!
        geometry = geometry.to_crs("epsg:4326")

        # Check if geodataframe contains points.
        if (geometry.iloc[0].geom_type == 'Point'):
            # Point geometry, initialize a new markercluster.
            marker_cluster = MarkerCluster(
                name=name,
                overlay=True,
                control=True,
                disableClusteringAtZoom=18
            )
            # icon = folium.Icon(color=line_color, icon="angry")
            # Add a circle marker to the markercluster for each point in the
            # data.
            for point in geodata[geom]:
                marker = folium.CircleMarker(
                    [point.y, point.x], radius=2, color=line_color)
                marker.add_to(marker_cluster)
            # Add the markercluster to the internal map.
            marker_cluster.add_to(self.map)
        else:
            # Add a chloropleth to the map for the line or polygon geometry and
            # add it to the internal map.
            folium.Choropleth(
                geometry,
                fill_color=fill_color,
                line_weight=1,
                line_color=line_color,
                fill_opacity=0.5,
                line_opacity=0.5,
                name=name
            ).add_to(self.map)

    def init_map(self):
        """
        Creates an interactive folium map that shows data for a given
        dataframe.

        Parameters
        ----------
        None

        Returns
        -------
        gwsw_kaart: folium.Map
            The map that shows the GWSW data given together with the WMS
            services.
        """
        # Add the wms layers to the map.
        # Somehow it was a problem to add wms layers first.
        # Therfore this solution to add them afterwards.
        for layer in self.wms_layers:
            layer.add_to(self.map)

        self.wms_layers = []

        folium.LayerControl().add_to(self.map)
        return self.map

    def ___find_center(self, gdf):
        """
        Returns the center of a given dataframe. (Not used yet.)

        Parameters
        ----------
        gdf: geodataframe
            The geometries to find the center for.

        Returns
        -------
        centroid: Array(2) of float
            A tuple with latitude and longitude in WGS84.
        """
        bbox = gdf.total_bounds
        x_c = (bbox[0] + bbox[2])/2
        y_c = (bbox[1] + bbox[3])/2
        return [y_c, x_c]

    def save_to_html(self, path):
        """
        Saves the internal folium map to html.

        Parameters
        ----------
        path: str
            The path to save the file to.

        Returns
        -------
            None
        """
        self.map.save(outfile=path)
