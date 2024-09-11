from math import hypot

# from pydantic import BaseModel
from typing import List
from settings import (
    MUNICIPALITIES,
    MAX_PUMP_TO_SEWERLINE_DISTANCE,
    MAX_SEWER_LINE_TO_PLOT_DISTANCE,
    MAX_SEWER_LINE_TO_SEWER_LINE_DISTANCE,
)
import logging
from pathlib import Path
import geopandas as gpd

# import sys
from shapely.geometry import MultiPoint, LineString, MultiPolygon, Point

# from shapely.ops import unary_union
# from shapely import get_coordinates, intersection
import shapefile

from tqdm import tqdm
from objects import Pump, Plot, WaterDeel, SewerLine

# als FORCE_RELOAD op True staat dan worden eerder gecreeerde analyse resultaten opnieuw gegenereerd
# dit duurt (veel) langer maar kan nodig zijn als bepaalde parameters veranderd zijn, bv de
# afstand tussen pomp en leiding of de grootte van de buffer rondom de leiding etc.
FORCE_RELOAD = False


def _rec_connect_closest_sewerline(
    pump: Pump,
    sewerlines: List[SewerLine],
    x: float,
    y: float,
    max_distance: float,
):
    for i, sl in enumerate(sewerlines):
        if pump.gml_id in sl.connected_pump_ids:
            continue
        dl1 = hypot(sl.x1 - x, sl.y1 - y)
        dl2 = hypot(sl.x2 - x, sl.y2 - y)
        if min([dl1, dl2]) < max_distance:
            sewerlines[i].connected_pump_ids.append(pump.gml_id)
            pump.connected_sewer_lines.append(sewerlines[i])
            _rec_connect_closest_sewerline(
                pump, sewerlines, sl.x1, sl.y1, MAX_SEWER_LINE_TO_SEWER_LINE_DISTANCE
            )
            _rec_connect_closest_sewerline(
                pump, sewerlines, sl.x2, sl.y2, MAX_SEWER_LINE_TO_SEWER_LINE_DISTANCE
            )


for municipality in MUNICIPALITIES:
    data_folder = f"data/{municipality}"
    analyse_folder = f"analyse/{municipality}"
    pumps = []
    plots = []
    waterdelen = []
    sewerlines = []

    Path(analyse_folder).mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=f"{analyse_folder}/{municipality.lower()}_analyse.log",
        level=logging.DEBUG,
        filemode="w",
    )
    logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)

    ##########################
    # STAP 1A - pompen laden #
    ##########################
    try:
        gdf_pumps = gpd.read_parquet(Path(data_folder) / "2b_pump_points.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file '2b_pump_points.parquet' in the given path '{data_folder}'"
        )
        continue

    for _, r in gdf_pumps.iterrows():
        if type(r.geometry) == MultiPoint:
            geoms = r.geometry.geoms
        else:
            geoms = [r.geometry]

        for geom in geoms:
            x, y = [
                n[0] for n in geom.xy
            ]  # coords won't work because of Point Z types mixed with Point types, xy will return 2 numpy arrays
            gml_id = r.id
            if gml_id is None:
                logging.info(f"Found pump without gml_id at x={x:.2f}, y={y:.2f}")
            else:
                pumps.append(Pump(gml_id=gml_id, x=x, y=y))

    logging.info(f"Found {len(pumps)} pumps.")

    #########################
    # STAP 1B - plots laden #
    #########################
    try:
        gdf_plots = gpd.read_parquet(Path(data_folder) / "1_plots.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file '1_plots.parquet' in the given path '{data_folder}'"
        )
        continue

    for _, r in gdf_plots.iterrows():
        gml_id = r.gml_id
        xs, ys = r.geometry.exterior.coords.xy
        if gml_id is None:
            logging.info(f"Found plot without gml_id at {r.geometry.centroid}")
        else:
            plots.append(
                Plot(
                    gml_id=gml_id,
                    polygon=[p for p in zip(xs, ys)],
                )
            )

    logging.info(f"Found {len(plots)} plots.")

    ##############################
    # STAP 1C - waterdelen laden #
    ##############################
    try:
        gdf_waterdelen = gpd.read_parquet(Path(data_folder) / "bgt_waterdeel.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file 'bgt_waterdeels.parquet' in the given path '{data_folder}'"
        )
        continue

    for _, r in gdf_waterdelen.iterrows():
        gml_id = r.gml_id
        xs, ys = r.geometry.exterior.coords.xy

        waterdelen.append(
            WaterDeel(
                gml_id=gml_id,
                polygon=[p for p in zip(xs, ys)],
            )
        )

    # gdf_waterdelen = gpd.GeoDataFrame(
    #     geometry=[wd.shapely_polygon for wd in waterdelen], crs="EPSG:28992"
    # )

    logging.info(f"Found {len(waterdelen)} waterparts.")

    #############################
    # STAP 1D - leidingen laden #
    #############################
    try:
        gdf_sewerlines = gpd.read_parquet(Path(data_folder) / "2a_sewer_lines.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file '2a_sewer_lines.parquet' in the given path '{data_folder}'"
        )
        continue

    for _, row in gdf_sewerlines.iterrows():
        for geom in row.geometry.geoms:
            for i in range(0, len(geom.coords) - 1, 2):
                x1, y1 = geom.coords[i][0], geom.coords[i][1]
                x2, y2 = geom.coords[i + 1][0], geom.coords[i + 1][1]
                if row.id is None:
                    logging.info(f"Found pump without gml_id at x={x1:.2f}, y={y1:.2f}")
                else:
                    sewerlines.append(
                        SewerLine(
                            gml_id=row.id,
                            x1=round(x1, 2),
                            y1=round(y1, 2),
                            x2=round(x2, 2),
                            y2=round(y2, 2),
                        )
                    )

    logging.info(f"Found {len(sewerlines)} sewerlines.")

    #######################################################
    # STAP 2 - Verwijder leidingen die waterdelen snijden #
    #######################################################

    #############################################
    # STAP 2A - maak 1 object van de waterdelen #
    #############################################
    waterdelen_filename = f"{analyse_folder}/02A_waterdelen_union.parquet"
    if not Path(waterdelen_filename).exists() or FORCE_RELOAD:
        logging.info(f"Creating combined waterdelen shape.")
        union_waterdelen = gdf_waterdelen.union_all()
        gdf_waterdelen.to_parquet(waterdelen_filename)
    else:
        logging.info(f"Reading combined waterdelen shape...")
        gdf_waterdelen = gpd.read_parquet(waterdelen_filename)
        union_waterdelen = gdf_waterdelen.union_all()

    ############################################################
    # STAP 2B - verwijder leidingen die snijden met waterdelen #
    ############################################################
    filtered_sewerlines_filename = f"{analyse_folder}/02B_filtered_sewerlines.shp"

    if not Path(filtered_sewerlines_filename).exists() or FORCE_RELOAD:
        logging.info("Filtering sewerlines that intersect with the waterdelen objects.")
        filtered_sewerlines = []
        print("Filtering sewerlines that intersect with the waterdelen objects....")
        for sewerline in tqdm(sewerlines):
            if sewerline.shapely_linestring.intersects(union_waterdelen):
                logging.info(
                    f"Sewerline '{sewerline.gml_id}' intersects with a waterpart so it is removed from the analysis."
                )
            else:
                filtered_sewerlines.append(sewerline)

        # create a shape for later access so we can skip this step
        w = shapefile.Writer(filtered_sewerlines_filename)
        w.field("id")
        for sewerline in filtered_sewerlines:
            w.line(
                [
                    [
                        (sewerline.x1, sewerline.y1),
                        (sewerline.x2, sewerline.y2),
                    ]
                ]
            )
            w.record(sewerline.gml_id)
        w.close()
        sewerlines = [o for o in filtered_sewerlines]
    else:
        logging.info(
            f"Reading sewerlines that do not intersect with the waterdelen objects."
        )
        sewerlines = []
        try:
            gdf_sewerlines_filtered = gpd.read_file(filtered_sewerlines_filename)
        except Exception as e:
            logging.error(
                f"Cannot find the file '{filtered_sewerlines_filename}' in the given path '{data_folder}'"
            )
            continue

        # reload the sewerlines, this time
        for _, row in gdf_sewerlines_filtered.iterrows():
            geometry = row["geometry"]
            coords = list(geometry.coords)
            sewerlines.append(
                SewerLine(
                    gml_id=row["id"],
                    x1=coords[0][0],
                    y1=coords[0][1],
                    x2=coords[1][0],
                    y2=coords[1][1],
                )
            )

    logging.info(f"Found {len(sewerlines)} sewerlines not crossing any waterparts.")

    ########################################
    # STAP 3 - Koppelen leiding met pompen #
    ########################################
    connected_sewerlines_filename = f"{analyse_folder}/03_connected_sewerlines.shp"
    if not Path(connected_sewerlines_filename).exists() or FORCE_RELOAD:
        w = shapefile.Writer(connected_sewerlines_filename)
        w.field("id")
        w.field("pump_id")

        print("Connecting pumps and sewerlines...")
        for pump in tqdm(pumps):
            _rec_connect_closest_sewerline(
                pump, sewerlines, pump.x, pump.y, MAX_PUMP_TO_SEWERLINE_DISTANCE
            )
            if len(pump.connected_sewer_lines) > 0:
                Path(f"{analyse_folder}/debug").mkdir(parents=True, exist_ok=True)
                sewerline_filename = f"{analyse_folder}/debug/03_pump_{pump.gml_id}.shp"
                w_pump = shapefile.Writer(sewerline_filename)
                w_pump.field("sewerline_id")
                for sl in pump.connected_sewer_lines:
                    w_pump.line([[(sl.x1, sl.y1), (sl.x2, sl.y2)]])
                    w_pump.record(sl.gml_id)
                w_pump.close()

            for sl in pump.connected_sewer_lines:
                w.line([[(sl.x1, sl.y1), (sl.x2, sl.y2)]])
                w.record(sl.gml_id, pump.gml_id)
        w.close()
    else:
        print("Reading connected sewerlines...")
        connected_sewerlines = []
        try:
            gdf_connected_sewerlines = gpd.read_file(connected_sewerlines_filename)
        except Exception as e:
            logging.error(
                f"Cannot find the file '{connected_sewerlines_filename}' in the given path '{data_folder}'"
            )
            continue

        # reload the sewerlines, this time with the connections
        for _, row in tqdm(gdf_connected_sewerlines.iterrows()):
            for sl in sewerlines:
                if sl.gml_id == row["id"]:
                    sl.connected_pump_ids.append(row["pump_id"])
                    for pump in pumps:
                        if pump.gml_id == row["pump_id"]:
                            pump.connected_sewer_lines.append(sl.gml_id)
                            break
                    break

    ##
    connected_sewerlines = [
        sewerline for sewerline in sewerlines if sewerline.connected
    ]
    logging.info(f"Found {len(connected_sewerlines)} connected sewerlines.")

    #################################################
    # STAP 4 - Per pomp aangesloten percelen vinden #
    #################################################
    logging.info(f"Creating a buffer around the pumps and connected sewerlines.")

    sewerline_data = {"id": [], "geometry": []}
    for sewerline in connected_sewerlines:
        ls = LineString([[sewerline.x1, sewerline.y1], [sewerline.x2, sewerline.y2]])
        buffer = ls.buffer(MAX_SEWER_LINE_TO_PLOT_DISTANCE)
        sewerline_data["id"].append(sewerline.gml_id)
        sewerline_data["geometry"].append(buffer)

    gdf = gpd.GeoDataFrame(sewerline_data, crs="EPSG:28992")
    gdf.to_file(Path(analyse_folder) / "04A_sewerlines_with_buffer.shp")

    # now combine the geometries of the pump sewerlines with the buffer to one polygon with
    # the pump id as the field
    dissolved_gdf = gdf.dissolve()
    dissolved_gdf.to_file(
        Path(analyse_folder) / "04B_sewerlines_with_buffer_dissolved.shp"
    )

    # subtract the water of the buffered lines
    dissolved_water_gdf = gdf_waterdelen.dissolve()
    dissolved_water_gdf.to_file(Path(analyse_folder) / "04C_waterdelen_combined.shp")

    subtracted_water_gdf = gpd.overlay(
        dissolved_gdf, gdf_waterdelen, how="difference", keep_geom_type=True
    )
    subtracted_water_gdf.to_file(
        Path(analyse_folder) / "04D_sewerlines_with_buffer_water_subtracted.shp"
    )

    # create individual polygons out of this dataframe
    individual_polygons = []
    for _, row in subtracted_water_gdf.iterrows():
        geom = row["geometry"]
        if isinstance(geom, MultiPolygon):
            for poly in geom.geoms:
                individual_polygons.append({"geometry": poly})
        else:
            individual_polygons.append({"geometry": geom})

    gdf_sewerlines = gpd.GeoDataFrame(individual_polygons, crs="EPSG:28992")
    gdf_sewerlines.to_file(
        Path(analyse_folder) / "04E_individual_sewerline_polygons.shp"
    )

    # remove all the polygons that do not contain an active pump
    active_pumps = [p for p in pumps if p.num_connected_sewer_lines > 0]

    # create a geodataframe for the active pumps
    gdf_active_pumps = gpd.GeoDataFrame(
        {"geometry": [Point(p.x, p.y) for p in active_pumps]},
        index=[p.gml_id for p in active_pumps],
        crs="EPSG:28992",
    )
    gdf_active_pumps.to_file(Path(analyse_folder) / "04F_active_pumps.shp")

    gdf_join = gpd.sjoin(
        gdf_sewerlines, gdf_active_pumps, how="inner", predicate="contains"
    )
    gdf_join.to_file(Path(analyse_folder) / "04G_active_areas.shp")

    # now find all plots that intersect with these active areas
    gdf_join = gdf_join.drop(["index_right"], axis=1)
    gdf_result = gpd.sjoin(gdf_plots, gdf_join, how="inner")

    gdf_result.to_file(Path(analyse_folder) / "05_plots_with_sewer.shp")
