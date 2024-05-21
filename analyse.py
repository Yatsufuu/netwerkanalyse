from math import hypot
from pydantic import BaseModel
from typing import Tuple, List
from settings import (
    MUNICIPALITIES,
    MAX_PUMP_TO_SEWERLINE_DISTANCE,
    MAX_SEWER_LINE_TO_PLOT_DISTANCE,
    MAX_SEWER_LINE_TO_SEWER_LINE_DISTANCE,
)
from pathlib import Path
import logging
import sys
import geopandas as gpd
from shapely.geometry import MultiPoint, LineString, MultiLineString, Point
from shapely.ops import unary_union, cascaded_union
from shapely import get_coordinates, intersection
import shapefile
from tqdm import tqdm

# als FORCE_RELOAD op True staat dan worden eerder gecreeerde analyse resultaten opnieuw gegenereerd
# dit duurt (veel) langer maar kan nodig zijn als bepaalde parameters veranderd zijn, bv de
# afstand tussen pomp en leiding of de grootte van de buffer rondom de leiding etc.
FORCE_RELOAD = False


class Plot(BaseModel):
    gml_id: str
    polygon: List[Tuple[float, float]] = []

    @property
    def shapely_polygon(self):
        return Polygon(self.polygon)


class SewerLine(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    gml_id: str
    x1: float
    y1: float
    x2: float
    y2: float
    connected_pump_ids: List[str] = []

    @property
    def shapely_linestring(self):
        return LineString([(self.x1, self.y1), (self.x2, self.y2)])


class Pump(BaseModel):
    gml_id: str
    x: float
    y: float

    connected_sewer_lines: List[SewerLine] = []
    connected_plots: List[Plot] = []

    @property
    def num_connected_plots(self):
        return len(self.connected_plots)

    @property
    def num_connected_sewer_lines(self):
        return len(self.connected_sewer_lines)


def rec_connect_closest_sewerline(
    sewerlines: List[SewerLine],
    pump: Pump,
    x: float,
    y: float,
    max_distance: float,
):
    for i, sl in enumerate(sewerlines):
        dl1 = hypot(sl.x1 - x, sl.y1 - y)
        dl2 = hypot(sl.x2 - x, sl.y2 - y)

        if pump.gml_id in sewerlines[i].connected_pump_ids:
            continue

        if min(dl1, dl2) < max_distance:
            pump.connected_sewer_lines.append(sewerlines[i])
            sewerlines[i].connected_pump_ids.append(pump.gml_id)
            rec_connect_closest_sewerline(
                sewerlines, pump, sl.x1, sl.y1, MAX_SEWER_LINE_TO_SEWER_LINE_DISTANCE
            )
            rec_connect_closest_sewerline(
                sewerlines, pump, sl.x2, sl.y2, MAX_SEWER_LINE_TO_SEWER_LINE_DISTANCE
            )


for municipality in MUNICIPALITIES:
    print(f"Loading and analyzing data for '{municipality}'")
    data_folder = f"data/{municipality}"
    analyse_folder = f"data/{municipality}/analyse"

    plots = []
    pumps = []
    sewerlines = []

    # create the data folder if it is not already done
    Path(f"./{analyse_folder}").mkdir(parents=True, exist_ok=True)

    # data inladen
    # get the plots
    try:
        gdf = gpd.read_parquet(Path(data_folder) / "1_plots.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file '1_plots.parquet' in the given path '{data_folder}'"
        )
        sys.exit(1)

    for _, r in gdf.iterrows():
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

    # get the pumps
    try:
        gdf = gpd.read_parquet(Path(data_folder) / "2b_pump_points.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file '2b_pump_points.parquet' in the given path '{data_folder}'"
        )
        sys.exit(1)

    for _, r in gdf.iterrows():
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

    logging.info(f"Found {len(pumps)} pumps...")
    logging.info(f"Found {len(plots)} plots...")

    try:
        gdf = gpd.read_parquet(Path(data_folder) / "bgt_waterdeel.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file 'bgt_waterdeel.parquet' in the given path '{data_folder}'"
        )
        sys.exit(1)

    # maak er een object van
    union_waterdelen = gdf.unary_union

    # haal sewerlines weg die snijden met het spoor of met het water
    try:
        gdf = gpd.read_parquet(Path(data_folder) / "2a_sewer_lines.parquet")
    except Exception as e:
        logging.error(
            f"Cannot find the file '2a_sewer_lines.parquet' in the given path '{data_folder}'"
        )
        sys.exit(1)

    # TODO check of we de data al hebben
    check_waterdeel_intersections = True
    if not FORCE_RELOAD:
        p = Path(analyse_folder) / "01_sewerlines_cut_with_waterlines.shp"
        if p.exists():
            gdf = gpd.read_file(p)
            for _, r in tqdm(gdf.iterrows()):
                x1, y1 = r.geometry.coords[0][0], r.geometry.coords[0][1]
                x2, y2 = r.geometry.coords[1][0], r.geometry.coords[1][1]
                gml_id = r.id
                sewerline = SewerLine(
                    gml_id=gml_id,
                    x1=round(x1, 2),
                    y1=round(y1, 2),
                    x2=round(x2, 2),
                    y2=round(y2, 2),
                )
                sewerlines.append(sewerline)

            check_waterdeel_intersections = False

    if check_waterdeel_intersections:
        print("removing sewerlines that cross the waterparts...")
        w = shapefile.Writer(
            str(Path(analyse_folder) / "01_sewerlines_cut_with_waterlines.shp")
        )
        w.field("id")

        for _, r in tqdm(gdf.iterrows()):
            # the current WFS returns each sewer line as a linestring but
            # sometimes a linestring has more than two points
            for geom in r.geometry.geoms:
                for i in range(0, len(geom.coords) - 1, 2):
                    x1, y1 = geom.coords[i][0], geom.coords[i][1]
                    x2, y2 = geom.coords[i + 1][0], geom.coords[i + 1][1]
                    gml_id = r.id
                    if gml_id is None:
                        logging.info(
                            f"Found pump without gml_id at x={x1:.2f}, y={y1:.2f}"
                        )
                    else:
                        sewerline = SewerLine(
                            gml_id=gml_id,
                            x1=round(x1, 2),
                            y1=round(y1, 2),
                            x2=round(x2, 2),
                            y2=round(y2, 2),
                        )

                        # check if we intersect with a waterline, in that case do no add the line
                        # we only need to check this if the source is 2a_sewer_lines.parquet
                        # if the source is 01_sewerlines_cut_with_waterlines.shp we already
                        # analysed this data
                        if check_waterdeel_intersections:
                            if sewerline.shapely_linestring.intersects(
                                union_waterdelen
                            ):
                                logging.info(
                                    f"Sewerline '{gml_id}' intersects with a waterpart so it is removed from the analysis."
                                )
                                sewerline = None

                        if sewerline is not None:
                            sewerlines.append(sewerline)

                            if check_waterdeel_intersections:
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

    # per pomp
    # zoek alle leidingen, let op dat een leiding onderdeel van meerdere pompen kan zijn
    # doe dit alleen als het nog niet gebeurd is
    # VRAAG onderscheid per type stelsel?
    p_csv = Path(analyse_folder) / "02_connected_sewerlines.csv"
    if p_csv.exists() and not FORCE_RELOAD:
        lines = open(p_csv, "r").readlines()[1:]
        for line in lines:
            args = [a.strip() for a in line.split(",")]

            for i in range(len(pumps)):
                if pumps[i].gml_id == args[1]:
                    for j in range(len(sewerlines)):
                        if sewerlines[j].gml_id == args[0]:
                            pumps[i].connected_sewer_lines.append(sewerlines[j])
                            sewerlines[j].connected_pump_ids.append(pumps[i].gml_id)
                            break

    else:
        w = shapefile.Writer(str(Path(analyse_folder) / "02_connected_sewerlines.shp"))
        w.field("id")
        w.field("pump_id")
        print("Connecting pumps and sewerlines...")
        for pump in tqdm(pumps):
            rec_connect_closest_sewerline(
                sewerlines, pump, pump.x, pump.y, MAX_PUMP_TO_SEWERLINE_DISTANCE
            )
            # log de info
            logging.info(
                f"Pump '{pump.gml_id}' has {len(pump.connected_sewer_lines)} connected sewerlines"
            )
            for sl in pump.connected_sewer_lines:
                w.line([[(sl.x1, sl.y1), (sl.x2, sl.y2)]])
                w.record(sl.gml_id, pump.gml_id)
        w.close()

        # stap 3 - bewaar de data voor later gebruik
        with open(p_csv, "w") as f:
            f.write("sewerline_gml_id,pump_gml_id\n")
            for sl in sewerlines:
                for pump_id in sl.connected_pump_ids:
                    f.write(f"{sl.gml_id},{pump_id}\n")

    # maak een polygoon met een buffer rondom de leidingen
    # VRAAG hoe om te gaan met water

    sewerline_data = {"pump_gmlid": [], "geometry": []}

    for pump in tqdm(pumps):
        # list with all sewerlines polygons of the added buffer
        logging.info(f"Handling pump '{pump.gml_id}'")
        if pump.num_connected_sewer_lines == 0:
            logging.info("This pump has no connected sewer lines.")
            continue

        # maak een polygon rondom de leidingen die de buffer aangeeft
        logging.info(
            f"Converting the {pump.num_connected_sewer_lines} connected sewerline(s) to a polygon with a buffer."
        )

        for sewerline in pump.connected_sewer_lines:
            ls = LineString(
                [[sewerline.x1, sewerline.y1], [sewerline.x2, sewerline.y2]]
            )
            buffer = ls.buffer(MAX_SEWER_LINE_TO_PLOT_DISTANCE)
            sewerline_data["pump_gmlid"].append(pump.gml_id)
            sewerline_data["geometry"].append(buffer)

    gdf = gpd.GeoDataFrame(sewerline_data, crs="EPSG:28992")
    gdf.to_file(str(Path(analyse_folder) / "03_sewerlines_with_buffer.shp"))

    # now combine the geometries of the pump sewerlines with the buffer to one polygon with
    # the pump id as the field
    dissolved_gdf = gdf.dissolve(by=["pump_gmlid"])
    dissolved_gdf.to_file(
        str(Path(analyse_folder) / "04_sewerlines_with_buffer_combined.shp")
    )

    # snij met de plots polygonen en kijk waar intersections zijn
    # bewaar deze als geconnecte plots aan de pomp
    break
