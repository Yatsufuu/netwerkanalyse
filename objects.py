from pydantic import BaseModel
from typing import List, Tuple
from shapely.geometry import LineString, Polygon


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
    def connected(self) -> bool:
        return len(self.connected_pump_ids) > 0

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


class WaterDeel(BaseModel):
    gml_id: str
    polygon: List[Tuple[float, float]] = []

    @property
    def shapely_polygon(self):
        return Polygon(self.polygon)
