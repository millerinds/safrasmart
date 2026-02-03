"""Flask port of NDVI Viewer."""
from __future__ import annotations

import io
import json
import os
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from datetime import date, datetime, timedelta
from typing import List, Tuple, Optional, Dict, Any
import threading
import uuid

import ee
from ee import oauth
from flask import Flask, render_template, request
from google.oauth2 import service_account
import geopandas as gpd
import pandas as pd
import fiona
from shapely.geometry import shape, box

app = Flask(__name__)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Simplification tolerance used only for the overlay geometry (UI), not for NDVI math.
CEM_OVERLAY_SIMPLIFY_TOLERANCE = float(os.getenv("CEM_OVERLAY_SIMPLIFY_TOLERANCE", "0.0001"))

CEM_BBOX_GEOJSON = os.getenv(
    "CEM_BBOX_GEOJSON",
    os.path.join(BASE_DIR, "out", "CEM_bbox_4326.geojson"),
)
CEM_TALHOES_GEOJSON = os.getenv(
    "CEM_TALHOES_GEOJSON",
    os.path.join(BASE_DIR, "out", "CEM_talhoes_simplified_4326.geojson"),
)

# In-memory job registry (ephemeral, per-process).
JOBS: Dict[str, Dict[str, Any]] = {}
# In-memory imports registry (ephemeral, per-process).
IMPORTS: Dict[str, Dict[str, Any]] = {}

REGION_DEFAULT_CROP = {
    "Norte": "Pastagem",
    "Nordeste": "Cana-de-açucar",
    "Centro-Oeste": "Soja/Milho",
    "Sudeste": "Cafe",
    "Sul": "Soja/Milho",
}

CROP_THRESHOLDS = {
    "Cafe": (0.45, 0.60),
    "Cana-de-açucar": (0.55, 0.70),
    "Soja/Milho": (0.50, 0.65),
    "Pastagem": (0.35, 0.55),
    "Citros": (0.45, 0.60),
}

CROP_OPTIONS = ["Cafe", "Cana-de-açucar", "Soja/Milho", "Pastagem", "Citros"]
# When clipping by many polygons, split in batches to reduce EE complexity.
CEM_CLIP_BATCH_THRESHOLD = int(os.getenv("CEM_CLIP_BATCH_THRESHOLD", "300"))


# ---------------------------
# Earth Engine init
# ---------------------------
_EE_READY = False


def ee_initialize() -> None:
    """Initialize Earth Engine once per process (supports Fly secrets)."""
    global _EE_READY
    if _EE_READY:
        return

    json_key = os.getenv("EE_SERVICE_ACCOUNT_JSON")
    json_file = os.getenv("EE_SERVICE_ACCOUNT_FILE")
    ee_project = os.getenv("EE_PROJECT")  # e.g. "earth-484311"

    try:
        if json_key:
            service_account_info = json.loads(json_key)
            if "client_email" not in service_account_info:
                raise ValueError("Service account email address missing in json key")

            creds = service_account.Credentials.from_service_account_info(
                service_account_info,
                scopes=oauth.SCOPES,
            )

            # Prefer passing the project explicitly (more reliable on Fly)
            if ee_project:
                ee.Initialize(creds, project=ee_project)
            else:
                ee.Initialize(creds)

        elif json_file:
            with open(json_file, "r", encoding="utf-8") as fh:
                service_account_info = json.load(fh)

            if "client_email" not in service_account_info:
                raise ValueError("Service account email address missing in json file")

            creds = service_account.Credentials.from_service_account_info(
                service_account_info,
                scopes=oauth.SCOPES,
            )

            if ee_project:
                ee.Initialize(creds, project=ee_project)
            else:
                ee.Initialize(creds)

        else:
            # Local dev / ADC fallback
            if ee_project:
                ee.Initialize(project=ee_project)
            else:
                ee.Initialize()

        _EE_READY = True

    except Exception:
        _EE_READY = False
        raise


def sat_collection(cloud_rate: int, initial_date: str, updated_date: str, aoi):
    # Build the Sentinel-2 collection, pre-filtered by cloudiness and AOI.
    try:
        collection = (
            ee.ImageCollection("COPERNICUS/S2_SR")
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloud_rate))
            .filterDate(initial_date, updated_date)
            .filterBounds(aoi)
        )

        collection_size = collection.size().getInfo()
        if collection_size == 0:
            return None, (
                "No Sentinel-2 L2A images found for the selected parameters. "
                "Try increasing the cloud threshold or expanding the date range."
            )

        def clip_collection(image):
            # Clip each image to AOI and scale reflectance.
            return image.clip(aoi).divide(10000)

        collection = collection.map(clip_collection)
        return collection, None

    except Exception as exc:
        return None, f"Earth Engine encountered an error while building Sentinel-2 Collection: {exc}"


# ---------------------------
# File parsers
# ---------------------------

COLUMN_SYNONYMS = {
    "x": ["x", "ln", "lon", "lons", "lng", "lngs", "longitude", "longitudes"],
    "y": ["y", "lt", "lat", "lats", "latitude", "latitudes"],
}


def _find_column(df, possible_col_name):
    # Find first matching column name for latitude/longitude heuristics.
    for c in possible_col_name:
        if c in df.columns:
            return c
    return None


def parse_geopackage(gpkg_path: str):
    # Parse a GeoPackage into EE geometries.
    geometry_list = []
    try:
        with fiona.open(gpkg_path) as fu:
            if len(fu) == 0:
                return [], "GeoPackage contains no features."

            for feat in fu:
                if not feat or not feat.get("geometry"):
                    continue
                try:
                    geom = shape(feat["geometry"])
                except Exception:
                    continue

                if geom.geom_type == "Polygon":
                    try:
                        geometry_list.append(ee.Geometry.Polygon(list(geom.exterior.coords)))
                    except Exception:
                        continue
                elif geom.geom_type == "MultiPolygon":
                    try:
                        coords = [list(poly.exterior.coords) for poly in geom.geoms]
                        geometry_list.append(ee.Geometry.MultiPolygon(coords))
                    except Exception:
                        continue

        if not geometry_list:
            return [], "No valid Polygon or MultiPolygon geometries found in GeoPackage."
        return geometry_list, None
    except fiona.errors.DriverError:
        return [], "Invalid or corrupted GeoPackage file."
    except Exception as exc:
        return [], f"Error processing GeoPackage file: {exc}"


def parse_csv(file_storage):
    # Parse CSV into polygons (either JSON coordinates or grouped vertices).
    try:
        try:
            df = pd.read_csv(file_storage)
        except Exception:
            return [], "Invalid CSV file. Could not be read."

        if df.empty:
            return [], "CSV file is empty."

        df.columns = df.columns.str.lower().str.strip()
        geometry_list = []

        if "coordinates" in df.columns:
            for _, row in df.iterrows():
                try:
                    coords = json.loads(row["coordinates"])
                    geometry_list.append(ee.Geometry.Polygon(coords))
                except Exception:
                    return [], "Invalid 'coordinates' column. Must contain valid JSON polygon coordinates."
            if not geometry_list:
                return [], "No valid geometries found in 'coordinates' column."
            return geometry_list, None

        required_cols = ["id", "vertex_index"]
        for col in required_cols:
            if col not in df.columns:
                return [], f"Missing required column '{col}'."

        x_col = _find_column(df, COLUMN_SYNONYMS["x"])
        y_col = _find_column(df, COLUMN_SYNONYMS["y"])

        if not x_col or not y_col:
            return [], "Could not detect longitude/latitude columns."

        for gid, group in df.groupby("id"):
            try:
                group = group.sort_values("vertex_index")
                coords = group[[x_col, y_col]].values.tolist()
                if len(coords) < 3:
                    return [], f"Polygon with id '{gid}' has fewer than 3 vertices."
                if coords[0] != coords[-1]:
                    coords.append(coords[0])
                geometry_list.append(ee.Geometry.Polygon(coords))
            except Exception:
                return [], f"Invalid polygon geometry for id '{gid}'."

        if not geometry_list:
            return [], "No valid polygon geometries could be constructed from CSV."
        return geometry_list, None

    except Exception as exc:
        return [], f"Error processing CSV file: {exc}"


def parse_zip_shapefile(file_storage):
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "uploaded.zip")
            try:
                with open(zip_path, "wb") as fh:
                    fh.write(file_storage.read())
            except Exception:
                return [], "Couldn't read uploaded Zip file."

            try:
                with zipfile.ZipFile(zip_path, "r") as zf:
                    zf.extractall(tmpdir)
            except zipfile.BadZipFile:
                return [], "Not a valid zip file."

            shp_files = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.lower().endswith(".shp")]
            if not shp_files:
                return [], "Zip file does not contain .SHP file."

            try:
                gdf = gpd.read_file(shp_files[0])
            except Exception as exc:
                return [], (
                    "Failed to read Shapefile with GeoPandas: Missing one or multiple companion files "
                    f"(.SHX, .DBF, .CPG, .PRJ) {exc}"
                )

            if gdf.empty:
                return [], "Shapefile contains no features."

            geometry_list = []
            for geom in gdf.geometry:
                try:
                    if geom.geom_type == "Polygon":
                        coords = [list(geom.exterior.coords)]
                        ee_geom = ee.Geometry.Polygon(coords)
                    elif geom.geom_type == "MultiPolygon":
                        coords = [list(p.exterior.coords) for p in geom.geoms]
                        ee_geom = ee.Geometry.MultiPolygon(coords)
                    else:
                        continue
                    geometry_list.append(ee_geom)
                except Exception:
                    continue

            if not geometry_list:
                return [], "No valid Polygon or MultiPolygon geometries found in Shapefile."
            return geometry_list, None

    except Exception as exc:
        return [], f"Unexpected error while processing Shapefile: {exc}"


def parse_zip_shapefile_full(file_storage):
    # Read ZIP -> GeoDataFrame -> build both overlay GeoJSON and EE geometry.
    def _fail(message: str):
        return None, None, None, None, None, None, None, None, message

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "uploaded.zip")
            try:
                with open(zip_path, "wb") as fh:
                    fh.write(file_storage.read())
            except Exception:
                return _fail("Couldn't read uploaded Zip file.")

            try:
                with zipfile.ZipFile(zip_path, "r") as zf:
                    zf.extractall(tmpdir)
            except zipfile.BadZipFile:
                return _fail("Not a valid zip file.")

            shp_files = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.lower().endswith(".shp")]
            if not shp_files:
                return _fail("Zip file does not contain .SHP file.")

            try:
                gdf = gpd.read_file(shp_files[0])
            except Exception as exc:
                return _fail(
                    "Failed to read Shapefile with GeoPandas: Missing companion files (.SHX, .DBF, .CPG, .PRJ) "
                    f"{exc}"
                )

            if gdf.empty:
                return _fail("Shapefile contains no features.")

            if gdf.crs is None:
                minx, miny, maxx, maxy = gdf.total_bounds
                if abs(minx) > 180 or abs(maxx) > 180 or abs(miny) > 90 or abs(maxy) > 90:
                    return _fail("Shapefile appears projected but has no CRS. Reproject to EPSG:4326.")
            else:
                gdf = gdf.to_crs(epsg=4326)

            def _coerce_property_value(value):
                if value is None:
                    return None
                try:
                    if pd.isna(value):
                        return None
                except Exception:
                    pass
                if isinstance(value, (datetime, date, pd.Timestamp)):
                    return value.isoformat()
                if hasattr(value, "item"):
                    try:
                        return value.item()
                    except Exception:
                        pass
                return value

            columns = list(gdf.columns)
            geom_idx = columns.index("geometry") if "geometry" in columns else -1
            field_names = [col for col in columns if col != "geometry"]
            rows_light = []
            for row_idx, row in enumerate(gdf.itertuples(index=False, name=None)):
                props = {}
                for col_idx, col_name in enumerate(columns):
                    if col_idx == geom_idx:
                        continue
                    value = _coerce_property_value(row[col_idx])
                    props[col_name] = value
                rows_light.append({"__fid": row_idx, "properties": props})

            minx, miny, maxx, maxy = gdf.total_bounds
            bbox_geom = box(minx, miny, maxx, maxy)
            bbox_ee = ee.Geometry.Rectangle([minx, miny, maxx, maxy])

            # Raw geometries for accurate EE clip (avoid simplification artifacts).
            talhoes_list_raw = []
            for geom in gdf.geometry:
                if geom is None or geom.is_empty:
                    continue
                try:
                    if geom.geom_type == "Polygon":
                        coords = [list(geom.exterior.coords)]
                        talhoes_list_raw.append(ee.Geometry.Polygon(coords))
                    elif geom.geom_type == "MultiPolygon":
                        coords = [list(p.exterior.coords) for p in geom.geoms]
                        talhoes_list_raw.append(ee.Geometry.MultiPolygon(coords))
                except Exception:
                    continue

            if not talhoes_list_raw:
                return _fail("No valid Polygon or MultiPolygon geometries in Shapefile.")

            # Overlay geometry is simplified for UI performance only.
            overlay_gdf = gdf.copy()
            overlay_gdf["geometry"] = overlay_gdf.geometry.simplify(
                tolerance=CEM_OVERLAY_SIMPLIFY_TOLERANCE, preserve_topology=True
            )
            try:
                overlay_gdf = overlay_gdf.explode(ignore_index=True)
            except TypeError:
                overlay_gdf = overlay_gdf.explode()
                overlay_gdf = overlay_gdf.reset_index(drop=True)

            # Exploded overlay drives the "talhoes" count and per-feature stats.
            talhoes_list = []
            talhoes_count = 0
            for geom in overlay_gdf.geometry:
                if geom is None or geom.is_empty:
                    continue
                try:
                    if geom.geom_type == "Polygon":
                        coords = [list(geom.exterior.coords)]
                        talhoes_list.append(ee.Geometry.Polygon(coords))
                        talhoes_count += 1
                    elif geom.geom_type == "MultiPolygon":
                        for part in geom.geoms:
                            if part is None or part.is_empty:
                                continue
                            coords = [list(part.exterior.coords)]
                            talhoes_list.append(ee.Geometry.Polygon(coords))
                            talhoes_count += 1
                except Exception:
                    continue

            if not talhoes_list:
                return _fail("No valid Polygon or MultiPolygon geometries in Shapefile.")

            talhoes_ee = ee.FeatureCollection([ee.Feature(g) for g in talhoes_list_raw]).geometry()
            overlay_geojson = json.loads(overlay_gdf[["geometry"]].to_json())
            # inject feature ids for popups
            for i, feat in enumerate(overlay_geojson.get("features", [])):
                props = feat.get("properties") or {}
                props["id"] = i
                feat["properties"] = props
            centroid = [bbox_geom.centroid.x, bbox_geom.centroid.y]

            bbox_bounds = [miny, minx, maxy, maxx]
            return (
                talhoes_ee,
                bbox_ee,
                overlay_geojson,
                centroid,
                talhoes_count,
                bbox_bounds,
                rows_light,
                field_names,
                None,
            )

    except Exception as exc:
        return _fail(f"Unexpected error while processing Shapefile: {exc}")


def parse_zip_shapefile_full_path(zip_path: str):
    try:
        class FileObj:
            def __init__(self, stream):
                self.stream = stream

            def read(self):
                return self.stream.read()

        with open(zip_path, "rb") as fh:
            return parse_zip_shapefile_full(FileObj(fh))
    except Exception as exc:
        return None, None, None, None, None, None, None, None, f"Unexpected error while processing Shapefile: {exc}"


def parse_kml(file_storage):
    try:
        file_storage.seek(0)
        try:
            tree = ET.parse(file_storage)
        except ET.ParseError:
            return [], "Invalid KML file. XML could not be parsed."

        root = tree.getroot()
        ns = {"kml": "http://www.opengis.net/kml/2.2"}
        placemarks = root.findall(".//kml:Placemark", ns)
        if not placemarks:
            return [], "No Placemark elements found in KML file."

        geometry_list = []
        for placemark in placemarks:
            coords_elem = placemark.find(".//kml:coordinates", ns)
            if coords_elem is None or not coords_elem.text:
                continue

            try:
                coords_raw = coords_elem.text.strip().split()
                coords = []

                for c in coords_raw:
                    lon, lat, *_ = map(float, c.split(","))
                    coords.append([lon, lat])

                if len(coords) < 3:
                    continue
                if coords[0] != coords[-1]:
                    coords.append(coords[0])

                ee_geom = ee.Geometry.Polygon([coords])
                geometry_list.append(ee_geom)
            except Exception:
                continue

        if not geometry_list:
            return [], "No valid polygon geometries could be constructed from KML."
        return geometry_list, None

    except Exception as exc:
        return [], f"Error processing KML file: {exc}"


def parse_topojson(file_storage):
    try:
        bytes_data = file_storage.read()
        try:
            topojson_data = json.loads(bytes_data)
        except json.JSONDecodeError:
            return [], "Invalid TopoJSON file. JSON could not be decoded."

        if topojson_data.get("type") != "Topology":
            return [], "File is not a valid TopoJSON (missing or invalid 'Topology' type)."

        arcs = topojson_data.get("arcs")
        objects = topojson_data.get("objects")
        transform = topojson_data.get("transform", {})
        scale = transform.get("scale", [1, 1])
        translate = transform.get("translate", [0, 0])

        if not arcs or not isinstance(arcs, list):
            return [], "TopoJSON file contains no valid arcs."
        if not objects or not isinstance(objects, dict):
            return [], "TopoJSON file contains no objects."

        decoded_arcs = []
        for arc in arcs:
            x, y = 0, 0
            points = []
            for dx, dy in arc:
                x += dx
                y += dy
                lon = x * scale[0] + translate[0]
                lat = y * scale[1] + translate[1]
                points.append([lon, lat])
            decoded_arcs.append(points)

        geometries = []
        for obj in objects.values():
            if obj.get("type") == "GeometryCollection":
                geometries.extend(obj.get("geometries", []))
            else:
                geometries.append(obj)

        if not geometries:
            return [], "TopoJSON file contains no geometries."

        geometry_list = []
        for geom_data in geometries:
            geom_type = geom_data.get("type")
            arcs_refs = geom_data.get("arcs")
            if not arcs_refs:
                continue

            try:
                if geom_type == "Polygon":
                    rings = []
                    for ring_refs in arcs_refs:
                        ring = []
                        for arc_ref in ring_refs:
                            arc_idx = abs(arc_ref)
                            arc_points = (
                                decoded_arcs[arc_idx]
                                if arc_ref >= 0
                                else list(reversed(decoded_arcs[arc_idx]))
                            )
                            ring.extend(arc_points)
                        rings.append(ring)
                    geometry_list.append(ee.Geometry.Polygon(rings))

                elif geom_type == "MultiPolygon":
                    polygons = []
                    for polygon_refs in arcs_refs:
                        rings = []
                        for ring_refs in polygon_refs:
                            ring = []
                            for arc_ref in ring_refs:
                                arc_idx = abs(arc_ref)
                                arc_points = (
                                    decoded_arcs[arc_idx]
                                    if arc_ref >= 0
                                    else list(reversed(decoded_arcs[arc_idx]))
                                )
                                ring.extend(arc_points)
                            rings.append(ring)
                        polygons.append(rings)
                    geometry_list.append(ee.Geometry.MultiPolygon(polygons))

            except Exception:
                continue

        if not geometry_list:
            return [], "No valid Polygon or MultiPolygon geometries could be constructed from TopoJSON."
        return geometry_list, None

    except Exception as exc:
        return [], f"Error processing TopoJSON file: {exc}"


def parse_geojson(file_storage):
    try:
        bytes_data = file_storage.read()
        try:
            geojson_data = json.loads(bytes_data)
        except json.JSONDecodeError:
            return [], "Invalid GeoJSON file. JSON could not be decoded."

        # Prefer GeoPandas to handle CRS and reproject to EPSG:4326
        try:
            gdf = gpd.read_file(io.BytesIO(bytes_data))
            if gdf.empty:
                return [], "GeoJSON contains no features."
            if gdf.crs is None:
                # Heuristic: if coordinates are out of lon/lat range, require CRS
                minx, miny, maxx, maxy = gdf.total_bounds
                if abs(minx) > 180 or abs(maxx) > 180 or abs(miny) > 90 or abs(maxy) > 90:
                    return [], "GeoJSON appears to be projected but has no CRS. Please reproject to EPSG:4326."
            else:
                gdf = gdf.to_crs(epsg=4326)

            geometry_list = []
            for geom in gdf.geometry:
                if geom is None or geom.is_empty:
                    continue
                if geom.geom_type == "Polygon":
                    coords = [list(geom.exterior.coords)]
                    geometry_list.append(ee.Geometry.Polygon(coords))
                elif geom.geom_type == "MultiPolygon":
                    coords = [list(p.exterior.coords) for p in geom.geoms]
                    geometry_list.append(ee.Geometry.MultiPolygon(coords))

            if not geometry_list:
                return [], "No valid Polygon or MultiPolygon geometries found in GeoJSON."
            return geometry_list, None

        except Exception:
            # Fallback: raw GeoJSON parsing
            if "features" in geojson_data and isinstance(geojson_data["features"], list):
                features = geojson_data["features"]
            elif "geometries" in geojson_data and isinstance(geojson_data["geometries"], list):
                features = [{"geometry": geo} for geo in geojson_data["geometries"]]
            else:
                return [], "Unsupported GeoJSON structure. Must contain a 'features' or 'geometries' field."

            geometry_list = []
            for feature in features:
                if "geometry" in feature and "coordinates" in feature["geometry"]:
                    coordinates = feature["geometry"]["coordinates"]
                    geometry_type = feature["geometry"]["type"]
                    geometry = (
                        ee.Geometry.Polygon(coordinates)
                        if geometry_type == "Polygon"
                        else ee.Geometry.MultiPolygon(coordinates)
                    )
                    geometry_list.append(geometry)

            if not geometry_list:
                return [], "No valid geometries in the GeoJSON file. Ensure it contains Polygon or MultiPolygon features."
            return geometry_list, None

    except Exception as exc:
        return [], f"Error processing GeoJSON file: {exc} Please verify the file is valid."


def load_geojson_path(path: str):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh), None
    except Exception as exc:
        return None, f"Failed to read GeoJSON from {path}: {exc}"


def geojson_to_ee_geometry(geojson_data):
    if "features" in geojson_data and isinstance(geojson_data["features"], list):
        features = geojson_data["features"]
    elif "geometries" in geojson_data and isinstance(geojson_data["geometries"], list):
        features = [{"geometry": geo} for geo in geojson_data["geometries"]]
    else:
        return None

    geometry_list = []
    for feature in features:
        if "geometry" in feature and "coordinates" in feature["geometry"]:
            coordinates = feature["geometry"]["coordinates"]
            geometry_type = feature["geometry"]["type"]
            geometry = (
                ee.Geometry.Polygon(coordinates)
                if geometry_type == "Polygon"
                else ee.Geometry.MultiPolygon(coordinates)
            )
            geometry_list.append(geometry)

    if not geometry_list:
        return None
    return ee.Geometry.MultiPolygon(geometry_list)


# ---------------------------
# Upload handling
# ---------------------------

def upload_files_proc(upload_files: List) -> Tuple[ee.Geometry, List[str], List[float] | None]:
    geometry_aoi_list = []
    error_messages = []
    last_centroid = None

    for upload_file in upload_files:
        file_name = (upload_file.filename or "").lower()
        upload_file.stream.seek(0)

        if file_name.endswith(".gpkg"):
            with tempfile.NamedTemporaryFile(suffix=".gpkg", delete=False) as tmp:
                tmp.write(upload_file.read())
                tmp.flush()
                gpkg_geoms, error = parse_geopackage(tmp.name)
            if error:
                error_messages.append(f"→ {file_name}: {error}")
            geometry_aoi_list.extend(gpkg_geoms)
            if gpkg_geoms:
                last_centroid = gpkg_geoms[0].centroid(maxError=1).getInfo()["coordinates"]
            continue

        if file_name.endswith(".csv"):
            csv_geoms, error = parse_csv(upload_file)
            if error:
                error_messages.append(f"→ {file_name}: {error}")
            geometry_aoi_list.extend(csv_geoms)
            if csv_geoms:
                last_centroid = csv_geoms[0].centroid(maxError=1).getInfo()["coordinates"]
            continue

        if file_name.endswith(".kml"):
            kml_geoms, error = parse_kml(upload_file)
            if error:
                error_messages.append(f"→ {file_name}: {error}")
            geometry_aoi_list.extend(kml_geoms)
            if kml_geoms:
                last_centroid = kml_geoms[0].centroid(maxError=1).getInfo()["coordinates"]
            continue

        if file_name.endswith(".zip"):
            shp_geoms, error = parse_zip_shapefile(upload_file)
            if error:
                error_messages.append(f"→ {file_name}: {error}")
            geometry_aoi_list.extend(shp_geoms)
            if shp_geoms:
                last_centroid = shp_geoms[0].centroid(maxError=1).getInfo()["coordinates"]
            continue

        if file_name.endswith(".topojson"):
            topojson_geoms, error = parse_topojson(upload_file)
            if error:
                error_messages.append(f"→ {file_name}: {error}")
            geometry_aoi_list.extend(topojson_geoms)
            if topojson_geoms:
                last_centroid = topojson_geoms[0].centroid(maxError=1).getInfo()["coordinates"]
            continue

        if file_name.endswith(".geojson") or file_name.endswith(".json"):
            geojson_geoms, error = parse_geojson(upload_file)
            if error:
                error_messages.append(f"→ {file_name}: {error}")
            geometry_aoi_list.extend(geojson_geoms)
            if geojson_geoms:
                last_centroid = geojson_geoms[0].centroid(maxError=1).getInfo()["coordinates"]
            continue

    if geometry_aoi_list:
        geometry_aoi = ee.Geometry.MultiPolygon(geometry_aoi_list)
    else:
        geometry_aoi = ee.Geometry.Point([16.25, 36.65])

    return geometry_aoi, error_messages, last_centroid


# ---------------------------
# Date utils
# ---------------------------

def date_input_proc(input_date: date, time_range: int) -> Tuple[str, str]:
    end_date = input_date
    start_date = input_date - timedelta(days=time_range)
    return start_date.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d")


# ---------------------------
# Palette
# ---------------------------

def get_palettes(accessibility: str):
    default_ndvi_palette = ["#ffffe5", "#f7fcb9", "#78c679", "#41ab5d", "#238443", "#005a32"]
    default_reclassified = ["#a50026", "#ed5e3d", "#f9f7ae", "#f4ff78", "#9ed569", "#229b51", "#006837"]

    if accessibility == "Deuteranopia":
        return (
            ["#fffaa1", "#f4ef8e", "#9a5d67", "#573f73", "#372851", "#191135"],
            ["#95a600", "#92ed3e", "#affac5", "#78ffb0", "#69d6c6", "#22459c", "#000e69"],
        )
    if accessibility == "Protanopia":
        return (
            ["#a6f697", "#7def75", "#2dcebb", "#1597ab", "#0c677e", "#002c47"],
            ["#95a600", "#92ed3e", "#affac5", "#78ffb0", "#69d6c6", "#22459c", "#000e69"],
        )
    if accessibility == "Tritanopia":
        return (
            ["#cdffd7", "#a1fbb6", "#6cb5c6", "#3a77a5", "#205080", "#001752"],
            ["#ed4700", "#ed8a00", "#e1fabe", "#99ff94", "#87bede", "#2e40cf", "#0600bc"],
        )
    if accessibility == "Achromatopsia":
        return (
            ["#407de0", "#2763da", "#394388", "#272c66", "#16194f", "#010034"],
            ["#004f3d", "#338796", "#66a4f5", "#3683ff", "#3d50ca", "#421c7f", "#290058"],
        )

    return default_ndvi_palette, default_reclassified


def detect_macroregion(lat: float, lon: float) -> str:
    # Heuristic by latitude/longitude bands (Brazil macroregions)
    if lat > 2:
        return "Norte"
    if lat > -10 and lon > -50:
        return "Nordeste"
    if lat > -20 and lon < -50:
        return "Centro-Oeste"
    if lat > -26 and lon > -55:
        return "Sudeste"
    return "Sul"


def classify_ndvi_status(ndvi: float, crop: str) -> str:
    low, high = CROP_THRESHOLDS.get(crop, (0.45, 0.60))
    if ndvi < low:
        return "Baixo"
    if ndvi < high:
        return "Medio"
    return "Alto"


def geojson_geom_to_ee(geom: Dict[str, Any]):
    gtype = geom.get("type")
    coords = geom.get("coordinates")
    if not gtype or coords is None:
        return None
    if gtype == "Polygon":
        return ee.Geometry.Polygon(coords)
    if gtype == "MultiPolygon":
        return ee.Geometry.MultiPolygon(coords)
    return None


def geojson_to_ee_featurecollection(geojson_data: Dict[str, Any]):
    features = geojson_data.get("features")
    if not isinstance(features, list):
        return None
    ee_features = []
    for feat in features:
        geom = feat.get("geometry")
        if not geom:
            continue
        ee_geom = geojson_geom_to_ee(geom)
        if ee_geom is None:
            continue
        ee_features.append(ee.Feature(ee_geom))
    if not ee_features:
        return None
    return ee.FeatureCollection(ee_features)


def clip_ndvi_with_batches(ndvi_image, overlay_geojson: Dict[str, Any]):
    # For many features, split into two clips to avoid EE geometry limits.
    features = overlay_geojson.get("features", []) if overlay_geojson else []
    if not features:
        return ndvi_image
    if len(features) <= CEM_CLIP_BATCH_THRESHOLD:
        fc = geojson_to_ee_featurecollection({"type": "FeatureCollection", "features": features})
        if fc is None:
            return ndvi_image
        return ndvi_image.clip(fc.geometry())
    mid = len(features) // 2
    fc1 = geojson_to_ee_featurecollection({"type": "FeatureCollection", "features": features[:mid]})
    fc2 = geojson_to_ee_featurecollection({"type": "FeatureCollection", "features": features[mid:]})
    if fc1 is None or fc2 is None:
        return ndvi_image
    img1 = ndvi_image.clip(fc1.geometry())
    img2 = ndvi_image.clip(fc2.geometry())
    return ee.ImageCollection([img1, img2]).mosaic()


def compute_tiles_for_date(job: Dict[str, Any], date_str: str):
    # Build EE tiles for RGB, raw NDVI, and reclassified NDVI.
    if not job.get("bbox_bounds") or not job.get("overlay_geojson"):
        return None, "Job sem geometria."

    bbox = job["bbox_bounds"]
    geometry_aoi = ee.Geometry.Rectangle([bbox[1], bbox[0], bbox[3], bbox[2]])
    cloud = job.get("cloud", 100)
    accessibility = job.get("accessibility", "Normal")
    ndvi_palette, _ = get_palettes(accessibility)

    target_date = datetime.strptime(date_str, "%Y-%m-%d").date()
    start_date, end_date = date_input_proc(target_date, 7)
    collection, err = sat_collection(cloud, start_date, end_date, geometry_aoi)
    if err:
        return None, err

    # Use median composite for stability against clouds/outliers.
    imagery = collection.median()
    tci_params = {"bands": ["B4", "B3", "B2"], "min": 0, "max": 1, "gamma": 1}

    ndvi_image = imagery.normalizedDifference(["B8", "B4"])
    # Clamp extreme negatives for a dedicated "NoData/Negative" gray band.
    ndvi_image = ndvi_image.max(-0.05)

    # Clip NDVI to talhoes overlay (batching reduces geometry complexity).
    overlay_geojson = job.get("overlay_geojson_filtered") or job.get("overlay_geojson")
    ndvi_image = clip_ndvi_with_batches(ndvi_image, overlay_geojson)

    ndvi_palette = ["#6b7280"] + ndvi_palette
    ndvi_params = {"min": -0.05, "max": 1, "palette": ndvi_palette}

    # Bucket NDVI values into 7 classes + gray for negative/no-data.
    ndvi_classified = (
        ee.Image(ndvi_image)
        .where(ndvi_image.lt(0), 0)
        .where(ndvi_image.gte(0).And(ndvi_image.lt(0.15)), 1)
        .where(ndvi_image.gte(0.15).And(ndvi_image.lt(0.25)), 2)
        .where(ndvi_image.gte(0.25).And(ndvi_image.lt(0.35)), 3)
        .where(ndvi_image.gte(0.35).And(ndvi_image.lt(0.45)), 4)
        .where(ndvi_image.gte(0.45).And(ndvi_image.lt(0.65)), 5)
        .where(ndvi_image.gte(0.65).And(ndvi_image.lt(0.75)), 6)
        .where(ndvi_image.gte(0.75), 7)
    )
    _, reclassified_ndvi_palette = get_palettes(accessibility)
    reclassified_ndvi_palette = ["#6b7280"] + reclassified_ndvi_palette
    ndvi_class_params = {"min": 0, "max": 7, "palette": reclassified_ndvi_palette}

    rgb_map = ee.Image(imagery).getMapId(tci_params)
    ndvi_map = ee.Image(ndvi_image).getMapId(ndvi_params)
    ndvi_class_map = ee.Image(ndvi_classified).getMapId(ndvi_class_params)

    return {
        "date": date_str,
        "rgb": rgb_map["tile_fetcher"].url_format,
        "ndvi": ndvi_map["tile_fetcher"].url_format,
        "ndvi_class": ndvi_class_map["tile_fetcher"].url_format,
    }, None


# ---------------------------
# Routes
# ---------------------------

@app.route("/", methods=["GET"])
def index():
    return render_template("landing.html")


def _process_job(job_id: str, zip_path: str, form: Dict[str, str]):
    # Background worker: parse geometry, compute initial tiles, and populate JOBS.
    try:
        ee_initialize()

        accessibility = form.get("accessibility", "Normal")
        cloud_pixel_percentage = int(form.get("cloud_pixel_percentage", 100))
        initial_date = datetime.strptime(form.get("initial_date"), "%Y-%m-%d").date()
        end_date = form.get("end_date")
        selected_date = (
            datetime.strptime(end_date, "%Y-%m-%d").date() if end_date else initial_date
        )

        (
            talhoes_geom,
            bbox_ee,
            overlay_geojson,
            last_centroid,
            talhoes_count,
            bbox_bounds,
            rows_light,
            field_names,
            shp_err,
        ) = parse_zip_shapefile_full_path(zip_path)
        if shp_err:
            JOBS[job_id]["errors"].append(shp_err)
            return
        JOBS[job_id]["shapefile_ok"] = True
        import_id = JOBS[job_id].get("import_id")
        if import_id:
            examples_by_field = {name: [] for name in (field_names or [])}
            examples_seen = {name: set() for name in (field_names or [])}
            sample_limit = min(len(rows_light or []), 5000)
            for row in (rows_light or [])[:sample_limit]:
                props = row.get("properties") or {}
                for field_name in examples_by_field:
                    if len(examples_by_field[field_name]) >= 5:
                        continue
                    value = props.get(field_name)
                    if value is None:
                        continue
                    try:
                        if pd.isna(value):
                            continue
                    except Exception:
                        pass
                    value_str = str(value).strip()
                    if not value_str:
                        continue
                    if value_str in examples_seen[field_name]:
                        continue
                    examples_seen[field_name].add(value_str)
                    examples_by_field[field_name].append(value_str)
                if all(len(examples_by_field[name]) >= 5 for name in examples_by_field):
                    break
            IMPORTS[import_id] = {
                "rows": rows_light or [],
                "fields": field_names or [],
                "examples": examples_by_field,
            }

        region = "Sudeste"
        if last_centroid:
            region = detect_macroregion(last_centroid[1], last_centroid[0])
        crop = REGION_DEFAULT_CROP.get(region, "Soja/Milho")
        JOBS[job_id]["talhoes_count"] = talhoes_count
        JOBS[job_id]["bbox_bounds"] = bbox_bounds
        JOBS[job_id]["talhoes_ee"] = talhoes_geom
        JOBS[job_id]["overlay_geojson"] = overlay_geojson
        JOBS[job_id]["overlay_geojson_full"] = overlay_geojson
        JOBS[job_id]["overlay_geojson_filtered"] = None
        JOBS[job_id]["filter_active"] = False
        JOBS[job_id]["region"] = region
        JOBS[job_id]["crop"] = crop
        JOBS[job_id]["base_date"] = initial_date.strftime("%Y-%m-%d")
        JOBS[job_id]["selected_date"] = selected_date.strftime("%Y-%m-%d")

        # Precompute initial tiles so dashboard loads quickly.
        tiles, tiles_err = compute_tiles_for_date(JOBS[job_id], JOBS[job_id]["selected_date"])
        if tiles_err:
            JOBS[job_id]["errors"].append(tiles_err)
            return
        JOBS[job_id]["tiles_full"] = tiles
        JOBS[job_id]["tiles"] = tiles
        JOBS[job_id]["tiles_filter_active"] = False
    except Exception as exc:
        err = f"Erro inesperado: {exc}"
        JOBS[job_id]["errors"].append(err)
    finally:
        JOBS[job_id]["done"] = True
        try:
            os.remove(zip_path)
        except Exception:
            pass


def _monthly_date_list(end_date: date, days: int) -> List[str]:
    return [
        (end_date - timedelta(days=i)).strftime("%Y-%m-%d")
        for i in range(days - 1, -1, -1)
    ]


def _process_monthly(job_id: str, end_str: str, days: int) -> None:
    # Background worker: compute NDVI tiles for a rolling window.
    job = JOBS.get(job_id)
    if not job:
        return

    try:
        ee_initialize()

        end_date = datetime.strptime(end_str, "%Y-%m-%d").date()
        dates = _monthly_date_list(end_date, days)

        job["monthly"] = {
            "status": "running",
            "total": len(dates),
            "completed": 0,
            "items": [{"date": d, "ok": False, "error": None} for d in dates],
            "results": [],
            "error": None,
        }

        for idx, d_str in enumerate(dates):
            tiles, err = compute_tiles_for_date(job, d_str)
            item = job["monthly"]["items"][idx]
            if err:
                item["error"] = err
            else:
                item["ok"] = True
                job["monthly"]["results"].append(tiles)
            job["monthly"]["completed"] = idx + 1

        job["monthly"]["status"] = "done"
    except Exception as exc:
        job["monthly"]["status"] = "error"
        job["monthly"]["error"] = f"Erro na análise mensal: {exc}"


@app.route("/submit", methods=["POST"])
def submit():
    upload_files = request.files.getlist("aoi_files")
    if not upload_files:
        return render_template("landing.html", form_errors=["Envie um arquivo .zip de Shapefile para processar."])

    shp_file = None
    for f in upload_files:
        if (f.filename or "").lower().endswith(".zip"):
            shp_file = f
            break
    if shp_file is None:
        return render_template("landing.html", form_errors=["Formato inválido. Envie um .zip contendo o Shapefile."])

    tmpdir = tempfile.mkdtemp()
    zip_path = os.path.join(tmpdir, "upload.zip")
    shp_file.save(zip_path)

    # Opção A: validar ZIP antes de criar job
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            has_shp = any(name.lower().endswith(".shp") for name in zf.namelist())
    except zipfile.BadZipFile:
        try:
            os.remove(zip_path)
        except Exception:
            pass
        return render_template("landing.html", form_errors=["Arquivo ZIP inválido."])
    if not has_shp:
        try:
            os.remove(zip_path)
        except Exception:
            pass
        return render_template("landing.html", form_errors=["ZIP sem shapefile (.shp)."])

    job_id = uuid.uuid4().hex
    import_id = uuid.uuid4().hex
    JOBS[job_id] = {
        "errors": [],
        "done": False,
        "shapefile_ok": False,
        "talhoes_count": None,
        "bbox_bounds": None,
        "overlay_geojson": None,
        "overlay_geojson_full": None,
        "overlay_geojson_filtered": None,
        "tiles": None,
        "tiles_full": None,
        "tiles_filter_active": False,
        "region": None,
        "crop": None,
        "cloud": None,
        "accessibility": None,
        "base_date": None,
        "selected_date": None,
        "filter_active": False,
        "monthly": {
            "status": "idle",
            "total": 0,
            "completed": 0,
            "items": [],
            "results": [],
            "error": None,
        },
        "import_id": import_id,
    }

    initial_date = request.form.get("initial_date") or (date.today() - timedelta(days=2)).strftime("%Y-%m-%d")
    end_date = request.form.get("end_date") or ""

    form = {
        "accessibility": request.form.get("accessibility", "Normal"),
        "cloud_pixel_percentage": request.form.get("cloud_pixel_percentage", "100"),
        "initial_date": initial_date,
        "end_date": end_date,
    }

    JOBS[job_id]["cloud"] = int(form["cloud_pixel_percentage"])
    JOBS[job_id]["accessibility"] = form["accessibility"]

    thread = threading.Thread(target=_process_job, args=(job_id, zip_path, form), daemon=True)
    thread.start()

    return render_template("loading.html", job_id=job_id)


@app.route("/result/<job_id>")
def result(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return "<div class='error'>Job não encontrado.</div>", 404

    if job["errors"]:
        errors = "".join(f"<div class='error'>{e}</div>" for e in job["errors"])
        return f"<div class='errors'>{errors}</div>"

    if job["done"]:
        return "OK"

    return "<div class='map-placeholder'>Processamento em andamento...</div>"


@app.route("/view/<job_id>")
def view(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return render_template("landing.html", form_errors=["Job não encontrado."])

    if job["errors"]:
        return render_template("landing.html", form_errors=job["errors"])

    return render_template(
        "smartgreen-dashboard.html",
        job_id=job_id,
        talhoes_count=job.get("talhoes_count"),
        region=job.get("region"),
        crop=job.get("crop"),
        selected_date=job.get("selected_date"),
        base_date=job.get("base_date"),
    )


@app.route("/api/job/<job_id>/meta")
def job_meta(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    return {
        "done": job.get("done"),
        "errors": job.get("errors") or [],
        "shapefile_ok": job.get("shapefile_ok"),
        "talhoes_count": job.get("talhoes_count"),
        "region": job.get("region"),
        "crop": job.get("crop"),
        "crop_options": CROP_OPTIONS,
        "bbox": job.get("bbox_bounds"),
        "date": job.get("selected_date"),
        "base_date": job.get("base_date"),
        "import_id": job.get("import_id"),
    }


@app.route("/api/import/<import_id>/fields")
def import_fields(import_id: str):
    imp = IMPORTS.get(import_id)
    if not imp:
        return {"error": "Import not found"}, 404

    field_names = imp.get("fields") or []
    if not field_names:
        return {"fields": []}
    examples_by_field = imp.get("examples") or {}

    fields_payload = []
    for field_name in field_names:
        fields_payload.append({"name": field_name, "examples": examples_by_field.get(field_name, [])})

    return {"fields": fields_payload}


@app.route("/api/import/<import_id>/filter-schema", methods=["POST"])
def import_filter_schema(import_id: str):
    imp = IMPORTS.get(import_id)
    if not imp:
        return {"error": "Import not found"}, 404

    payload = request.get_json(silent=True) or {}
    farm_field = payload.get("farm_field")
    field_field = payload.get("field_field")
    id_field = payload.get("id_field")

    if not farm_field or not field_field:
        return {"error": "farm_field and field_field are required"}, 400
    if farm_field == field_field:
        return {"error": "farm_field must be different from field_field"}, 400
    if id_field:
        if id_field == farm_field or id_field == field_field:
            return {"error": "id_field must be different from farm_field and field_field"}, 400

    imp["filter_schema"] = {
        "farm_field": farm_field,
        "field_field": field_field,
        "id_field": id_field,
    }
    return {"ok": True}


def _normalize_value(value: Any) -> Optional[str]:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    value_str = str(value).strip()
    return value_str or None


def _build_field_label(talhao_value: Optional[str], id_value: Optional[str]) -> Optional[str]:
    if talhao_value and id_value:
        return f"Talhão {talhao_value} (ID {id_value})"
    if talhao_value:
        return f"Talhão {talhao_value}"
    if id_value:
        return f"ID {id_value}"
    return None


@app.route("/api/import/<import_id>/farms")
def import_farms(import_id: str):
    imp = IMPORTS.get(import_id)
    if not imp:
        return {"error": "Import not found"}, 404

    schema = imp.get("filter_schema")
    if not schema:
        return {"error": "filter_schema not set"}, 400

    farm_field = schema.get("farm_field")
    field_field = schema.get("field_field")
    id_field = schema.get("id_field")
    rows = imp.get("rows") or []
    q = request.args.get("q", "").strip().lower()

    farms_map: Dict[str, List[str]] = {}
    for row in rows:
        props = row.get("properties") or {}
        farm_value = _normalize_value(props.get(farm_field))
        if not farm_value:
            continue
        if q and q not in farm_value.lower():
            continue
        if farm_value not in farms_map:
            farms_map[farm_value] = []
        if len(farms_map[farm_value]) >= 5:
            continue
        talhao_value = _normalize_value(props.get(field_field))
        id_value = _normalize_value(props.get(id_field)) if id_field else None
        label = _build_field_label(talhao_value, id_value)
        if label and label not in farms_map[farm_value]:
            farms_map[farm_value].append(label)

    farms_payload = [
        {"value": farm_value, "examples": farms_map[farm_value]}
        for farm_value in sorted(farms_map.keys(), key=lambda v: v.lower())
    ]
    return {"farms": farms_payload}


@app.route("/api/import/<import_id>/fields-by-farm")
def import_fields_by_farm(import_id: str):
    imp = IMPORTS.get(import_id)
    if not imp:
        return {"error": "Import not found"}, 404

    schema = imp.get("filter_schema")
    if not schema:
        return {"error": "filter_schema not set"}, 400

    farm = request.args.get("farm")
    if not farm:
        return {"error": "farm is required"}, 400

    farm_norm = _normalize_value(farm)
    if not farm_norm:
        return {"fields": []}

    farm_field = schema.get("farm_field")
    field_field = schema.get("field_field")
    id_field = schema.get("id_field")
    rows = imp.get("rows") or []
    q = request.args.get("q", "").strip().lower()

    results = []
    for row in rows:
        props = row.get("properties") or {}
        row_farm = _normalize_value(props.get(farm_field))
        if row_farm != farm_norm:
            continue
        talhao_value = _normalize_value(props.get(field_field))
        id_value = _normalize_value(props.get(id_field)) if id_field else None
        if q:
            talhao_match = talhao_value and q in talhao_value.lower()
            id_match = id_value and q in id_value.lower()
            if not (talhao_match or id_match):
                continue
        label = _build_field_label(talhao_value, id_value) or "Talhão"
        results.append(
            {
                "feature_id": row.get("__fid"),
                "label": label,
                "talhao_value": talhao_value,
                "id_value": id_value,
            }
        )
        if len(results) >= 500:
            break

    return {"fields": results}


def _filter_overlay_by_feature_ids(overlay_geojson: Dict[str, Any], feature_ids: List[int]):
    if not overlay_geojson or not feature_ids:
        return None
    id_set = {int(fid) for fid in feature_ids if fid is not None}
    if not id_set:
        return None
    features = overlay_geojson.get("features", [])
    if not isinstance(features, list):
        return None
    filtered = [feat for feat in features if int((feat.get("properties") or {}).get("id", -1)) in id_set]
    if not filtered:
        return {"type": "FeatureCollection", "features": []}
    return {"type": "FeatureCollection", "features": filtered}


@app.route("/api/job/<job_id>/filter-selection", methods=["POST"])
def job_filter_selection(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    if not job.get("overlay_geojson_full"):
        return {"error": "Job sem geometria."}, 400

    payload = request.get_json(silent=True) or {}
    feature_ids = payload.get("feature_ids") or []
    if not isinstance(feature_ids, list):
        return {"error": "feature_ids must be a list"}, 400

    if not feature_ids:
        job["overlay_geojson_filtered"] = None
        job["filter_active"] = False
        if job.get("tiles_full"):
            job["tiles"] = job["tiles_full"]
            job["tiles_filter_active"] = False
        return {"ok": True, "filtered": False}

    filtered_overlay = _filter_overlay_by_feature_ids(job.get("overlay_geojson_full"), feature_ids)
    if filtered_overlay is None:
        return {"error": "Invalid feature_ids"}, 400

    job["overlay_geojson_filtered"] = filtered_overlay
    job["filter_active"] = True
    date_str = job.get("selected_date")
    if not date_str:
        return {"ok": True, "filtered": True}
    tiles, err = compute_tiles_for_date(job, date_str)
    if err:
        return {"error": err}, 400
    job["tiles"] = tiles
    job["tiles_filter_active"] = True
    return {"ok": True, "filtered": True, "tiles": tiles}


@app.route("/api/job/<job_id>/talhoes")
def job_talhoes(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    if not job.get("overlay_geojson"):
        if job.get("errors"):
            return {"error": job["errors"][0]}, 400
        return {"error": "Talhoes not found"}, 404
    return job["overlay_geojson"]


@app.route("/api/job/<job_id>/tiles")
def job_tiles(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    date_str = request.args.get("date")
    cached_tiles = job.get("tiles")
    if date_str and cached_tiles:
        cached_date = cached_tiles.get("date") or job.get("selected_date")
        cached_filter = job.get("tiles_filter_active")
        if cached_date == date_str and cached_filter == job.get("filter_active"):
            return cached_tiles
    if not date_str:
        if job.get("tiles"):
            return job["tiles"]
        return {"error": "Date required"}, 400

    tiles, err = compute_tiles_for_date(job, date_str)
    if err:
        return {"error": err}, 400
    job["tiles"] = tiles
    job["tiles_filter_active"] = job.get("filter_active")
    return tiles


@app.route("/api/job/<job_id>/monthly/start", methods=["POST"])
def monthly_start(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    if not job.get("done"):
        return {"error": "Job not ready"}, 400
    if job.get("shapefile_ok") is not True:
        return {"error": "Shapefile não encontrado no ZIP."}, 400

    monthly = job.get("monthly") or {}
    if monthly.get("status") == "running":
        return {"error": "Monthly analysis already running"}, 409

    end_str = request.args.get("end") or request.form.get("end")
    if not end_str:
        end_str = job.get("selected_date") or job.get("base_date")
    if not end_str:
        return {"error": "End date required"}, 400

    try:
        days = int(request.args.get("days") or request.form.get("days") or "31")
    except ValueError:
        return {"error": "Invalid days"}, 400
    if days < 1:
        return {"error": "Days must be >= 1"}, 400

    thread = threading.Thread(
        target=_process_monthly, args=(job_id, end_str, days), daemon=True
    )
    thread.start()
    return {"ok": True, "status": "running", "end": end_str, "days": days}


@app.route("/api/job/<job_id>/monthly/status")
def monthly_status(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    monthly = job.get("monthly") or {}
    total = int(monthly.get("total") or 0)
    completed = int(monthly.get("completed") or 0)
    percent = int((completed / total) * 100) if total > 0 else 0
    payload = {
        "status": monthly.get("status", "idle"),
        "total": total,
        "completed": completed,
        "percent": percent,
        "items": monthly.get("items") or [],
    }
    if monthly.get("status") == "done":
        payload["results"] = monthly.get("results") or []
    if monthly.get("status") == "error":
        payload["error"] = monthly.get("error") or "Erro na análise mensal."
    return payload


@app.route("/api/job/<job_id>/series")
def job_series(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    days = int(request.args.get("days", "30"))
    end_str = request.args.get("end")
    if not end_str:
        end_str = job.get("base_date") or job.get("selected_date")
    if not end_str:
        return {"error": "End date required"}, 400

    end_date = datetime.strptime(end_str, "%Y-%m-%d").date()
    series = []
    for i in range(days):
        d = end_date - timedelta(days=i)
        d_str = d.strftime("%Y-%m-%d")
        tiles, err = compute_tiles_for_date(job, d_str)
        if err:
            continue
        series.append(tiles)

    return {"series": list(reversed(series))}


@app.route("/api/job/<job_id>/talhao_stats")
def talhao_stats(job_id: str):
    job = JOBS.get(job_id)
    if not job or not job.get("overlay_geojson"):
        return {"error": "Job not found"}, 404

    idx = request.args.get("index")
    date_str = request.args.get("date")
    if idx is None or date_str is None:
        return {"error": "index and date required"}, 400

    try:
        idx = int(idx)
    except ValueError:
        return {"error": "invalid index"}, 400

    features = job["overlay_geojson"].get("features", [])
    if idx < 0 or idx >= len(features):
        return {"error": "index out of range"}, 400

    geom = geojson_geom_to_ee(features[idx]["geometry"])
    if geom is None:
        return {"error": "invalid geometry"}, 400

    # recompute NDVI image for stats
    bbox = job["bbox_bounds"]
    geometry_aoi = ee.Geometry.Rectangle([bbox[1], bbox[0], bbox[3], bbox[2]])
    cloud = job.get("cloud", 100)
    start_date, end_date = date_input_proc(datetime.strptime(date_str, "%Y-%m-%d").date(), 7)
    collection, err = sat_collection(cloud, start_date, end_date, geometry_aoi)
    if err:
        return {"error": err}, 400
    imagery = collection.median()
    ndvi_image = imagery.normalizedDifference(["B8", "B4"])
    ndvi_image = ndvi_image.max(0)

    # Pixel validity metrics
    # Considera válidos apenas pixels com máscara ativa e NDVI > 0 (evita nuvem/baixa iluminação que vira cinza).
    valid_mask = ndvi_image.mask().And(ndvi_image.gt(0))

    # Conta de pixels válidos (usa count em imagem constante mascarada pelos válidos).
    valid_pixels = ee.Image.constant(1).updateMask(valid_mask).reduceRegion(
        reducer=ee.Reducer.count(),
        geometry=geom,
        scale=10,
        maxPixels=1e13,
    )

    # Total de pixels no talhão (sem máscara).
    total_pixels = ee.Image.constant(1).reduceRegion(
        reducer=ee.Reducer.count(),
        geometry=geom,
        scale=10,
        maxPixels=1e13,
    )

    valid_pixels_value = ee.Number(valid_pixels.get("constant", 0))
    total_pixels_value = ee.Number(total_pixels.get("constant", 0))

    invalid_pixels_value = total_pixels_value.subtract(valid_pixels_value).max(0)
    percent_valid = ee.Algorithms.If(
        total_pixels_value.eq(0),
        ee.Number(0),
        valid_pixels_value.divide(total_pixels_value).multiply(100),
    )
    percent_valid = ee.Number(percent_valid).min(100)

    pixels_validos = int(valid_pixels_value.getInfo() or 0)
    pixels_totais = int(total_pixels_value.getInfo() or 0)
    pixels_invalidos = int(invalid_pixels_value.getInfo() or 0)
    percentual_validos = float(percent_valid.getInfo() or 0.0)
    mensagem = "Sem dados disponíveis" if pixels_totais == 0 else None

    stats = ndvi_image.updateMask(valid_mask).reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=geom,
        scale=10,
        maxPixels=1e9,
    ).getInfo()
    ndvi_mean = stats.get("nd", 0)

    area_ha = geom.area(maxError=1).getInfo() / 10000.0
    crop = job.get("crop", "Soja/Milho")
    status = classify_ndvi_status(ndvi_mean, crop)

    return {
        "ndvi": ndvi_mean,  # compatibilidade: chave existente
        "ndvi_medio": ndvi_mean,  # novo alias solicitado
        "area_ha": area_ha,
        "crop": crop,
        "status": status,
        "pixels_validos": pixels_validos,
        "pixels_totais": pixels_totais,
        "pixels_invalidos": pixels_invalidos,
        "percentual_validos": percentual_validos,
        "mensagem": mensagem,
    }


@app.route("/api/job/<job_id>/set_crop", methods=["POST"])
def set_crop(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        return {"error": "Job not found"}, 404
    crop = request.args.get("crop") or request.form.get("crop")
    if crop not in CROP_OPTIONS:
        return {"error": "Invalid crop"}, 400
    job["crop"] = crop
    return {"ok": True, "crop": crop}


if __name__ == "__main__":
    ee_initialize()
    app.run(host="0.0.0.0", port=int(os.getenv("PORT", "8000")), debug=True)
