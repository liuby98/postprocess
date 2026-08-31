#!/usr/bin/env python3
"""Convert MOD10CM Collection 6.1 HDF-EOS2 granules to one CF-NetCDF.

The script reads monthly HDF4 granules one at a time, extracts
``Snow_Cover_Monthly_CMG`` and ``Snow_Spatial_QA``, masks invalid or non-good
quality cells, subsets a regular latitude/longitude region, sorts granules by
their acquisition month, and writes one time-ordered NetCDF file. No
intermediate per-month NetCDF files are created.
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

try:
    from netCDF4 import Dataset, date2num, num2date
except ImportError:  # pragma: no cover
    Dataset = date2num = num2date = None

try:
    from pyhdf.error import HDF4Error
    from pyhdf.SD import SD, SDC
except ImportError:  # pragma: no cover
    SD = SDC = None
    HDF4Error = Exception


warnings.filterwarnings(
    "ignore",
    message="Setting the shape on a NumPy array has been deprecated.*",
    category=DeprecationWarning,
)


DEFAULT_ROOT = Path(
    "/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/"
    "Snow_Cover_Fraction/MODIS_MOD10CM"
)
DEFAULT_INPUT_DIR = DEFAULT_ROOT / "raw_hdf"
DEFAULT_OUTPUT = DEFAULT_ROOT / "MOD10CM_SCF_2001_2017.nc"
DEFAULT_PATTERN = "MOD10CM.A*.061.*.hdf"
DEFAULT_BBOX = (65.0, 5.0, 145.0, 60.0)

SNOW_NAME = "Snow_Cover_Monthly_CMG"
QA_NAME = "Snow_Spatial_QA"
HDF4_MAGIC = b"\x0e\x03\x13\x01"
OUTPUT_FILL = 255
EXPECTED_SHAPE = (3600, 7200)
FILENAME_RE = re.compile(
    r"^MOD10CM\.A(?P<year>\d{4})(?P<doy>\d{3})\."
    r"(?P<collection>\d{3})\.(?P<production>\d{13})\.hdf$"
)


@dataclass(frozen=True)
class Granule:
    path: Path
    month: datetime
    collection: str
    production: str


@dataclass(frozen=True)
class GridSubset:
    row_slice: slice
    col_slice: slice
    latitude: object
    longitude: object


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert MODIS Terra MOD10CM C6.1 monthly HDF4 granules to one "
            "quality-filtered, time-sorted CF-NetCDF file."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--pattern", default=DEFAULT_PATTERN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--end-year", type=int, default=2017)
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        default=DEFAULT_BBOX,
        metavar=("WEST", "SOUTH", "EAST", "NORTH"),
        help="Geographic subset; default: 65 5 145 60.",
    )
    parser.add_argument(
        "--ignore-qa",
        action="store_true",
        help="Keep valid 0-100 snow values regardless of Snow_Spatial_QA.",
    )
    parser.add_argument(
        "--strict-coverage",
        action="store_true",
        help="Fail if any requested calendar month is missing.",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=4,
        choices=range(1, 10),
        metavar="1-9",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Inspect every input and report coverage without writing.",
    )
    return parser


def require_dependencies() -> None:
    missing = []
    if np is None:
        missing.append("numpy")
    if Dataset is None:
        missing.append("netCDF4")
    if SD is None:
        missing.append("pyhdf")
    if missing:
        names = " ".join(missing)
        raise RuntimeError(
            f"Missing Python packages: {', '.join(missing)}. Install with:\n"
            f"  conda install -c conda-forge {names}"
        )


def parse_granule(path: Path) -> Granule:
    match = FILENAME_RE.match(path.name)
    if match is None:
        raise ValueError(
            f"Unexpected filename {path.name}; expected "
            "MOD10CM.AYYYYDDD.061.YYYYDDDhhmmss.hdf"
        )
    year = int(match.group("year"))
    doy = int(match.group("doy"))
    try:
        acquisition = datetime.strptime(f"{year:04d}{doy:03d}", "%Y%j").replace(
            tzinfo=timezone.utc
        )
    except ValueError as exc:
        raise ValueError(f"Invalid acquisition date in {path.name}: {exc}") from exc
    if acquisition.day != 1:
        raise ValueError(
            f"{path.name} resolves to {acquisition:%Y-%m-%d}, not the first "
            "day of a month"
        )
    return Granule(
        path=path,
        month=acquisition,
        collection=match.group("collection"),
        production=match.group("production"),
    )


def discover_granules(args: argparse.Namespace) -> list[Granule]:
    if args.start_year > args.end_year:
        raise ValueError("--start-year must not exceed --end-year")
    if not args.input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {args.input_dir}")
    paths = list(args.input_dir.glob(args.pattern))
    if not paths:
        raise FileNotFoundError(
            f"No files match {args.input_dir / args.pattern}. Confirm that "
            "downloads have normal .hdf filenames."
        )
    candidates = [parse_granule(path) for path in paths]
    candidates = [
        item
        for item in candidates
        if args.start_year <= item.month.year <= args.end_year
    ]
    if not candidates:
        raise ValueError(f"No granules fall within {args.start_year}-{args.end_year}")

    by_month: dict[tuple[int, int], Granule] = {}
    for item in candidates:
        key = (item.month.year, item.month.month)
        previous = by_month.get(key)
        current_rank = (item.collection, item.production, item.path.name)
        previous_rank = (
            (previous.collection, previous.production, previous.path.name)
            if previous is not None
            else None
        )
        if previous is None or current_rank > previous_rank:
            if previous is not None:
                print(
                    f"[WARN] Duplicate {item.month:%Y-%m}; keeping "
                    f"{item.path.name} instead of {previous.path.name}"
                )
            by_month[key] = item
        else:
            print(f"[WARN] Duplicate {item.month:%Y-%m}; ignoring {item.path.name}")
    return sorted(by_month.values(), key=lambda item: item.month)


def expected_months(start_year: int, end_year: int) -> list[tuple[int, int]]:
    return [
        (year, month)
        for year in range(start_year, end_year + 1)
        for month in range(1, 13)
    ]


def report_coverage(granules: list[Granule], args: argparse.Namespace) -> list[str]:
    available = {(item.month.year, item.month.month) for item in granules}
    missing = [
        f"{year:04d}-{month:02d}"
        for year, month in expected_months(args.start_year, args.end_year)
        if (year, month) not in available
    ]
    print(
        f"[INFO] Selected {len(granules)} unique monthly granules: "
        f"{granules[0].month:%Y-%m} to {granules[-1].month:%Y-%m}"
    )
    if missing:
        message = f"Missing {len(missing)} month(s): {', '.join(missing)}"
        if args.strict_coverage:
            raise ValueError(message)
        print(f"[WARN] {message}; missing months will not be filled with zero")
    else:
        print("[OK] Complete monthly coverage for the requested years")
    return missing


def verify_hdf4_signature(path: Path) -> None:
    with path.open("rb") as handle:
        signature = handle.read(4)
    if signature != HDF4_MAGIC:
        raise ValueError(
            f"{path} is not HDF4; it may be an HTML error or incomplete download"
        )


def inspect_hdf(path: Path) -> tuple[int, int]:
    verify_hdf4_signature(path)
    hdf = SD(str(path), SDC.READ)
    try:
        absent = sorted({SNOW_NAME, QA_NAME} - set(hdf.datasets()))
        if absent:
            raise KeyError(f"{path} lacks HDF datasets: {', '.join(absent)}")
        snow = hdf.select(SNOW_NAME)
        qa = hdf.select(QA_NAME)
        try:
            snow_shape = tuple(snow.info()[2])
            qa_shape = tuple(qa.info()[2])
        finally:
            snow.endaccess()
            qa.endaccess()
    finally:
        hdf.end()
    if snow_shape != qa_shape:
        raise ValueError(f"Snow/QA shapes differ in {path}: {snow_shape} vs {qa_shape}")
    if snow_shape != EXPECTED_SHAPE:
        raise ValueError(
            f"Unexpected grid in {path}: {snow_shape}; expected {EXPECTED_SHAPE}"
        )
    return snow_shape


def build_subset(shape: tuple[int, int], bbox: list[float]) -> GridSubset:
    west, south, east, north = map(float, bbox)
    if not (-180.0 <= west < east <= 180.0):
        raise ValueError("--bbox requires -180 <= WEST < EAST <= 180")
    if not (-90.0 <= south < north <= 90.0):
        raise ValueError("--bbox requires -90 <= SOUTH < NORTH <= 90")
    nlat, nlon = shape
    dlat = 180.0 / nlat
    dlon = 360.0 / nlon
    if not np.isclose(dlat, 0.05) or not np.isclose(dlon, 0.05):
        raise ValueError(f"Grid resolution is {dlon:g} x {dlat:g}, not 0.05")
    latitude_desc = 90.0 - (np.arange(nlat, dtype=np.float64) + 0.5) * dlat
    longitude = -180.0 + (np.arange(nlon, dtype=np.float64) + 0.5) * dlon
    rows = np.where((latitude_desc >= south) & (latitude_desc <= north))[0]
    cols = np.where((longitude >= west) & (longitude <= east))[0]
    if rows.size == 0 or cols.size == 0:
        raise ValueError(f"No cells fall inside bbox {bbox}")
    row_slice = slice(int(rows.min()), int(rows.max()) + 1)
    col_slice = slice(int(cols.min()), int(cols.max()) + 1)
    return GridSubset(
        row_slice=row_slice,
        col_slice=col_slice,
        latitude=latitude_desc[row_slice][::-1].astype(np.float32),
        longitude=longitude[col_slice].astype(np.float32),
    )


def read_hdf_subset(path: Path, subset: GridSubset) -> tuple[object, object]:
    start = (subset.row_slice.start, subset.col_slice.start)
    count = (
        subset.row_slice.stop - subset.row_slice.start,
        subset.col_slice.stop - subset.col_slice.start,
    )
    hdf = SD(str(path), SDC.READ)
    try:
        snow_sds = hdf.select(SNOW_NAME)
        qa_sds = hdf.select(QA_NAME)
        try:
            snow = np.asarray(snow_sds.get(start=start, count=count), dtype=np.uint8)
            qa = np.asarray(qa_sds.get(start=start, count=count), dtype=np.uint8)
        finally:
            snow_sds.endaccess()
            qa_sds.endaccess()
    finally:
        hdf.end()
    # Source rows run north-to-south; output latitude is ascending.
    return snow[::-1, :], qa[::-1, :]


def next_month(value: datetime) -> datetime:
    if value.month == 12:
        return datetime(value.year + 1, 1, 1, tzinfo=timezone.utc)
    return datetime(value.year, value.month + 1, 1, tzinfo=timezone.utc)


def create_output(
    path: Path,
    granules: list[Granule],
    subset: GridSubset,
    missing: list[str],
    args: argparse.Namespace,
) -> None:
    chunk_y = min(256, len(subset.latitude))
    chunk_x = min(256, len(subset.longitude))
    with Dataset(path, "w", format="NETCDF4") as nc_obj:
        nc_obj.createDimension("time", len(granules))
        nc_obj.createDimension("lat", len(subset.latitude))
        nc_obj.createDimension("lon", len(subset.longitude))
        nc_obj.createDimension("bnds", 2)

        time_var = nc_obj.createVariable("time", "f8", ("time",))
        bounds_var = nc_obj.createVariable("time_bnds", "f8", ("time", "bnds"))
        lat_var = nc_obj.createVariable("lat", "f4", ("lat",))
        lon_var = nc_obj.createVariable("lon", "f4", ("lon",))
        snow_var = nc_obj.createVariable(
            SNOW_NAME,
            "u1",
            ("time", "lat", "lon"),
            fill_value=np.uint8(OUTPUT_FILL),
            zlib=True,
            complevel=args.compression_level,
            shuffle=True,
            chunksizes=(1, chunk_y, chunk_x),
        )
        qa_var = nc_obj.createVariable(
            QA_NAME,
            "u1",
            ("time", "lat", "lon"),
            fill_value=np.uint8(OUTPUT_FILL),
            zlib=True,
            complevel=args.compression_level,
            shuffle=True,
            chunksizes=(1, chunk_y, chunk_x),
        )

        time_var.standard_name = "time"
        time_var.long_name = "middle of monthly averaging interval"
        time_var.units = "days since 1970-01-01 00:00:00"
        time_var.calendar = "proleptic_gregorian"
        time_var.bounds = "time_bnds"
        bounds_var.long_name = "monthly time interval bounds"
        lat_var.standard_name = "latitude"
        lat_var.long_name = "latitude of grid-cell center"
        lat_var.units = "degrees_north"
        lat_var.axis = "Y"
        lon_var.standard_name = "longitude"
        lon_var.long_name = "longitude of grid-cell center"
        lon_var.units = "degrees_east"
        lon_var.axis = "X"
        lat_var[:] = subset.latitude
        lon_var[:] = subset.longitude

        snow_var.long_name = "MODIS monthly mean snow-covered area"
        snow_var.units = "percent"
        snow_var.valid_min = np.uint8(0)
        snow_var.valid_max = np.uint8(100)
        snow_var.coordinates = "time lat lon"
        snow_var.grid_mapping = "latitude_longitude"
        snow_var.cell_methods = "time: mean area: mean"
        snow_var.comment = "Values 211/250/253/254/255 are missing. " + (
            "Snow_Spatial_QA was not used for filtering."
            if args.ignore_qa
            else "Only Snow_Spatial_QA=0 is retained."
        )
        qa_var.long_name = "MOD10CM monthly snow spatial quality flag"
        qa_var.units = "1"
        qa_var.flag_values = np.asarray([0, 1, 252, 254], dtype=np.uint8)
        qa_var.flag_meanings = "good_quality other_quality antarctica_mask water_mask"
        qa_var.coordinates = "time lat lon"

        crs = nc_obj.createVariable("latitude_longitude", "i4")
        crs.grid_mapping_name = "latitude_longitude"
        crs.longitude_of_prime_meridian = 0.0
        crs.semi_major_axis = 6378137.0
        crs.inverse_flattening = 298.257223563
        crs.epsg_code = "EPSG:4326"

        starts = [item.month for item in granules]
        ends = [next_month(item.month) for item in granules]
        start_num = date2num(starts, units=time_var.units, calendar=time_var.calendar)
        end_num = date2num(ends, units=time_var.units, calendar=time_var.calendar)
        bounds_var[:, 0] = start_num
        bounds_var[:, 1] = end_num
        time_var[:] = (np.asarray(start_num) + np.asarray(end_num)) / 2.0

        nc_obj.Conventions = "CF-1.8, ACDD-1.3"
        nc_obj.title = "MODIS Terra MOD10CM C6.1 monthly snow-cover fraction"
        nc_obj.summary = "Quality-filtered subset merged chronologically from HDF-EOS2."
        nc_obj.source = "MODIS/Terra MOD10CM Collection 6.1 HDF-EOS2 granules"
        nc_obj.product_doi = "10.5067/MODIS/MOD10CM.061"
        nc_obj.references = "https://doi.org/10.5067/MODIS/MOD10CM.061"
        nc_obj.history = (
            f"{datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}: created by "
            "preprocess_mod10cm.py"
        )
        nc_obj.processing_level = "quality-filtered geographic subset"
        nc_obj.source_granule_count = len(granules)
        nc_obj.missing_months = ",".join(missing) if missing else "none"
        nc_obj.geospatial_lat_min = float(subset.latitude.min())
        nc_obj.geospatial_lat_max = float(subset.latitude.max())
        nc_obj.geospatial_lon_min = float(subset.longitude.min())
        nc_obj.geospatial_lon_max = float(subset.longitude.max())
        nc_obj.geospatial_lat_resolution = 0.05
        nc_obj.geospatial_lon_resolution = 0.05
        nc_obj.time_coverage_start = f"{starts[0]:%Y-%m-%dT00:00:00Z}"
        nc_obj.time_coverage_end = f"{ends[-1]:%Y-%m-%dT00:00:00Z}"

        for index, item in enumerate(granules):
            snow, qa = read_hdf_subset(item.path, subset)
            valid = snow <= 100
            if not args.ignore_qa:
                valid &= qa == 0
            cleaned = np.full(snow.shape, OUTPUT_FILL, dtype=np.uint8)
            cleaned[valid] = snow[valid]
            snow_var[index, :, :] = cleaned
            qa_var[index, :, :] = qa
            if index == 0 or (index + 1) % 12 == 0 or index + 1 == len(granules):
                valid_fraction = 100.0 * np.count_nonzero(valid) / valid.size
                print(
                    f"[WRITE] {index + 1:3d}/{len(granules)} "
                    f"{item.month:%Y-%m}: valid={valid_fraction:.1f}%"
                )


def validate_output(path: Path, granules: list[Granule], subset: GridSubset) -> None:
    with Dataset(path) as nc_obj:
        expected = (
            len(granules),
            len(subset.latitude),
            len(subset.longitude),
        )
        actual = nc_obj.variables[SNOW_NAME].shape
        if actual != expected:
            raise ValueError(f"Output shape is {actual}; expected {expected}")
        time_var = nc_obj.variables["time"]
        values = np.asarray(time_var[:], dtype=float)
        if values.size > 1 and not np.all(np.diff(values) > 0):
            raise ValueError("Output time coordinate is not strictly increasing")
        dates = num2date(
            values,
            units=time_var.units,
            calendar=time_var.calendar,
            only_use_cftime_datetimes=True,
        )
        actual_months = [(item.year, item.month) for item in dates]
        expected_month_values = [
            (item.month.year, item.month.month) for item in granules
        ]
        if actual_months != expected_month_values:
            raise ValueError("Output months do not match sorted input granules")
    print(
        f"[OK] Output validated: {path}\n"
        f"     time={expected[0]}, lat={expected[1]}, lon={expected[2]}"
    )


def validate_all_inputs(granules: list[Granule]) -> tuple[int, int]:
    shape = None
    for index, item in enumerate(granules):
        current = inspect_hdf(item.path)
        if shape is None:
            shape = current
        elif current != shape:
            raise ValueError(f"Grid mismatch in {item.path}: {current} vs {shape}")
        if index == 0 or (index + 1) % 25 == 0 or index + 1 == len(granules):
            print(f"[CHECK] {index + 1:3d}/{len(granules)} {item.path.name}")
    return shape


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        require_dependencies()
        granules = discover_granules(args)
        missing = report_coverage(granules, args)
        first_shape = inspect_hdf(granules[0].path)
        subset = build_subset(first_shape, args.bbox)
        print(
            f"[INFO] Subset: lat={len(subset.latitude)} "
            f"({subset.latitude[0]:.3f}..{subset.latitude[-1]:.3f}), "
            f"lon={len(subset.longitude)} "
            f"({subset.longitude[0]:.3f}..{subset.longitude[-1]:.3f})"
        )
        if args.check_only:
            validate_all_inputs(granules)
            print("[OK] All HDF4 inputs passed structural validation")
            return 0

        args.output.parent.mkdir(parents=True, exist_ok=True)
        if args.output.exists() and not args.overwrite:
            raise FileExistsError(
                f"Output exists: {args.output}. Use --overwrite to replace it."
            )
        # Keep the temporary file outside the ``*.nc*`` glob consumed by the
        # plotting scripts, so an in-progress conversion cannot be evaluated.
        partial = args.output.with_name(f".{args.output.stem}.part")
        if partial.exists():
            if not args.overwrite:
                raise FileExistsError(
                    f"Partial output exists: {partial}. Inspect it or use --overwrite."
                )
            partial.unlink()
        try:
            create_output(partial, granules, subset, missing, args)
            validate_output(partial, granules, subset)
            partial.replace(args.output)
        except Exception:
            if partial.exists():
                partial.unlink()
            raise
        print(f"[DONE] Merged NetCDF: {args.output}")
        return 0
    except (
        FileNotFoundError,
        FileExistsError,
        HDF4Error,
        KeyError,
        RuntimeError,
        ValueError,
    ) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
