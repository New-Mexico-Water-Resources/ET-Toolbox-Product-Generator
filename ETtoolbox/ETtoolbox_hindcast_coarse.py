import logging
import sys
from datetime import datetime, timedelta, date
from os.path import join, abspath, dirname, exists
from typing import List, Callable, Union

import cl
from GEDI import GEDICanopyHeight
from GEOS5FP import GEOS5FP
from LANCE_GEOS5FP_NRT import LANCE_GEOS5FP_NRT, LANCENotAvailableError, GEOS5FPNotAvailableError
from ORNL.MODISCI import MODISCI
from PTJPLSM import PTJPLSM, DEFAULT_PREVIEW_QUALITY, DEFAULT_RESAMPLING
from SRTM import SRTM
from SoilGrids import SoilGrids
from raster import Raster, RasterGrid
from sentinel import sentinel_tile_grid

logger = logging.getLogger(__name__)

DEFAULT_MODEL_NAME = "PTJPLSM"
DEFAULT_DOWNSCALE_AIR = True
DEFAULT_DOWNSCALE_HUMIDITY = False
DEFAULT_DOWNSCALE_MOISTURE = True
DEFAULT_MESO_CELL_SIZE = 500
DEFAULT_COARSE_CELL_SIZE = 27375

DEFAULT_TARGET_VARIABLES = [
    "Rn",
    "LE",
    "ETc",
    "ETi",
    "ETs",
    "ET",
    "ESI",
    "WUE",
    "SM",
    "Ta",
    "RH"
]


def ET_toolbox_hindcast_coarse_tile(
        tile: str,
        geometry: RasterGrid = None,
        present_date: Union[date, str] = None,
        water: Raster = None,
        model: PTJPLSM = None,
        model_name: str = "PTJPLSM",
        working_directory: str = None,
        static_directory: str = None,
        LANCE_download: str = None,
        output_directory: str = None,
        SRTM_connection: SRTM = None,
        SRTM_download: str = None,
        GEOS5FP_connection: GEOS5FP = None,
        GEOS5FP_download: str = None,
        GEOS5FP_products: str = None,
        GEDI_connection: GEDICanopyHeight = None,
        GEDI_download: str = None,
        ORNL_connection: MODISCI = None,
        CI_directory: str = None,
        soil_grids_connection: SoilGrids = None,
        soil_grids_download: str = None,
        intermediate_directory: str = None,
        preview_quality: int = DEFAULT_PREVIEW_QUALITY,
        ANN_model: Callable = None,
        ANN_model_filename: str = None,
        spacetrack_credentials_filename: str = None,
        resampling: str = DEFAULT_RESAMPLING,
        meso_cell_size: float = DEFAULT_MESO_CELL_SIZE,
        coarse_cell_size: float = DEFAULT_COARSE_CELL_SIZE,
        downscale_air: bool = DEFAULT_DOWNSCALE_AIR,
        downscale_humidity: bool = DEFAULT_DOWNSCALE_HUMIDITY,
        downscale_moisture: bool = DEFAULT_DOWNSCALE_MOISTURE,
        save_intermediate: bool = False,
        show_distribution: bool = True,
        load_previous: bool = True,
        target_variables: List[str] = None):
    if present_date is None:
        present_date = datetime.utcnow().date()

    logger.info(
        f"generating ET Toolbox hindcast and forecast at tile {cl.place(tile)} centered on present date: {cl.time(present_date)}")

    if geometry is None:
        geometry = sentinel_tile_grid.grid(tile, cell_size=meso_cell_size)

    if target_variables is None:
        target_variables = DEFAULT_TARGET_VARIABLES

    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP(
            working_directory=working_directory,
            download_directory=GEOS5FP_download,
            products_directory=GEOS5FP_products
        )

    if SRTM_connection is None:
        # FIXME fix handling of credentials here
        SRTM_connection = SRTM(
            working_directory=working_directory,
            download_directory=SRTM_download,
            offline_ok=True
        )

    bundled_spacetrack_credentials_filename = join(abspath(dirname(__file__)), "spacetrack")

    if spacetrack_credentials_filename is None and exists(bundled_spacetrack_credentials_filename):
        spacetrack_credentials_filename = bundled_spacetrack_credentials_filename

    missing_dates = []

    for relative_days in range(-7, 1):
        target_date = present_date + timedelta(days=relative_days)
        logger.info(f"LANCE GEOS-5 FP target date: {cl.time(target_date)} ({cl.time(relative_days)} days)")

        time_solar = datetime(target_date.year, target_date.month, target_date.day, 13, 30)
        logger.info(f"LANCE target time solar: {cl.time(time_solar)}")

        try:
            LANCE_GEOS5FP_NRT(
                target_date=target_date,
                geometry=geometry,
                target=tile,
                working_directory=working_directory,
                static_directory=static_directory,
                SRTM_connection=SRTM_connection,
                SRTM_download=SRTM_download,
                GEOS5FP_connection=GEOS5FP_connection,
                GEOS5FP_download=GEOS5FP_download,
                GEOS5FP_products=GEOS5FP_products,
                GEDI_connection=GEDI_connection,
                GEDI_download=GEDI_download,
                ORNL_connection=ORNL_connection,
                CI_directory=CI_directory,
                soil_grids_connection=soil_grids_connection,
                soil_grids_download=soil_grids_download,
                LANCE_download_directory=LANCE_download,
                LANCE_output_directory=output_directory,
                intermediate_directory=intermediate_directory,
                preview_quality=preview_quality,
                ANN_model=ANN_model,
                ANN_model_filename=ANN_model_filename,
                model=model,
                model_name=model_name,
                water=water,
                coarse_cell_size=coarse_cell_size,
                target_variables=target_variables,
                downscale_air=downscale_air,
                downscale_humidity=downscale_humidity,
                downscale_moisture=downscale_moisture,
                resampling=resampling,
                show_distribution=show_distribution,
                load_previous=load_previous,
                save_intermediate=save_intermediate,
                spacetrack_credentials_filename=spacetrack_credentials_filename
            )
        except (LANCENotAvailableError, GEOS5FPNotAvailableError) as e:
            logger.warning(e)
            logger.warning(f"LANCE GEOS-5 FP cannot be processed for date: {target_date}")
            missing_dates.append(target_date)
            continue
        except Exception as e:
            logger.exception(e)
            logger.warning(f"LANCE GEOS-5 FP cannot be processed for date: {target_date}")
            missing_dates.append(target_date)
            continue

    logger.info("missing LANCE GEOS-5 FP dates: " + ", ".join(cl.time(d) for d in missing_dates))


def main(argv=sys.argv):
    tile = argv[1]

    if "--working" in argv:
        working_directory = argv[argv.index("--working") + 1]
    else:
        working_directory = "."

    if "--static" in argv:
        static_directory = argv[argv.index("--static") + 1]
    else:
        static_directory = join(working_directory, "PTJPL_static")

    if "--SRTM" in argv:
        SRTM_download = argv[argv.index("--SRTM") + 1]
    else:
        SRTM_download = join(working_directory, "SRTM_download_directory")

    if "--LANCE" in argv:
        LANCE_download_directory = argv[argv.index("--LANCE") + 1]
    else:
        LANCE_download_directory = join(working_directory, "LANCE_download_directory")

    if "--GEOS5FP" in argv:
        GEOS5FP_download = argv[argv.index("--GEOS5FP") + 1]
    else:
        GEOS5FP_download = join(working_directory, "GEOS5FP_download_directory")

    if "--output" in argv:
        output_directory = argv[argv.index("--output") + 1]
    else:
        output_directory = join(working_directory, "output_directory")

    if "--spacetrack" in argv:
        spacetrack_credentials_filename = argv[argv.index("--spacetrack") + 1]
    else:
        spacetrack_credentials_filename = None

    ET_toolbox_hindcast_coarse_tile(
        tile=tile,
        working_directory=working_directory,
        static_directory=static_directory,
        output_directory=output_directory,
        SRTM_download=SRTM_download,
        LANCE_download=LANCE_download_directory,
        GEOS5FP_download=GEOS5FP_download,
        spacetrack_credentials_filename=spacetrack_credentials_filename
    )
