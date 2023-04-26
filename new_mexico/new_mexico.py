import sys
from datetime import datetime
from os.path import join
from time import sleep

from ETtoolbox.ETtoolbox_hindcast_coarse import ET_toolbox_hindcast_coarse_tile
import logging
import cl

logger = logging.getLogger(__name__)

def new_mexico_VIIRS_server(
        working_directory: str = None,
        static_directory: str = None,
        LANCE_download: str = None,
        SRTM_download: str = None,
        GEOS5FP_download: str = None):
    if working_directory is None:
        working_directory = "/new_mexico_VIIRS"

    if static_directory is None:
        static_directory = "/static"

    if LANCE_download is None:
        LANCE_download = "/LANCE"

    if SRTM_download is None:
        SRTM_download = "/SRTM"

    if GEOS5FP_download is None:
        GEOS5FP_download = "/GEOS5FP"

    logger.info("starting New Mexico VIIRS data production")

    logger.info(f"working directory: {working_directory}")
    logger.info(f"static directory: {static_directory}")
    logger.info(f"LANCE directory: {LANCE_download}")
    logger.info(f"SRTM directory: {SRTM_download}")
    logger.info(f"GEOS-5 FP directory: {GEOS5FP_download}")

    tiles = ["h08v05", "h09v05"]

    while(True):
        runtime = datetime.utcnow()
        logger.info(f"running New Nexico VIIRS data production at time {runtime} UTC")

        for tile in tiles:
            ET_toolbox_hindcast_coarse_tile(
                tile=tile,
                working_directory=working_directory,
                static_directory=static_directory,
                SRTM_download=SRTM_download,
                LANCE_download=LANCE_download,
                GEOS5FP_download=GEOS5FP_download,
            )

        while(datetime.utcnow().hour % 3 != 0):
            sleep(60)


def main(argv=sys.argv):
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

    return new_mexico_VIIRS_server()


if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
