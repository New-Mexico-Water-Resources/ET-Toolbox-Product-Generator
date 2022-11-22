"""
This module calculates the Terra/Aqua orbit and segments VIIRS swaths.
"""
import logging
from datetime import timedelta, datetime, date
from typing import Union
from dateutil import parser
import pytz
from raster import Polygon, Point
import numpy as np
import geopandas as gpd
import cl
from .EOS_orbits import SpaceTrackCredentials, SwathGroundTrack, Sensor, OrbitHistory, Satellite, TLEDatabase

__author__ = "Gregory Halverson"

SEARCH_RADIUS_BUFFER_MINUTES = 180
DEFAULT_VIIRS_SWATH_SAMPLE_SIZE = 2

DEFAULT_SEARCH_RADIUS_HOURS = 12

# this determines the number of points located along a swath ground track
DEFAULT_SWATH_SAMPLE_SIZE = 2

# making the cross-track swath width of VIIRS constant at 2330 km for now
# haven't been able to dynamically calculate swath width from ephemerally calculated elevation so far
VIIRS_CROSS_TRACK_SWATH_WIDTH_METERS = 2330000

# VIIRS swaths are always five minutes of orbit
VIIRS_SWATH_DURATION_MINUTES = 5

NPP_SATELLITE = Satellite('NPP', 37849)

DRIFT_CORRECTION = 0

logger = logging.getLogger(__name__)

def longitude_to_offset(longitude_deg):
    return timedelta(hours=(np.radians(longitude_deg) / np.pi * 12))


def UTC_to_solar(datetime_UTC, longitude_deg):
    if isinstance(datetime_UTC, str):
        datetime_UTC = parser.parse(datetime_UTC)

    return datetime_UTC + longitude_to_offset(longitude_deg)


def solar_to_UTC(datetime_solar, longitude_deg):
    if isinstance(datetime_solar, str):
        datetime_solar = parser.parse(datetime_solar)

    longitude_deg = float(longitude_deg)

    datetime_UTC = pytz.UTC.localize(datetime_solar.replace(tzinfo=None) - longitude_to_offset(longitude_deg))

    return datetime_UTC


def get_longitudes(geometry):
    if isinstance(geometry, Polygon):
        return list(geometry.exterior.coords.xy[0])
    if isinstance(geometry, Point):
        return [geometry.x]

class VIIRSOrbit:
    """
    class encapsulating the activity of the VIIRS sensor on both the Aqua and Terra satellites
    """

    def __init__(self, satellite=NPP_SATELLITE, TLE_database=None, lock=None):
        """
        This class encapsulates the activity of the VIIRS sensor on either the Aqua or Terra satellites.
        :param satellite:
            This is an object of type eos_orbits.Satellite
            that describes the satellite the VIIRS intrument is attached to.
            Two constants, TERRA_SATELLITE and AQUA_SATELLITE are provided to support VIIRS calculations.
        :param TLE_database:
            This is an object of type eos_orbits.TLEDatabase that connects to the spacetrack repository
            of orbital parameters.
        """
        self.logger = logging.getLogger(__name__)
        self.lock = lock

        if satellite is None:
            satellite = NPP_SATELLITE

        # validate satellite: only NPP is supported
        if satellite != NPP_SATELLITE:
            raise ValueError('Satellite used for VIIRS orbit must be NPP.')

        # initialize satellite attribute
        self.satellite = satellite

        # initialize sensor attribute from satellite, cross-track swath width, and swath duration
        self.sensor = Sensor(
            'VIIRS',
            satellite,
            cross_track_swath_width_meters=VIIRS_CROSS_TRACK_SWATH_WIDTH_METERS,
            swath_duration_minutes=VIIRS_SWATH_DURATION_MINUTES
        )

        if TLE_database is None:
            TLE_database = TLEDatabase(lock=lock)

        # initialize TLE database attribute
        self.TLE_database = TLE_database
        self.TLE_database.check_database(self.satellite)

    def orbit_history(
            self,
            start_UTC,
            end_UTC,
            sample_size=DEFAULT_SWATH_SAMPLE_SIZE,
            swath_duration=VIIRS_SWATH_DURATION_MINUTES):
        # attach sensor to empty orbit history object to store results
        orbit_history = OrbitHistory(sensor=self.sensor)

        start_time = datetime(
            start_UTC.year,
            start_UTC.month,
            start_UTC.day,
            start_UTC.hour,
            5 * (start_UTC.minute // 5),
            0
        )

        # round end time_handling up to end of last swath
        end_time = datetime(
            end_UTC.year,
            end_UTC.month,
            end_UTC.day,
            end_UTC.hour,
            5 * (end_UTC.minute // 5),
            0
        ) + timedelta(minutes=5)

        logger.info(f"generating VIIRS orbit history from {cl.time(start_UTC)} UTC to {cl.time(end_UTC)} UTC")

        orbit_history_duration_minutes = int((end_time - start_time).total_seconds() / 60)

        # increment in steps of swath_duration minutes
        for swath_start_offset_minutes in range(0, orbit_history_duration_minutes, swath_duration):
            swath_time = start_time + timedelta(minutes=swath_start_offset_minutes)
            swath_start = swath_time - timedelta(minutes=DRIFT_CORRECTION)
            swath_end = swath_start + timedelta(minutes=swath_duration)

            # generate start and end timestamps
            timestamp_start = swath_start.strftime('%Y-%m-%d %H:%M:%S')
            timestamp_end = swath_end.strftime('%Y-%m-%d %H:%M:%S')

            # generate swath name from hour and minutes
            swath_name = '%02d%02d' % (swath_time.hour, swath_time.minute)

            # calculate swath ground track
            swath_ground_track = SwathGroundTrack(
                self.satellite,
                swath_start,
                swath_end,
                self.TLE_database,
                sample_size=sample_size,
                additional_properties={
                    'name': swath_name,
                    'time_start': timestamp_start,
                    'time_end': timestamp_end,
                },
                lock=self.lock
            )

            # append swath to orbit history
            orbit_history.add_swath(swath_ground_track)

        return orbit_history

    def temporal_search_radius(
            self,
            datetime_UTC,
            radius_hours=DEFAULT_SEARCH_RADIUS_HOURS,
            sample_size=DEFAULT_SWATH_SAMPLE_SIZE):

        start_datetime_UTC = datetime_UTC - timedelta(hours=radius_hours)
        end_datetime_UTC = datetime_UTC + timedelta(hours=radius_hours)

        orbit_history = self.orbit_history(
            start_datetime_UTC,
            end_datetime_UTC,
            sample_size=sample_size
        )

        return orbit_history

    def intersect_polygon_and_time(
            self,
            target_polygon,
            datetime_UTC,
            search_radius_hours=DEFAULT_SEARCH_RADIUS_HOURS,
            sample_size=DEFAULT_SWATH_SAMPLE_SIZE,
            use_split_geometry=False,
            day_filter=True):

        swaths = self.temporal_search_radius(
            datetime_UTC,
            search_radius_hours,
            sample_size=sample_size
        )

        # print(len(swaths))

        if day_filter:
            swaths = swaths.day_scenes()

        # print(len(swaths))

        swaths = swaths.intersect_polygon(target_polygon, use_split_geometry=use_split_geometry)

        return swaths

    def swath_polygon(
            self,
            observation_datetime_UTC,
            sample_size=DEFAULT_SWATH_SAMPLE_SIZE,
            additional_properties=None,
            use_split_geometry=False,
            lock=None):
        logger = logging.getLogger(__name__)

        if additional_properties is None:
            additional_properties = {}

        # round start time to start time of swath
        swath_time = datetime(
            observation_datetime_UTC.year,
            observation_datetime_UTC.month,
            observation_datetime_UTC.day,
            observation_datetime_UTC.hour,
            5 * (observation_datetime_UTC.minute // 5),
            0
        )

        swath_start = swath_time - timedelta(minutes=DRIFT_CORRECTION)

        swath_end = swath_start + timedelta(minutes=5)

        additional_properties.update({
            'name': swath_time.strftime('%H%M'),
            'time_start': str(swath_start),
            'time_end': str(swath_end),
        })

        # calculate swath ground track
        swath_ground_track = SwathGroundTrack(
            self.satellite,
            swath_start,
            swath_end,
            self.TLE_database,
            sample_size=sample_size,
            additional_properties=additional_properties,
            lock=lock
        )

        polygon = swath_ground_track.polygon(self.sensor, use_split_geometry=use_split_geometry)

        return polygon


def find_VIIRS_swaths(
        geometry: Polygon,
        date_solar: Union[date, str],
        spacetrack_credentials: SpaceTrackCredentials = None,
        use_split_geometry: bool = True) -> gpd.GeoDataFrame:
    """
    Spatio-temporally intersect target spatial extent and observation time_handling
    with VIIRS orbital swaths for NPP.
    """

    # create empty region identifier list
    region_identifier_list = []

    longitudes = get_longitudes(geometry)

    minimum_longitude = min(longitudes)
    maximum_longitude = max(longitudes)

    search_radius = timedelta(
        hours=(np.radians(maximum_longitude - minimum_longitude) / np.pi * 12)).total_seconds() / 3600.0 / 2.0 + (
                            SEARCH_RADIUS_BUFFER_MINUTES / 60.0)

    TLE_database = TLEDatabase(spacetrack_credentials)
    VIIRS_orbit = VIIRSOrbit(TLE_database=TLE_database)

    centroid_longitude = Polygon(geometry).centroid.latlon.x

    if isinstance(date_solar, str):
        date_solar = parser.parse(date_solar).date()
    
    datetime_solar = datetime(date_solar.year, date_solar.month, date_solar.day, 10, 30, 0)
    datetime_UTC = solar_to_UTC(datetime_solar, centroid_longitude)

    intersecting_swath_orbit_history = VIIRS_orbit.intersect_polygon_and_time(
        geometry.buffer(1.0),
        datetime_UTC,
        search_radius_hours=search_radius,
        sample_size=DEFAULT_VIIRS_SWATH_SAMPLE_SIZE,
        use_split_geometry=use_split_geometry
    )
    
    times = []
    names = []
    geometries = []

    for swath in intersecting_swath_orbit_history.swath_ground_tracks:
        times.append(swath.start_datetime)
        names.append(swath.start_datetime.strftime("%H%M"))
        geometries.append(swath._polygon)
    
    gdf = gpd.GeoDataFrame({"time": times, "name": names}, geometry=geometries)

    return gdf
