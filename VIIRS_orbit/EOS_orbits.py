"""
This package calculates satellite positions.
"""
import base64
import itertools
import json
import logging
import sqlite3
import warnings
from datetime import date, datetime, timedelta
from itertools import tee
from os import makedirs
from os.path import exists, join, dirname, abspath
from time import sleep

import ephem
import fiona
import numpy
import shapely
import shapely.wkt
from dateutil import parser as dateparser, parser
from fiona.crs import from_string
from numpy.core.umath import arctan2, pi, cos, sin, tan, arccos, radians, degrees
from pyproj import Proj, Transformer
from shapely.geometry import mapping, Point, LineString, LinearRing, Polygon, shape
from shapely.ops import transform as shapely_transform
from six import string_types
from six.moves import zip
from spacetrack import SpaceTrackClient

__author__ = "Gregory Halverson"

from VIIRS_orbit.spacetrack_credentials import get_spacetrack_credentials
from transform import transform_point

DATABASE_UP_TO_DATE_THRESHOLD = timedelta(days=2)

DEFAULT_LINE_PROPERTIES = (
    'node', 'bearing', 'lat_start', 'lon_start', 'elev_start', 'lat_end', 'lon_end', 'elev_end', 'daynight')

DEFAULT_SHAPEFILE_PROPERTIES = [
    ('name', 'str'),
    ('time_start', 'str'),
    ('time_end', 'str'),
    ('node', 'str'),
    ('bearing', 'float'),
    ('lat_start', 'float'),
    ('lon_start', 'float'),
    ('elev_start', 'float'),
    ('lat_end', 'float'),
    ('lon_end', 'float'),
    ('elev_end', 'float'),
    ('daynight', 'str')
]

SHAPEFILE_DRIVER_NAME = 'ESRI Shapefile'

WGS84_PROJ4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

WGS84 = WGS84_PROJ4

WGS84_FIONA_CRS = fiona.crs.from_string(WGS84_PROJ4)

DATABASE_DIRECTORY_NAME = 'eos_orbits'
DEFAULT_SQLITE3_DB_FILENAME = 'tle.sqlite3'
DEFAULT_SWATH_DB_FILENAME = 'swath.sqlite3'

SQLITE3_TIMEOUT_SECONDS = 30

logger = logging.getLogger(__name__)

class DummyLock:
    def __enter__(self):
        pass

    def __exit__(self, exception_type, exception_value, exception_traceback):
        pass


# calculate compass bearing from one lat/lon coordinate to another
def bearing(coord_start, coord_end):
    """
    Calculates compass bearing from one lat/lon coordinate to another.
    :param coord_start:
        starting point of bearing given as shapely point
        with degrees longitude as x and degrees latitude as y
        (shapely.geometry.Point)
    :param coord_end:
        end point of bearing given as shapely point
        with degrees longitude as x and degrees latitude as y
        (shapely.geometry.Point)
    :return:
        compass bearing from start point to end point in degrees
        clockwise from north
        (float)
    """

    # convert starting latitude from degrees to radians
    lat_start = radians(coord_start.y)

    # convert ending latitude from degrees to radians
    lat_end = radians(coord_end.y)

    # calculate change in longitude and convert from degrees to radians
    lon_delta = radians(coord_end.x - coord_start.x)

    # calculate bearing and normalize degrees
    bearing = (degrees(arctan2(
        sin(lon_delta) * cos(lat_end),
        cos(lat_start) * sin(lat_end) - (sin(lat_start) * cos(lat_end) * cos(lon_delta))
    )) + 360) % 360

    return bearing


def center_aeqd_proj(center_coord: Point) -> Proj:
    return Proj(center_aeqd_proj4(center_coord))


def center_aeqd_proj4(center_coord: Point) -> str:
    return '+proj=aeqd +lat_0=%f +lon_0=%f' % (
        center_coord.y,
        center_coord.x
    )


def centered_aeqd_coord(target_coord, center_coord):
    return transform_point(target_coord, WGS84, center_aeqd_proj4(center_coord))


# calculate a series of times
def time_series(start_datetime, end_datetime, sample_size):
    """
    Calculate a series of times.
    :param start_datetime: Start time_handling as datetime object.
    :param end_datetime: End time_handling as datetime object.
    :param sample_size: Number of times in time_handling series including start and end times.
    :return: list of sample_size datetime objects starting at start_datetime and ending at end_datetime
    """
    # calculate duration between start and end times
    duration = end_datetime - start_datetime

    # calculate timedelta between elements of the series
    delta = duration / (sample_size - 1)

    # calculate series
    series = [
        start_datetime + delta * i
        for i
        in range(sample_size)
    ]

    return series


def split_geometry(geometry):
    def to_polar(lon, lat):
        phi = numpy.pi / 180. * lon
        radius = numpy.pi / 180. * (90. - sign * lat)

        # nudge points at +/- 180 out of the way so they don't intersect the testing wedge
        phi = numpy.sign(phi) * numpy.where(numpy.abs(phi) > numpy.pi - 1.5 * epsilon, numpy.pi - 1.5 * epsilon,
                                            numpy.abs(phi))
        radius = numpy.where(radius < 1.5 * epsilon, 1.5 * epsilon, radius)

        x = radius * numpy.sin(phi)
        y = radius * numpy.cos(phi)
        if (isinstance(lon, list)):
            x = x.tolist()
            y = y.tolist()
        elif (isinstance(lon, tuple)):
            x = tuple(x)
            y = tuple(y)

        return (x, y)

    def from_polar(x, y):
        radius = numpy.sqrt(numpy.array(x) ** 2 + numpy.array(y) ** 2)
        phi = numpy.arctan2(x, y)

        # close up the tiny gap
        radius = numpy.where(radius < 2 * epsilon, 0., radius)
        phi = numpy.sign(phi) * numpy.where(numpy.abs(phi) > numpy.pi - 2 * epsilon, numpy.pi, numpy.abs(phi))

        lon = 180. / numpy.pi * phi
        lat = sign * (90. - 180. / numpy.pi * radius)

        if (isinstance(x, list)):
            lon = lon.tolist()
            lat = lat.tolist()
        elif (isinstance(x, tuple)):
            lon = tuple(lon)
            lat = tuple(lat)
        return (lon, lat)

    epsilon = 1e-14

    antimeridian_wedge = shapely.geometry.Polygon([
        (epsilon, -numpy.pi),
        (epsilon ** 2, -epsilon),
        (0, epsilon),
        (-epsilon ** 2, -epsilon),
        (-epsilon, -numpy.pi),
        (epsilon, -numpy.pi)
    ])

    feature_shape = shapely.geometry.shape(geometry)
    sign = 2. * (0.5 * (feature_shape.bounds[1] + feature_shape.bounds[3]) >= 0.) - 1.

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        polar_shape = shapely_transform(to_polar, feature_shape)

    if not polar_shape.intersects(antimeridian_wedge):
        return geometry

    else:
        pass

    difference = polar_shape.difference(antimeridian_wedge)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        output_shape = shapely_transform(from_polar, difference)

    return output_shape


# class holding the name and ID of a satellite
class Satellite:
    def __init__(self, satellite_name, norad_id):
        """
        Name and ID of a satellite.
        :param satellite_name: Name of satellite used by spacetrack as string.
        :param norad_id: NORAD ID of satellite as integer.
        :return:
        """
        self.satellite_name = str(satellite_name)
        self.norad_id = int(norad_id)

    def __str__(self):
        return self.satellite_name

    def __repr__(self):
        return "<Satellite name: {} ID: {}>".format(self.satellite_name, self.norad_id)


class SpaceTrackCredentials:
    def __init__(self, username, password):
        self.username = str(username)
        self.password = str(password)


# class to hold Two-Line Element orbital parameters
class TLE:
    def __init__(self, epoch_datetime, tle_text):
        """
        Two-Line Element orbital parameters.
        :param epoch_datetime: Epoch time_handling of TLE as datetime object.
        :param tle_text: Three-line TLE as single string.
        """
        if isinstance(epoch_datetime, datetime):
            self.epoch_datetime = epoch_datetime
        elif isinstance(epoch_datetime, string_types):
            self.epoch_datetime = dateparser.parse(epoch_datetime)
        else:
            raise ValueError('epoch must be datetime object or string')

        self.tle_text = str(tle_text)

    def __str__(self):
        return self.tle_text

    def __repr__(self):
        return self.__str__()

    # epoch time_handling of TLE
    def epoch(self):
        """
        :return: Epoch time_handling of Two-Line Element.
        """
        return self.epoch_datetime

    # TLE text split into three lines
    def as_three_lines(self):
        """
        :return: Two-Line Element text split into three lines.
        """
        return str(self.tle_text).split('\n')


# local cache of spacetrack database
class TLEDatabase:
    def __init__(
            self,
            spacetrack_credentials=None,
            sqlite3_db_location=DEFAULT_SQLITE3_DB_FILENAME,
            lock=None):
        """
        Local cache of spacetrack database.
        :param sqlite3_db_location:
            Name of local sqlite3 database file.
        :param spacetrack_credentials:
            Object of class SpaceTrackCredentials containing username and password to access spacetrack TLEs.
        """

        if lock is None:
            lock = DummyLock()

        self.lock = lock

        if sqlite3_db_location is None:
            self.db_location = DEFAULT_SQLITE3_DB_FILENAME
        else:
            self.db_location = sqlite3_db_location

        if spacetrack_credentials is None:
            credentials = get_spacetrack_credentials()
            username = credentials["username"]
            password = credentials["password"]
            spacetrack_credentials = SpaceTrackCredentials(username, password)

        self.spacetrack_credentials = spacetrack_credentials

    def _parse_tle_query_text_to_records(self, tle_query_text):
        lines = tle_query_text.split('\n')
        records = []

        for one, two, three in zip(*[itertools.islice(lines, i, None, 3) for i in range(3)]):
            date_text = two.split()[3]
            year = 2000 + int(date_text[:2])
            doy = int(date_text[2:5])
            d = date(year, 1, 1) + timedelta(days=(doy - 1))
            day_fraction = float(date_text[5:])
            dt = datetime(d.year, d.month, d.day) + timedelta(days=day_fraction)
            tle = '\n'.join((one, two, three))
            records.append((dt, tle))

        return records

    def _query_tle_record(self, satellite, target_datetime):
        TRIES = 3
        WAIT_SECONDS = 30

        for i in range(1, TRIES + 1):
            try:
                with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
                    cursor = connection.cursor()
                    target_timestamp = target_datetime.strftime('%Y-%m-%d %H:%M:%S')

                    # generate SQL string
                    tle_query_sql_string = 'SELECT * FROM %s WHERE dt <= ? ORDER BY dt DESC LIMIT 1 ' % (
                        satellite.satellite_name
                    )

                    # query database
                    records = list(
                        cursor.execute(
                            tle_query_sql_string,
                            [
                                target_timestamp
                            ]
                        )
                    )

                    # check if record exists
                    if len(records) == 0:
                        record = None
                    else:
                        record = records[0]

                return record
            except:
                logger.warning(
                    "access to db at {} failed on attempt {}, waiting {} seconds".format(
                        self.db_location,
                        i,
                        WAIT_SECONDS
                    ))

                sleep(WAIT_SECONDS)

                continue

            raise IOError("all attempts to access db at {} failed".format(self.db_location))

    @property
    def database_exists(self):
        return exists(self.db_location)

    def table_exists(self, satellite):
        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()

            # generate table existence SQL query
            table_existence_sql_string = "SELECT * FROM sqlite_master WHERE type='table' AND name='%s'" % (
                satellite.satellite_name
            )

            result = len(list(cursor.execute(table_existence_sql_string))) != 0

        return result

    def update_tle_table(
            self,
            satellite):

        last_access_filename = join(dirname(self.db_location), "tle_last_access.txt")
        now = datetime.utcnow()

        if exists(last_access_filename):
            try:
                with open(last_access_filename, "r") as f:
                    last = parser.parse(f.read())

                delta = now - last
                delta_seconds = delta.total_seconds()

                logger.info("seconds since last access to Spacetrack: {}".format(delta_seconds))

                if delta_seconds < 21600:
                    logger.warning("preventing premature access to Spacetrack")
                    return
            except:
                pass

        now_timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

        logger.info("accessing Spacetrack at {} UTC".format(now_timestamp))

        directory = dirname(last_access_filename)

        if not exists(directory) and directory != "":
            makedirs(directory)

        with open(last_access_filename, "w") as f:
            f.write(now_timestamp)

        # generate table existence SQL query
        table_existence_sql_string = "SELECT * FROM sqlite_master WHERE type='table' AND name='%s'" % (
            satellite.satellite_name
        )

        # generate table create SQL query
        table_creation_sql_string = 'CREATE TABLE %s (dt TEXT UNIQUE, tle TEXT)' % (
            satellite.satellite_name
        )

        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()

            # create table for satellite if it doesn't exist
            if (len(list(cursor.execute(table_existence_sql_string))) == 0):
                try:
                    cursor.execute(table_creation_sql_string)
                except sqlite3.OperationalError as e:
                    pass

        # connect to Space Track
        st = SpaceTrackClient(
            identity=self.spacetrack_credentials.username,
            password=self.spacetrack_credentials.password
        )

        # query all TLEs for satellite
        tle_query_text = st.tle(
            norad_cat_id=satellite.norad_id,
            epoch='<%s' % (datetime.utcnow().date() + timedelta(1)).strftime('%Y-%m-%d'),
            orderby='epoch',
            format='3le'
        )

        # parse retrieved TLE
        records = self._parse_tle_query_text_to_records(tle_query_text)

        # generate insertion SQL string
        tle_insertion_sql_string = 'INSERT OR REPLACE INTO %s VALUES(?, ?)' % satellite.satellite_name

        try:
            with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
                cursor = connection.cursor()

                # store TLEs in local database
                for dt, tle in records:
                    cursor.execute(tle_insertion_sql_string, (
                        dt.strftime('%Y-%m-%d %H:%M:%S'),
                        tle
                    ))
        except Exception as e:
            logger.error(e)
            raise IOError("unable to connect to database {}".format(self.db_location))

    def record_count(self, satellite):
        sql_query_string = 'SELECT count(*) FROM %s' % satellite.satellite_name

        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()
            cursor.execute(sql_query_string)

    def most_recent_entry(self, satellite):
        if not self.database_exists:
            self.update_tle_table(satellite)

        now = datetime.utcnow()
        most_recent = self._query_tle_record(satellite, now)

        if most_recent is None:
            return None
        else:
            return (TLE(*most_recent))

    def age_of_last_entry(self, satellite):
        most_recent = self.most_recent_entry(satellite)

        if most_recent is None:
            self.update_tle_table(satellite)
            most_recent = self.most_recent_entry(satellite)

            if most_recent is None:
                raise Exception('cannot access most recent TLE')

        age = datetime.utcnow() - most_recent.epoch()

        return age

    def check_database(self, satellite):
        if not self.database_exists:
            self.update_tle_table(satellite)
            return

        # logger.info("checking if table for {} exists".format(satellite))

        if not self.table_exists(satellite):
            # logger.info('table does not exist, updating')
            self.update_tle_table(satellite)
            return

        logger.info("checking if orbit table is up to date for {}".format(satellite))

        age = self.age_of_last_entry(satellite)

        logger.info("age of last entry: {}".format(age))

        if age > DATABASE_UP_TO_DATE_THRESHOLD:
            # logger.info('table is not up to date, updating')
            self.update_tle_table(satellite)
            return

    def get_tle(self, satellite, target_datetime):
        try:
            record = self._query_tle_record(satellite, target_datetime)
        except:
            self.check_database(satellite)
            record = self._query_tle_record(satellite, target_datetime)

        if record is None:
            print('no record found for {} at {}'.format(satellite, target_datetime))
            return None
        else:
            return TLE(*record)

    def get_satellite_history(self, satellite):
        self.check_database(satellite)

        # generate SQL string
        tle_query_sql_string = 'SELECT * FROM %s' % (
            satellite.satellite_name
        )

        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()
            satellite_history = list(cursor.execute(tle_query_sql_string))

        return satellite_history


class SwathDatabase:
    # TODO cache swath polygons here the same way as TLEs

    def __init__(
            self,
            satellite_name,
            sensor_name,
            swath_db_location=DEFAULT_SWATH_DB_FILENAME,
            lock=None):

        if lock is None:
            lock = DummyLock()

        self.lock = lock

        self.db_location = swath_db_location
        self.satellite_name = str(satellite_name)
        self.sensor_name = str(sensor_name)

    @property
    def database_exists(self):
        return exists(self.db_location)

    @property
    def table_name(self):
        return "{}_{}".format(self.satellite_name, self.sensor_name)

    @property
    def table_exists(self):
        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()

            # generate table existence SQL query
            table_existence_sql_string = "SELECT * FROM sqlite_master WHERE type='table' AND name='%s'" % self.table_name
            result = len(list(cursor.execute(table_existence_sql_string))) != 0

        return result

    def update_table(self):
        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()
            table_name = self.table_name

            # generate table existence SQL query
            table_existence_sql_string = "SELECT * FROM sqlite_master WHERE type='table' AND name='%s'" % table_name

            # generate table create SQL query
            table_creation_sql_string = 'CREATE TABLE %s (dt TEXT UNIQUE, poly TEXT)' % table_name

            # create table for satellite if it doesn't exist
            if (len(list(cursor.execute(table_existence_sql_string))) == 0):
                cursor.execute(table_creation_sql_string)

    def put(self, dt, poly):
        if isinstance(poly, string_types):
            poly = shape(poly)

        if isinstance(dt, string_types):
            dt = parser.parse(dt)

        wkt = poly.wkt
        dt_string = str(dt.strftime('%Y-%m-%d %H:%M:%S'))

        self.update_table()

        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()

            # generate insertion SQL string
            insertion_sql_string = 'INSERT OR REPLACE INTO %s VALUES(?, ?)' % self.table_name

            cursor.execute(insertion_sql_string, (
                dt_string,
                wkt
            ))

    def get(self, dt):
        if not self.table_exists:
            return None

        if isinstance(dt, string_types):
            dt = parser.parse(dt)

        dt_string = str(dt.strftime('%Y-%m-%d %H:%M:%S'))

        with sqlite3.connect(self.db_location, timeout=SQLITE3_TIMEOUT_SECONDS) as connection:
            cursor = connection.cursor()

            # generate SQL string
            query_sql_string = 'SELECT * FROM %s WHERE dt == ? ORDER BY dt DESC LIMIT 1 ' % self.table_name

            # query database
            records = list(cursor.execute(query_sql_string, [dt_string]))

            # check if record exists
            if len(records) == 0:
                record = None
            else:
                record = records[0]

        if record is None:
            return None

        poly_string = str(record[1])
        poly = shapely.wkt.loads(poly_string)

        return poly


# class to calculate the position of a satellite
class SatellitePosition:
    def __init__(self, overpass_datetime, TLE):
        """
        Calculates the position of a satellite.
        :param overpass_datetime: time_handling of overpass as datetime
        :param TLE: tle as TLE object
        """

        # initialize attributes
        self.overpass_datetime = overpass_datetime
        self.TLE = TLE

        # create ephemeris from TLE
        self.ephemeris = ephem.readtle(*TLE.as_three_lines())

        # calculate ephemeris at time_handling of overpass
        self.ephemeris.compute(self.overpass_datetime)

        # calculate latitude at nadir in degrees
        # calculate longitude at nadir in degrees
        # retrieve elevation in meters
        # pack x, y, z coordinate into point object
        self.coordinate = Point(
            round(self.ephemeris.sublong * 57.2958, 6),
            round(self.ephemeris.sublat * 57.2958, 6),
            self.ephemeris.elevation
        )

    # get latitude in degrees
    def latitude(self):
        return self.coordinate.y

    # get longitude in degrees
    def longitude(self):
        return self.coordinate.x

    # get elevation in meters
    def elevation(self):
        return self.coordinate.z

    def day_night(self):
        solar_time = self.overpass_datetime + timedelta(hours=(radians(self.coordinate.x) / pi * 12))
        solar_hour = solar_time.hour + solar_time.minute / 3600.0 + solar_time.second / 86400.0

        day_of_year = self.overpass_datetime.timetuple().tm_yday
        day_angle_rad = (2 * pi * (day_of_year - 1)) / 365
        solar_declination_deg = (0.006918 - 0.399912 * cos(day_angle_rad) + 0.070257 * sin(
            day_angle_rad) - 0.006758 * cos(2 * day_angle_rad) + 0.000907 * sin(2 * day_angle_rad) - 0.002697 * cos(
            3 * day_angle_rad) + 0.00148 * sin(3 * day_angle_rad)) * (180 / pi)
        sunrise_cosine = -tan(pi * self.coordinate.y / 180) * tan(pi * solar_declination_deg / 180)

        if sunrise_cosine >= 1:
            sha = 0
        elif sunrise_cosine <= -1:
            sha = 180
        else:
            sha = arccos(sunrise_cosine) * (180 / pi)

        sunrise_hour = 12 - (sha / 15)
        daylight_hours = (2.0 / 15.0) * sha

        if sunrise_hour < solar_hour < sunrise_hour + daylight_hours:
            return 'day'
        else:
            return 'night'


class FieldOfView:
    def __init__(self, satellite_position, next_position, sensor):
        next_coord_aeqd = centered_aeqd_coord(next_position.coordinate, satellite_position.coordinate)
        direction = arctan2(next_coord_aeqd.y, next_coord_aeqd.x)

        right_direction = direction - pi / 2
        left_direction = direction + pi / 2

        self.cross_track_swath_width = sensor.cross_track_swath_width(satellite_position)

        radius = self.cross_track_swath_width / 2.0

        right_coord_aeqd = Point(radius * cos(right_direction), radius * sin(right_direction))
        left_coord_aeqd = Point(radius * cos(left_direction), radius * sin(left_direction))

        transformer = Transformer.from_crs(center_aeqd_proj4(satellite_position.coordinate), WGS84)
        self.right_coord = transform_point(right_coord_aeqd, transformer=transformer)
        self.left_coord = transform_point(left_coord_aeqd, transformer=transformer)


# class to calculate ground track of satellite swath
class SwathGroundTrack:
    def __init__(
            self,
            satellite,
            start_datetime,
            end_datetime,
            tle_db,
            sample_size=2,
            additional_properties=None,
            lock=None):
        """
        Class to calculate ground track of satellite swath.
        :param satellite:
            satellite name string to query in TLE database
        :param start_datetime:
            start time_handling of swath as datetime
        :param end_datetime:
            end time_handling of swath as datetime
        :param tle_db:
            TLE database as TLEDatabase object
        :param sample_size:
            number of points to calculate along swath
            including start and end point
            defaults to 2 for just start and end
        """

        if lock is None:
            lock = DummyLock()

        self.lock = lock

        if additional_properties is None:
            additional_properties = {}

        # initialize attributes
        self.satellite = satellite
        self.start_datetime = start_datetime
        self.end_datetime = end_datetime
        self.tle_db = tle_db
        self.additional_properties = additional_properties

        # logger.info("calculating satellite positions")

        position_times = time_series(
            self.start_datetime,
            self.end_datetime,
            sample_size
        )

        self.position_list = []

        for overpass_time in position_times:
            TLE = self.tle_db.get_tle(
                self.satellite,
                overpass_time
            )

            if TLE is None:
                continue

            position = SatellitePosition(
                overpass_time,
                TLE
            )

            if position is not None:
                self.position_list.append(position)

        # calculate bearing
        self.bearing = bearing(self.start_position().coordinate, self.end_position().coordinate)

        # round off floating point error from bearing
        self.bearing = round(self.bearing, 6)

        # determine node from bearing
        if self.bearing < 90 or self.bearing > 270:
            self.node = 'ascending'
        else:
            self.node = 'descending'

        self._polygon = None

    def __str__(self):
        return "<SwathGroundTrack start: %s end: %s>" % (
            self.start_datetime,
            self.end_datetime
        )

    # get start position of swath
    def start_position(self):
        return self.position_list[0]

    # get start coordinate of swath
    def start_coord(self):
        return self.start_position().coordinate

    # get start time_handling of swath
    def start_time(self):
        return self.start_position().overpass_datetime

    # get end position of swath
    def end_position(self):
        return self.position_list[-1]

    # get end coordinate of swath
    def end_coord(self):
        return self.end_position().coordinate

    # get end time_handling of swath
    def end_time(self):
        return self.end_position().overpass_datetime

    # get time_handling duration of swath
    def duration(self):
        return self.end_time() - self.start_time()

    # get number of points calculated in swath ground track
    def sample_size(self):
        return len(self.position_list)

    # get time_handling between points in swath ground track
    def time_delta(self):
        return self.duration() / self.sample_size()

    def day_night(self):
        if self.start_position().day_night() == 'day':  # and self.end_position().day_night() == 'day':
            return 'day'
        else:
            return 'night'

    # get list of coordinates in swath ground track
    def points(self):
        return [
            position.coordinate
            for position
            in self.position_list
        ]

    # get swath ground track coordinates as line
    def line(self):
        return LineString(self.points())

    # get swath footprint polygon for sensor on swath track
    def polygon(self, sensor, use_split_geometry=True, use_database=False):

        if self._polygon is not None:
            return self._polygon

        if use_database:
            swath_database = SwathDatabase(
                satellite_name=self.satellite,
                sensor_name=sensor,
                lock=self.lock
            )

            polygon = swath_database.get(self.start_time())

            if polygon is not None:
                self._polygon = polygon
                return self._polygon

        last_time = self.end_time() + self.time_delta()
        last_tle = self.tle_db.get_tle(self.satellite, last_time)
        last_position = SatellitePosition(last_time, last_tle)

        a, b = tee(self.position_list + [last_position])
        next(b, None)
        position_list = zip(a, b)

        fov_list = [
            FieldOfView(position, next_position, sensor)
            for position, next_position
            in position_list
        ]

        coordinates = [
                          (self.start_coord().x, self.start_coord().y)
                      ] + [
                          (fov.left_coord.x, fov.left_coord.y)
                          for fov
                          in fov_list
                      ] + [
                          (self.end_coord().x, self.end_coord().y)
                      ] + [
                          (fov.right_coord.x, fov.right_coord.y)
                          for fov
                          in list(reversed(fov_list))
                      ]

        polygon = Polygon(LinearRing(coordinates))

        if use_split_geometry:
            polygon = split_geometry(polygon)

        if use_database:
            swath_database.put(self.start_time(), polygon)

        self._polygon = polygon

        return self._polygon

    def add_properties(self, properties_dict):
        self.additional_properties.update(properties_dict)

    def get_properties(self,
                       requested_properties=DEFAULT_LINE_PROPERTIES,
                       additional_properties=None):

        # create dictionary to hold record properties
        properties = {}

        # query requested properties for record
        for property_name in requested_properties:
            if property_name == 'node':
                properties['node'] = self.node
            elif property_name == 'bearing':
                properties['bearing'] = self.bearing
            elif property_name == 'lat_start':
                properties['lat_start'] = self.start_position().latitude()
            elif property_name == 'lon_start':
                properties['lon_start'] = self.start_position().longitude()
            elif property_name == 'elev_start':
                properties['elev_start'] = self.start_position().elevation()
            elif property_name == 'lat_end':
                properties['lat_end'] = self.end_position().latitude()
            elif property_name == 'lon_end':
                properties['lon_end'] = self.end_position().longitude()
            elif property_name == 'elev_end':
                properties['elev_end'] = self.end_position().elevation()
            elif property_name == 'daynight':
                properties['daynight'] = self.day_night()

        # include properties specified at object creation
        if not self.additional_properties is None:
            properties.update(self.additional_properties)

        # include properties specified at record request
        if not additional_properties is None:
            properties.update(additional_properties)

        return properties

    # generate record for storage in shapefile or geojson
    def line_record(
            self,
            requested_properties=DEFAULT_LINE_PROPERTIES,
            additional_properties=None):
        """
        Generates record for storage in shapefile.
        :param requested_properties:
            list of properties to include in record as string
            node: ascending or descening node
            bearing: compass bearing from start point to end point
            in degrees clockwise from north
            lat_start: starting latitude in degrees
            lon_start: starting longitude in degrees
            elev_start: starting elevation in meters
            lat_end: ending latitude in degrees
            lon_end: ending longitude in degrees
            elev_end: ending elevation in meters
        :param additional_properties:
            any additional attributes to append to the record properties
        :return:
            python dictionary in spatial format containing geometry and properties
        """

        # create record from mapped geometry and properties
        record = {
            'geometry': mapping(self.line()),
            'properties': self.get_properties(requested_properties, additional_properties)
        }

        return record

    # generate record for storage in shapefile or geojson
    def polygon_record(
            self,
            sensor,
            requested_properties=DEFAULT_LINE_PROPERTIES,
            additional_properties=None):
        """
        Generates record for storage in shapefile.
        :param requested_properties:
            list of properties to include in record as string
            node: ascending or descening node
            bearing: compass bearing from start point to end point
            in degrees clockwise from north
            lat_start: starting latitude in degrees
            lon_start: starting longitude in degrees
            elev_start: starting elevation in meters
            lat_end: ending latitude in degrees
            lon_end: ending longitude in degrees
            elev_end: ending elevation in meters
        :param additional_properties:
            any additional attributes to append to the record properties
        :return:
            python dictionary in spatial format containing geometry and properties
        """

        # create record from mapped geometry and properties
        record = {
            'geometry': mapping(self.polygon(sensor)),
            'properties': self.get_properties(requested_properties, additional_properties)
        }

        return record


# class to represent a satellite sensor
class Sensor:
    def __init__(
            self,
            name,
            satellite,
            cross_track_swath_width_meters,
            swath_duration_minutes):
        self.name = name
        self.satellite = satellite
        self._cross_track_swath_width = cross_track_swath_width_meters
        self.swath_duration_minutes = swath_duration_minutes

    def __str__(self):
        return self.name

    def __repr__(self):
        return """<Sensor
        satellite: {}
        cross_track_swath_width: {}
        swath_duration_minutes: {}
        >""".format(
            self.satellite,
            self._cross_track_swath_width,
            self.swath_duration_minutes
        )

    def cross_track_swath_width(self, satellite_position):
        # this is the place to take elevation into account in calculating cross-track swath width
        return self._cross_track_swath_width


# class to contain a set of swath ground tracks from a satellite sensor
class OrbitHistory:
    def __init__(self, sensor, swath_ground_track_list=None, lock=None):
        if lock is None:
            lock = DummyLock()

        self.lock = lock
        self.swath_ground_tracks = []
        self.sensor = sensor

        if not swath_ground_track_list is None:
            for swath_ground_track in swath_ground_track_list:
                self.add_swath(swath_ground_track)

    def __str__(self):
        return "<OrbitHistory\n%s\n>" % (
            "\n".join([
                "\t%s" % str(swath)
                for swath
                in self.swath_ground_tracks
            ])
        )

    def add_swath(self, swath_ground_track):
        self.swath_ground_tracks.append(swath_ground_track)

    def remove_swath(self, swath_ground_track):
        self.swath_ground_tracks.remove(swath_ground_track)

    def swaths_intersecting_polygon(self, target_polygon_latlon, use_split_geometry=False):
        swaths = []

        for swath in self.swath_ground_tracks:
            swath_polygon = swath.polygon(self.sensor, use_split_geometry=use_split_geometry)

            

            if target_polygon_latlon.intersects(swath_polygon):
                swaths.append(swath)

        return swaths

    def write_line_shapefile(self, filename, shapefile_properties=DEFAULT_SHAPEFILE_PROPERTIES):
        line_schema = {
            'geometry': 'LineString',
            'properties': shapefile_properties
        }

        line_records = [
            swath.line_record()
            for swath
            in self.swath_ground_tracks
        ]

        # write records to shapefile
        with fiona.open(filename, 'w', driver=SHAPEFILE_DRIVER_NAME, crs=WGS84_FIONA_CRS, schema=line_schema) as f:
            for line_record in line_records:
                f.write(line_record)

    def write_polygon_shapefile(self, filename, shapefile_properties=DEFAULT_SHAPEFILE_PROPERTIES):
        polygon_schema = {
            'geometry': 'Polygon',
            'properties': shapefile_properties
        }

        polygon_records = [
            swath.polygon_record(self.sensor)
            for swath
            in self.swath_ground_tracks
        ]

        # write records to shapefile
        with fiona.open(filename, 'w', driver=SHAPEFILE_DRIVER_NAME, crs=WGS84_FIONA_CRS,
                        schema=polygon_schema) as f:
            for polygon_record in polygon_records:
                f.write(polygon_record)

    def day_scenes(self):
        return OrbitHistory(
            self.sensor,
            [
                swath
                for swath
                in self.swath_ground_tracks
                if swath.day_night() == 'day'
            ]
        )

    def intersect_polygon(self, target_polygon, use_split_geometry=False):
        return OrbitHistory(self.sensor, self.swaths_intersecting_polygon(
            target_polygon,
            use_split_geometry=use_split_geometry
        ))
