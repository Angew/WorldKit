import argparse
from collections import namedtuple
import datetime
import enum
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import os.path
import skyfield.almanac
import skyfield.api
from skyfield.api import utc, wgs84
import skyfield.searchlib
import svgwrite
import sys


# ToDo and notes
#
# * Invent tide formula
#   o Validate against some db
#   (Tides by https://moon.nasa.gov/moon-in-motion/tides/)
#
# * Rearrange left box

# Notes on Tides
#
# I won't get the landmass-dependent parts right anyway, so super-fine
# accuracy is unnecessary. Let's set high tide to "50 minutes before
# Earth's centre - Observer - Moon form a line in the equatorial plane."
# Note that there may not be a need for actual projection: longitude and
# right ascension are both along the equator, and I don't care about
# distances, only about directions. However, the vernal equinox (right
# ascension 0h) is fixed in space, it doesn't rotate with the Earth.
# So, IIUC, "form a line" means that the Observer's RA is the same as
# the Moon's. And since it's all wibbly-wobbly anyway, that sounds good
# enough for me. Note that Skyfield has some helpers like "which geo is below
# a given RA", which might be just what I need.


def adjacent_filter(it, pred):
    it = iter(it)
    try:
        last = next(it)
        for current in it:
            if pred(last, current):
                yield last
            last = current
    except StopIteration:
        pass


def lookahead(it, scope=1):
    """
    `scope+1`-sized sliding window iterating over `it`
    """
    it = iter(it)
    view = []
    try:
        for i in range(scope):
            view.append(next(it))
        for item in it:
            view.append(item)
            yield tuple(view)
            del view[0]
    except StopIteration:
        pass


# tt = Terrestrial Time (as a Julian date)
# Julian date = days since 1 Jan -4713 noon
# Time in Julian date is UT (=UT1), which is successor to GMT, so it's centred on Greenwhich.

# Note that since tt is from noon, tt+0.5 is distance from starting midnight

script_dir = sys.path[0]

timescale = skyfield.api.load.timescale(builtin=True)
sf_load = skyfield.api.Loader(os.path.join(script_dir, "skyfield_data", "loaded"))
planets = sf_load("de406.bsp")  # goes from -3000 to +3000



# Data

EPOCH_MIDNIGHT = 2110762.5 # Midnight on start of Solstice day

MILES_PER_LAT = 69


class Landmarks:
    Brasel = wgs84.latlon(55.6 * skyfield.api.N, 0)

    north_tip = wgs84.latlon(60.165 * skyfield.api.N, 1.25 * skyfield.api.E)

    north_of_Brasel = wgs84.latlon(north_tip.latitude.degrees, Brasel.longitude.degrees)

    reference_Portsmouth = wgs84.latlon(50.80256 * skyfield.api.N, 1.11175 * skyfield.api.W)
    reference_antiPortsmouth = wgs84.latlon(50.80256 * skyfield.api.S, 181.11175 * skyfield.api.W)

    @staticmethod
    def get_midpoint(a, b):
        return wgs84.latlon(
            (a.latitude.degrees + b.latitude.degrees)/2,
            (a.longitude.degrees + b.longitude.degrees)/2
        )

    @staticmethod
    def lat_distance(a, b):
        return abs(a.latitude.degrees - b.latitude.degrees) * MILES_PER_LAT


def epsilon_minutes(m):
    return m/24/60


class Season(enum.IntEnum):
    """
    My enum to match values skyfield uses.
    """
    Spring = 0
    Summer = 1
    Autumn = 2
    Winter = 3


class Twilight(enum.IntEnum):
    """
    My enum to match return values of skyfield.almanac.dark_twilight_day
    """
    Night = 0
    Astronomical = 1
    Nautical = 2
    Civil = 3
    Day = 4


# Computations

class Compute:
    @classmethod
    def day_length(cls, point, day):
        goal_func = skyfield.almanac.sunrise_sunset(planets, point)
        t = skyfield.almanac.find_discrete(
            timescale.tt_jd(EPOCH_MIDNIGHT+day),
            timescale.tt_jd(EPOCH_MIDNIGHT+day+1),
            goal_func,
            epsilon=epsilon_minutes(1)
        )[0]
        return t[1] - t[0]

    @classmethod
    def day_lengths(cls, point, when=None):
        when = when or range(366)
        return [cls.day_length(point, day) for day in when]

    @classmethod
    def moonrise(cls, point, day):
        moon = planets["Moon"]
        goal_func = skyfield.almanac.risings_and_settings(
            planets, moon, point, radius_degrees=0.25)
        for t, up in zip(*skyfield.almanac.find_discrete(
            timescale.tt_jd(EPOCH_MIDNIGHT+day),
            timescale.tt_jd(EPOCH_MIDNIGHT+day+1),
            goal_func,
            epsilon=epsilon_minutes(1)
        )):
            if up:
                return t
        return None

    @classmethod
    def moonrises(cls, point, when=None):
        when = when or range(366)
        return [cls.moonrise(point, day) for day in when]

    @classmethod
    def north_dawn_distance(cls, day_start):
        t0 = day_start
        t1 = day_start + datetime.timedelta(days=1)
        north = Landmarks.north_of_Brasel
        south = Landmarks.Brasel
        while True:
            point = Landmarks.get_midpoint(north, south)
            dist = Landmarks.lat_distance(north, south)
            if dist < 1:
                return round(Landmarks.lat_distance(Landmarks.Brasel, point))
            goal_func = skyfield.almanac.dark_twilight_day(planets, point)
            times, lights = skyfield.almanac.find_discrete(
                t0,
                t1,
                goal_func,
                epsilon=epsilon_minutes(1)
            )
            if lights[0] > Twilight.Nautical:
                north = point
            else:
                south = point


# Utilities

MeasuredTides = namedtuple("MeasuredTides", "high_tides, low_tides")

def load_tides():
    tide_data = {}
    with open("data/tide.data.py") as tide_file:
        exec(tide_file.read(), {}, tide_data)
    return MeasuredTides(tide_data["measured_high_tides"], tide_data["measured_low_tides"])


# Commands

def find_epoch(options):
    """
    List candidates suitable as epoch of the game.
    """
    print("Started")
    seasons = skyfield.almanac.seasons(planets)
    def isWinter(time):
        return seasons(time) == Season.Winter
    def isSafe(time):
        """Is the time point safely within the day wrt. epsilon precision?"""
        fract = math.modf(time.tt + 0.5)[0]
        return epsilon < fract < 1-epsilon
    isWinter.rough_period = 365.0
    epsilon = 1/24/2 # Half-hour precision is good enough
    midpoint = 2_100_000 # Roughly +1000 AD
    period = 500 * 365 # Search 500 year interval
    event_time, event_is_winter = skyfield.almanac.find_discrete(
        timescale.tt_jd(midpoint - period/2)
        , timescale.tt_jd(midpoint + period/2)
        , isWinter
        , epsilon = epsilon
    )
    solstices = (p[0] for p in zip(event_time, event_is_winter) if p[1])
    candidates = (p for p in lookahead(solstices) if
        isSafe(p[0]) and isSafe(p[1])
        and math.floor(p[1].tt + 0.5) - math.floor(p[0].tt + 0.5) > 365 # leap year
    )
    print("Julian starting winter solstice (dates for verification)")
    print("\n".join(f"{p.tt} ({p.utc_jpl()}) - ({s.utc_jpl()})" for p, s in candidates))


class SvgUnits:
    def cm(self, val):
        return f"{val}cm"

    def mm(self, val):
        return f"{val}mm"


class DiaryPage:
    """
    Offset manager for one diary page (the drawn rectangle)
    """
    def __init__(self, top_left, size, units=SvgUnits.mm):
        self._size = size
        self._top_left = top_left
        self._bottom_right = tuple(p+s for p, s in zip(top_left, size))
        self._units = lambda v: units(SvgUnits, v)

    def u(self, val):
        """Convert `val` using page's units"""
        return self._units(val)

    def from_tl(self, x, y):
        """Convert units offset from page's top-left corner. Positive offset into page"""
        return (self._units(self._top_left[0]+x), self._units(self._top_left[1]+y))

    def from_tr(self, x, y):
        """Convert units offset from page's top-right corner. Positive offset into page"""
        return (self._units(self._bottom_right[0]-x), self._units(self._top_left[1]+y))

    def from_bl(self, x, y):
        """Convert units offset from page's bottom-left corner. Positive offset into page"""
        return (self._units(self._top_left[0]+x), self._units(self._bottom_right[1]-y))

    def from_br(self, x, y):
        """Convert units offset from page's bottom-right corner. Positive offset into page"""
        return (self._units(self._bottom_right[0]-x), self._units(self._bottom_right[1]-y))

    def top_left(self):
        """Page's top-left corner converted in units"""
        return map(self._units, self._top_left)

    def bottom_right(self):
        """Page's bottom-right corner converted in units"""
        return map(self._units, self._bottom_right)

    def size(self):
        """Page's size converted in units"""
        return map(self._units, self._size)

    def subpage(self, top_left, size):
        """Return a new page offset from page's top left corner"""
        return self.__class__(tuple(p+q for p, q in zip(self._top_left, top_left)), size)


class Diary:
    # Information in one row:
    # Party notes
    # World notes
    # Weather
    # Varkania date
    # Earth date
    # Moon phase
    # Dawn time: Brasel + northtip OR Brasel + miles to north dawn point
    # Sunrise time: Brasel + northtip
    # Sunset time: Brasel + northtip
    # Dusk time: Brasel + northtip OR Brasel + miles to north dusk point
    # Tide times
    # ...
    # Lunar eclipse, if any
    # Solar eclipse, if any
    def __init__(self, options):
        self.options = options
        self.ROW_HEIGHT = 28.5
        if os.pardir in options.output:
            raise ValueError(f"{os.pardir} in output path is not supported yet.")
        os.makedirs(options.output, exist_ok=True)
        self.fill_stroke = {
            "fill": "none",
            "stroke": "black"
        }

    def generate(self):
        self.compute_data()
        u = SvgUnits()
        # ToDo: svg pages
        # ToDo: enough pages for entire year
        dwg = svgwrite.Drawing("generated-diary.svg", (u.mm(297), u.mm(210)))
        page_size = (100+self.ROW_HEIGHT, 7*self.ROW_HEIGHT)
        self.add_page(dwg, DiaryPage((15, 5), page_size))
        self.add_page(dwg, DiaryPage((297/2+15, 5), page_size))
        dwg.saveas(os.path.join(self.options.output, dwg.filename))

    def add_page(self, dwg, page):
        """
        Draw one page into `dwg` (sized as `page`)
        """
        # Notes to change: 7mm x 3mm is enough for raw time text
        u = SvgUnits()
        # Boundary
        dwg.add(dwg.rect(
            insert=page.top_left(),
            size=page.size(),
            **self.fill_stroke
        ))
        # Splitters
        for row in range(1, 7):
            y = row*self.ROW_HEIGHT
            dwg.add(dwg.line(
                start=page.from_tl(0, y),
                end=page.from_tr(0, y),
                stroke_width=2,
                **self.fill_stroke
            ))
        for x in self.ROW_HEIGHT, self.ROW_HEIGHT+60:
            dwg.add(dwg.line(
                start=page.from_tl(x, 0),
                end=page.from_bl(x, 0),
                **self.fill_stroke
            ))
        # Rows
        for row_num in range(7):
            row = page.subpage((0, row_num * self.ROW_HEIGHT), (page._size[0], self.ROW_HEIGHT))
            square = row.subpage((0, 0), (self.ROW_HEIGHT, self.ROW_HEIGHT))
            # Lines
            y = 0
            for offset in 7, 7, 5:
                y += offset
                dwg.add(dwg.line(
                    start=square.from_bl(0, y),
                    end=square.from_br(0, y),
                    **self.fill_stroke
                ))
            for x in 7, 7 + (self.ROW_HEIGHT-7)/2:
                dwg.add(dwg.line(
                    start=square.from_bl(x, 14),
                    end=square.from_bl(x, 0),
                    **self.fill_stroke
                ))
            dwg.add(dwg.line(
                start=row.from_br(40, 10),
                end=row.from_br(0, 10),
                **self.fill_stroke
            ))
            # Content
            self.draw_sun(dwg, square.subpage((0, self.ROW_HEIGHT-14), (7, 7)))
            # ToDo: actual varying content here
            dwg.add(dwg.text(
                "23:59",
                insert=square.from_bl(7.5, 7.5),
                font_family="sans-serif",
                font_size=10,
            ))

    def draw_sun(self, dwg, area):
        """
        Draw sun icon in the middle of `area`

        `area` must be square
        """
        c = area._size[0]/2
        dwg.add(dwg.circle(
            center=(area.from_tl(c, c)),
            r=area.u(c*0.3),
            **self.fill_stroke
        ))
        r1 = c*0.4
        r2 = c*0.8
        for i in range(12):
            angle = i * math.pi/6
            dwg.add(dwg.line(
                start=area.from_tl(c + math.cos(angle)*r1, c + math.sin(angle)*r1),
                end=area.from_tl(c + math.cos(angle)*r2, c + math.sin(angle)*r2),
                **self.fill_stroke
            ))


def generate_diary(options):
    """
    Generate all necessary SVG file(s) for one year of a diary
    """
    Diary(options).generate()


def analyse(options):
    """
    Display various analyses (e.g. twilight times)

    Used to guide DM's decisions on which things must be present in the diary
    and which can be interpolated/eyeballed.
    """
    # Is there a day when the northmost point of Varkania experiences no full night?
    # Known results: days 142-224
    if options.twilight:
        print("Twilight")
        print("========")
        def goal_func(time):
            return skyfield.almanac.dark_twilight_day(planets, Landmarks.north_tip)(time) < Twilight.Nautical
        goal_func.rough_period = 0.5

        for day in range(366):
            times, events = skyfield.almanac.find_discrete(
                timescale.tt_jd(EPOCH_MIDNIGHT+day),
                timescale.tt_jd(EPOCH_MIDNIGHT+day+1),
                goal_func,
                epsilon=epsilon_minutes(1)
            )

            if not times:
                print("No full night on day", day)

    # Is there a day when Brasel experiences no full night?
    # Known results: days 165-202
    if options.twilight_brasel:
        print("Twilight in Brasel")
        print("==================")
        def goal_func(time):
            return skyfield.almanac.dark_twilight_day(planets, Landmarks.Brasel)(time) < Twilight.Nautical
        goal_func.rough_period = 0.5

        for day in range(366):
            times, events = skyfield.almanac.find_discrete(
                timescale.tt_jd(EPOCH_MIDNIGHT+day),
                timescale.tt_jd(EPOCH_MIDNIGHT+day+1),
                goal_func,
                epsilon=epsilon_minutes(1)
            )

            if not times:
                print("No full night on day", day)

    # Is there a day when the northmost point of Varkania experiences nothing but daylight?
    # Known results: no
    if options.dawn:
        print("Dawn")
        print("====")
        def goal_func(time):
            return skyfield.almanac.dark_twilight_day(planets, Landmarks.north_tip)(time) < Twilight.Day
        goal_func.rough_period = 0.5

        for day in range(366):
            times, events = skyfield.almanac.find_discrete(
                timescale.tt_jd(EPOCH_MIDNIGHT+day),
                timescale.tt_jd(EPOCH_MIDNIGHT+day+1),
                goal_func,
                epsilon=epsilon_minutes(1)
            )

            if not times:
                print("No twilight on day", day)

    # Variation in day lengths between Brasel and Varkania's northmost point
    if options.day_variation:
        print("Day length variation")
        print("====================")
        days = {}
        for point in "north_tip", "Brasel":
            days[point] = Compute.day_lengths(getattr(Landmarks, point))
        variations = [abs(n - b)*24 for n, b in zip(*days.values())]
        for v in variations:
            print("Variation", v)
        m, M = min(variations), max(variations)
        print("Largest variation [h]:", M)
        print("Smallest variation [h]:", m)

    # Difference between midpoint day length interpolation and computation
    # Known data: difference is at most 0.088 h => linearity's fine
    if options.day_linearity:
        print("Day length linearity")
        print("====================")
        days = {}
        for point in "north_tip", "Brasel":
            days[point] = Compute.day_lengths(getattr(Landmarks, point))
        computed = Compute.day_lengths(Landmarks.get_midpoint(Landmarks.Brasel, Landmarks.north_tip))
        linear = [(a+b)/2 for a, b in zip(*days.values())]
        diffs = [abs(c-l)*24 for c, l in zip(computed, linear)]
        m, M = min(diffs), max(diffs)
        print("Largest diff [h]:", M)
        print("Smallest diff [h]:", m)

    # Variation in Moonrise time between Brasel and northtip
    # Known result: 29min
    if options.moonrise_variation:
        print("Moonrise variation")
        print("==================")
        rises = {}
        for point in "north_of_Brasel", "Brasel":
            rises[point] = Compute.moonrises(getattr(Landmarks, point))

        variations = [(n - b)*24*60 for n, b in zip(*rises.values()) if n is not None and b is not None]
        for v in variations:
            print("Variation", v)
        variations = list(map(abs, variations))
        m, M = min(variations), max(variations)
        print("Largest abs variation [min]:", M)
        print("Smallest abs variation [min]:", m)

    # Difference between moonrise interpolation and computation
    # Known data: difference is at most 1.3 min => linearity's fine
    if options.moonrise_linearity:
        print("Moonrise linearity")
        print("==================")
        rises = {}
        for point in "north_of_Brasel", "Brasel":
            rises[point] = Compute.moonrises(getattr(Landmarks, point))
        computed = Compute.moonrises(Landmarks.get_midpoint(Landmarks.Brasel, Landmarks.north_of_Brasel))
        def avg(a, b):
            if a is None or b is None:
                return None
            return timescale.tt_jd((a.tt+b.tt)/2)
        linear = [avg(a, b) for a, b in zip(*rises.values())]
        diffs = [abs(c-l)*24*60 for c, l in zip(computed, linear) if c is not None and l is not None]
        m, M = min(diffs), max(diffs)
        print("Largest diff [min]:", M)
        print("Smallest diff [min]:", m)

    # Difference between tides computed using my formula and those from official database
    if options.tide_accuracy:
        print("Tide accuracy")
        print("=============")

        measured = load_tides()
        measured_high_tides = measured.high_tides

        year = 2020 # FuWo: compute from tide data
        # Tide computation algorithm (idea):
        # * As MVP, ignore 50' difference
        # * High tide happens when observer longitude is beneath Moon's RA
        # * For a find_discrete, I can use sign of "observer.lon - Moon.lon_below" (mod. 180 somehow)
        #   * Better, there's find_minima, which should work even better for me
        # Let's try that
        earth = planets["Earth"]
        moon = planets["Moon"]
        def distance_from_high_tide(time):
            m = earth.at(time).observe(moon)
            tide_lon = wgs84.latlon_of(m)[1]
            dist = np.vstack([
                abs(tide_lon.degrees - Landmarks.reference_Portsmouth.longitude.degrees),
                abs(tide_lon.degrees - Landmarks.reference_antiPortsmouth.longitude.degrees)
            ])
            return dist.min(axis=0)
        distance_from_high_tide.rough_period = 0.5
        computed_high_tides = skyfield.searchlib.find_minima(
            timescale.utc(year),
            timescale.utc(year+1),
            distance_from_high_tide,
            epsilon=epsilon_minutes(1)
        )
        prefix = 100
        fig, ax = plt.subplots()
        measured = [t[0] for t in measured_high_tides[:prefix]]
        computed = [t.tt for t in computed_high_tides[0][:prefix]]
        difference = [(m - c)*24*60 for m, c in zip(measured, computed)]
        ax.plot(difference, "bs", label="Difference [min]")
        ax.legend()
        ax.grid(True)
        fig.set_figwidth(50)
        fig.savefig("analyse/high-tides-diffs.png")
        # * plot both measured and computed data
        # * write differences to csv file (I want them analysable manually)
        # * print summary: largest difference, smallest difference, mean difference


def preprocess(options):
    if options.tide:
        preprocess_tide(options)
    if options.diary:
        preprocess_diary(options)


def preprocess_tide(options):
    path, in_file_name = os.path.split(options.tide)
    if not in_file_name:
        raise ValueError("Tide file database specified incorrectly")
    out_file_name = "tide.data.py"

    high_tide_idxs = []
    low_tide_idxs = []
    tide_data = []
    MeasurementTuple = namedtuple("MeasurementTuple", "time, level")

    # Read raw data into MeasurementTuple list
    RawDataTuple = namedtuple("RawDataTuple", "number, date, time, level, residue_")
    def Measurement(data):
        date_string = f'{data.date.replace("/", "-")}.{data.time}+00:00'
        return MeasurementTuple(timescale.from_datetime(datetime.datetime.fromisoformat(date_string)).tt, float(data.level.strip("T")))
    last_measurement = None
    rising = True
    with open(options.tide) as in_file:
        for line in filter(None, in_file):
            try:
                raw_data = RawDataTuple(*line.split())
            except TypeError:
                continue
            if not raw_data.number.endswith(")"):
                continue
            measurement = Measurement(raw_data)
            tide_data.append(measurement)
            if not last_measurement:
                last_measurement = measurement
                continue
            if measurement.level < last_measurement.level:
                if rising:
                    high_tide_idxs.append(len(tide_data) - 2)
                    rising = False
            elif measurement.level > last_measurement.level:
                if not rising:
                    low_tide_idxs.append(len(tide_data) - 2)
                    rising = True
            last_measurement = measurement

    # Remove false extrema
    def remove_false_extrema(extrema_idxs, left_is_extremer, window_reach=8):
        try:
            idx_idx_ex = 0
            # No for loop, since extrema_idxs is modified inside
            # Reaching beyon the end of extrema_idxs is used as stop condition
            while True:
                idx_ex = extrema_idxs[idx_idx_ex]
                is_extremum = True
                ex_level = tide_data[idx_ex].level
                test_range = (
                    list(range(max(idx_ex - window_reach, 0), idx_ex)) +
                    list(range(idx_ex + 1, min(idx_ex + window_reach + 1, len(tide_data))))
                )
                ex_value_idxs = [idx_ex]
                for idx_test in test_range:
                    test_level = tide_data[idx_test].level
                    if left_is_extremer(test_level, ex_level):
                        is_extremum = False
                        break
                    elif test_level == ex_level:
                        ex_value_idxs.append(idx_test)
                if is_extremum:
                    if ex_value_idxs[-1] - ex_value_idxs[0] > len(ex_value_idxs)-1:
                        raise ValueError(f"Ambiguous extremum at data nrs: {(i+1 for i in ex_value_idxs)}")
                    next_valid = idx_ex + window_reach + 1
                    idx_idx_ex += 1
                    while extrema_idxs[idx_idx_ex] < next_valid:
                        del extrema_idxs[idx_idx_ex]
                else:
                    del extrema_idxs[idx_idx_ex]
        except IndexError:
            pass
    remove_false_extrema(high_tide_idxs, lambda l, r: l > r)
    remove_false_extrema(low_tide_idxs, lambda l, r: l < r)

    # Write out preprocessed data
    with open(os.path.join(path, out_file_name), 'w') as out_file:
        for tides in "high_tide", "low_tide":
            out_file.write(f"measured_{tides}s = [\n")
            data = (tide_data[i] for i in locals()[tides+"_idxs"])
            out_file.write(",\n".join((f"    ({x.time}, {x.level})" for x in data)))
            out_file.write("\n]\n")


class DayData:
    """
    All relevant diary data for one day
    """
    def __init__(self, jd_start):
        self.t_start = timescale.tt_jd(jd_start)
        self.t_end = timescale.tt_jd(jd_start+1)
        for tip in "north", "south":
            for event in "dawn", "sunrise", "sunset", "dusk":
                setattr(self, f"{tip}_{event}", None)
        self.north_dawn_distance = None


def compute_diary_data(day_data):
    # Regular Sun
    points = {
        "south": Landmarks.Brasel,
        "north": Landmarks.north_of_Brasel
    }
    for tip, point in points.items():
        goal_func = skyfield.almanac.dark_twilight_day(planets, point)
        times, events = skyfield.almanac.find_discrete(
            day_data[0].t_start,
            day_data[-1].t_end,
            goal_func,
            epsilon=epsilon_minutes(1)
        )
        last_light = goal_func(day_data[0].t_start).item()
        for time, light in zip(times, events):
            day = day_data[int(time.tt - EPOCH_MIDNIGHT)]
            if last_light < light:
                if light == Twilight.Nautical:
                    setattr(day, f"{tip}_dawn", time)
                elif light == Twilight.Day:
                    setattr(day, f"{tip}_sunrise", time)
            else:
                if light == Twilight.Civil:
                    setattr(day, f"{tip}_sunset", time)
                elif light == Twilight.Astronomical:
                    setattr(day, f"{tip}_dusk", time)
            last_light = light
    # Northmost Sun
    for day in day_data:
        if day.north_dawn is None:
            day.north_dawn_distance = Compute.north_dawn_distance(day.t_start)


def preprocess_diary(options):
    day_data = [DayData(EPOCH_MIDNIGHT+day) for day in range(366)]
    compute_diary_data(day_data)
    # Debug
    with open(os.path.join(script_dir, "data", "diary.csv"), "w") as f:
        f.write("Day,Dawn,Sunrise,Sunset,Dusk\n")
        for day in day_data:
            f.write(day.t_start.utc_strftime("%d.%b,"))
            if day.north_dawn_distance is not None:
                f.write(f"{day.north_dawn_distance} miles,")
            else:
                f.write(day.north_dawn.utc_strftime("%H:%M,"))
            f.write(day.north_sunrise.utc_strftime("%H:%M,"))
            f.write(day.north_sunset.utc_strftime("%H:%M,"))
            if day.north_dusk is not None:
                f.write(day.north_dusk.utc_strftime("%H:%M"))
            else:
                f.write("--")
            f.write("\n")


def verify(options):
    out_dir = os.path.join(script_dir, "verify")
    os.makedirs(out_dir, exist_ok=True)

    # Computed sun event data
    if options.sun:
        t_start = timescale.tt_jd(EPOCH_MIDNIGHT)
        def write_sun_data(file_name, events, initial_value):
            with open(os.path.join(out_dir, file_name), "w") as file:
                file.write(f"Initial value:,{Twilight(initial_value).name}")
                last_day = -1
                for time, light in zip(*events):
                    day = time.tt_calendar()[2]
                    if day != last_day:
                        file.write(f"\n{time.tt_strftime('%d.%b')}")
                        last_day = day
                    file.write(f",{Twilight(light).name},{time.tt_strftime('%H:%M')}")
                file.write("\n")
        # North Sun data
        goal_func = skyfield.almanac.dark_twilight_day(planets, Landmarks.north_of_Brasel)
        coarse_events = skyfield.searchlib.find_discrete(
            t_start,
            timescale.tt_jd(t_start.tt+366),
            goal_func,
            epsilon=epsilon_minutes(1)
        )
        write_sun_data("north-sun.coarse.csv", coarse_events, goal_func(t_start))
        fine_events = skyfield.searchlib.find_discrete(
            t_start,
            timescale.tt_jd(t_start.tt+366),
            goal_func
        )
        write_sun_data("north-sun.fine.csv", fine_events, goal_func(t_start))

    # Measured tide data differences
    if options.measured_tide:
        data = [d[0] for d in load_tides().high_tides]
        diffs = [(t[1] - t[0])*24*60 for t in lookahead(data)]
        fig, ax = plt.subplots()
        #d-data[0] for d in data[1:]]
        ax.plot(diffs, label="Difference [min]")
        ax.legend()
        ax.grid(True)
        fig.set_figwidth(50)
        fig.savefig(os.path.join(out_dir, "measured-tide-diff.high.png"))


def play(options):
    """
    Development command to test/run whatever I need in the script.

    This function exists solely for quick Python code prototyping and can thus change wildly.
    """
    t = timescale.utc(2020)
    print(t)
    print(t.tt)


def main(args):
    parser = argparse.ArgumentParser(
        prog=os.path.basename(args[0]),
        description="Diary processing tool",
    )
    parser.set_defaults(process=lambda _: parser.parse_args(['-h']))
    subcommands = parser.add_subparsers()

    command_find_epoch = subcommands.add_parser(
        "find-epoch",
        description="Find a suitable epoch time based on leap year requirements",
    )
    command_find_epoch.set_defaults(process=find_epoch)

    command_generate_diary = subcommands.add_parser(
        "generate-diary",
        description="Generate SVG file(s) for a diary",
    )
    command_generate_diary.set_defaults(process=generate_diary)
    command_generate_diary.add_argument(
        "-o", "--output",
        help="Output directory",
        metavar="OUTPUT_DIR",
        default=".",
    )

    command_analyse = subcommands.add_parser(
        "analyse",
        description="Different analyses (e.g. twilight time)",
    )
    command_analyse.set_defaults(process=analyse)
    command_analyse.add_argument(
        "--twilight",
        help="Days when north tip has no full night",
        action="store_true",
    )
    command_analyse.add_argument(
        "--twilight-brasel",
        help="Days when Brasel has no full night",
        action="store_true",
    )
    command_analyse.add_argument(
        "--dawn",
        help="Days when north tip has full day",
        action="store_true",
    )
    command_analyse.add_argument(
        "--day-variation",
        help="Day length differences between north and Brasel",
        action="store_true",
    )
    command_analyse.add_argument(
        "--day-linearity",
        help="Day length differences from linearity at midpoint",
        action="store_true",
    )
    command_analyse.add_argument(
        "--moonrise-variation",
        help="Moonrise differences between north and Brasel",
        action="store_true",
    )
    command_analyse.add_argument(
        "--moonrise-linearity",
        help="Moonrise differences from linearity at midpoint",
        action="store_true",
    )
    command_analyse.add_argument(
        "--tide-accuracy",
        help="Difference between tide formula and official tide data",
        action="store_true",
    )

    command_compute = subcommands.add_parser(
        "preprocess",
        description="Pre-compute some data",
    )
    command_compute.set_defaults(process=preprocess)
    command_compute.add_argument(
        "--tide",
        help="Precompute tide from given database FILE",
        metavar="FILE",
    )
    command_compute.add_argument(
        "--diary",
        help="Precompute master diary data",
        action="store_true",
    )

    command_verify = subcommands.add_parser(
        "verify",
        description="Verify diary data",
    )
    command_verify.set_defaults(process=verify)
    command_verify.add_argument(
        "--sun",
        help="Dawn, sunrise, sunset, dusk times (or miles)",
        action="store_true",
    )
    command_verify.add_argument(
        "--measured-tide",
        help="Differences between measured tide points",
        action="store_true",
    )

    command_play = subcommands.add_parser(
        "play",
        description="Development command to test/run something in the script",
    )
    command_play.set_defaults(process=play)

    options = parser.parse_args(args[1:])
    options.process(options)


if __name__ == "__main__":
    main(sys.argv)
