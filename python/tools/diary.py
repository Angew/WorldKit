import argparse
import enum
import os
import os.path
import math
import skyfield.almanac
import skyfield.api
from skyfield.api import wgs84
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
# enough for me.


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

class Landmarks:
    Brasel = wgs84.latlon(55.6 * skyfield.api.N, 0)

    north_tip = wgs84.latlon(60.165 * skyfield.api.N, 1.25 * skyfield.api.E)

    north_of_Brasel = wgs84.latlon(north_tip.latitude.degrees, Brasel.longitude.degrees)

    @staticmethod
    def get_midpoint(a, b):
        return wgs84.latlon(
            (a.latitude.degrees + b.latitude.degrees)/2,
            (a.longitude.degrees + b.longitude.degrees)/2
        )


EPOCH_MIDNIGHT = 2110762.5 # Midnight on start of Solstice day


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
            epsilon=1/24/60 # 1 minute precision
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
            epsilon=1/24/60 # 1 minute precision
        )):
            if up:
                return t
        return None


    @classmethod
    def moonrises(cls, point, when=None):
        when = when or range(366)
        return [cls.moonrise(point, day) for day in when]


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
                epsilon=1/24/60 # 1 minute precision
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
                epsilon=1/24/60 # 1 minute precision
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
                epsilon=1/24/60 # 1 minute precision
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


def preprocess(options):
    if options.tide:
        preprocess_tide(options)


def preprocess_tide(options):
    path, in_file_name = os.path.split(options.tide)
    if not in_file_name:
        raise ValueError("Tide file database specified incorrectly")
    out_file_name = os.path.splitext(in_file_name)[0] + ".data.py"
    # * accumulate maxima in streaming fashion
    # * write them out


def play(options):
    """
    Development command to test/run whatever I need in the script.

    This function exists solely for quick Python code prototyping and can thus change wildly.
    """
    print(planets)


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

    command_play = subcommands.add_parser(
        "play",
        description="Development command to test/run something in the script",
    )
    command_play.set_defaults(process=play)

    options = parser.parse_args(args[1:])
    options.process(options)


if __name__ == "__main__":
    main(sys.argv)
