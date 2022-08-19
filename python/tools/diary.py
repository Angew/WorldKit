import argparse
import enum
import os
import os.path
import math
import skyfield.almanac
import skyfield.api
import svgwrite
import sys


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


class Season(enum.IntEnum):
    """
    My enum to match values skyfield uses.
    """
    Spring = 0
    Summer = 1
    Autumn = 2
    Winter = 3


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
        return self._units(val)

    def from_tl(self, x, y):
        return (self._units(self._top_left[0]+x), self._units(self._top_left[1]+y))

    def from_tr(self, x, y):
        return (self._units(self._bottom_right[0]-x), self._units(self._top_left[1]+y))

    def from_bl(self, x, y):
        return (self._units(self._top_left[0]+x), self._units(self._bottom_right[1]-y))

    def from_br(self, x, y):
        return (self._units(self._bottom_right[0]-x), self._units(self._bottom_right[1]-y))

    def top_left(self):
        return map(self._units, self._top_left)

    def bottom_right(self):
        return map(self._units, self._bottom_right)

    def size(self):
        return map(self._units, self._size)

    def subpage(self, top_left, size):
        return self.__class__(tuple(p+q for p, q in zip(self._top_left, top_left)), size)


class Diary:
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
        dwg = svgwrite.Drawing("diary.svg", (u.mm(297), u.mm(210)))
        page_size = (100+self.ROW_HEIGHT, 7*self.ROW_HEIGHT)
        self.add_page(dwg, DiaryPage((15, 5), page_size))
        self.add_page(dwg, DiaryPage((297/2+15, 5), page_size))
        dwg.saveas(os.path.join(self.options.output, dwg.filename))

    def add_page(self, dwg, page):
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

    def draw_sun(self, dwg, area):
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


def play(options):
    """
    Development command to test/run whatever I need in the script.

    This function exists solely for quick Python code prototyping and can thus change wildly.
    """
    print(timescale.tt_jd(2110763).utc_jpl())


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

    command_play = subcommands.add_parser(
        "play",
        description="Development command to test/run something in the script",
    )
    command_play.set_defaults(process=play)

    options = parser.parse_args(args[1:])
    options.process(options)


if __name__ == "__main__":
    main(sys.argv)
