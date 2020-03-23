import argparse
import enum
import os.path
import math
import skyfield.almanac
import skyfield.api
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


# Julian date = days since 1 Jan -4713 noon
# Time in Julian date is UT (=UT1), which is successor to GMT, so it's centred on Greenwhich.

script_dir = sys.path[0]

timescale = skyfield.api.load.timescale(builtin=True)
sf_load = skyfield.api.Loader(os.path.join(script_dir, "skyfield_data", "loaded"))
planets = sf_load("de406.bsp")  # goes from -3000 to +3000


class Season(enum.IntEnum):
    Spring = 0
    Summer = 1
    Autumn = 2
    Winter = 3


def find_epoch(options):
    print("Started")
    seasons = skyfield.almanac.seasons(planets)
    def isWinter(time):
        return seasons(time) == Season.Winter
    def isSafe(time):
        fract = math.modf(time.tt + 0.5)[0]
        return epsilon < fract < 1-epsilon
    isWinter.rough_period = 365.0
    epsilon = 1/24/2 # Half-hour precision is good enough
    midpoint = 2_100_000 # Roughly +1000
    period = 500 * 365
    t, s = skyfield.almanac.find_discrete(
        # timescale.tt_jd(1_000_000)  # Roughly -2000
        # , timescale.tt_jd(2_500_000)  # Roughly +2000
        timescale.tt_jd(midpoint - period/2)
        , timescale.tt_jd(midpoint + period/2)
        , isWinter
        , epsilon = epsilon
    )
    solstices = (p[0] for p in zip(t, s) if p[1])
    candidates = (p for p in lookahead(solstices) if
        isSafe(p[0]) and isSafe(p[1])
        and math.floor(p[1].tt + 0.5) - math.floor(p[0].tt + 0.5) > 365
    )
    # last = next(solstices)
    # candidates = []
    # for s in solstices:
        # if math.floor(s.tt + 0.5) - math.floor(last.tt + 0.5) > 365:
            # candidates.append((last, s))
        # last = s
    # candidates = filter(lambda c: isSafe(c[0]) and isSafe(c[1]), candidates)
    print("\n".join(f"{p.tt} ({p.utc_jpl()}) - ({s.utc_jpl()})" for p, s in candidates))


def play(options):
    for i in lookahead([1, 2, 3, 4, 5, 6]):
        print(i)


def main(args):
    parser = argparse.ArgumentParser(
        prog=os.path.basename(args[0]),
        description="Diary processing tool",
    )
    parser.set_defaults(process=lambda _: parser.parse_args(['-h']))
    subcommands = parser.add_subparsers()

    command_find_epoch = subcommands.add_parser(
        "find-epoch"
    )
    command_find_epoch.set_defaults(process=find_epoch)

    command_play = subcommands.add_parser(
        "play"
    )
    command_play.set_defaults(process=play)

    options = parser.parse_args(args[1:])
    options.process(options)


if __name__ == "__main__":
    main(sys.argv)
