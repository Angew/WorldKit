import argparse
import os.path
import skyfield.almanac
import skyfield.api
import sys


# Julian date = days since 1 Jan -4713

timescale = skyfield.api.load.timescale(builtin=True)
sf_load = skyfield.api.Loader(os.path.join(os.path.dirname(sys.path[0]), "skyfield_data", "loaded"))
planets = sf_load("de406.bsp")  # goes from -3000 to +3000


def find_epoch(options):
    print("Started")
    t, s = skyfield.almanac.find_discrete(
        # timescale.tt_jd(1_000_000)  # Roughly -2000
        # , timescale.tt_jd(2_500_000)  # Roughly +2000
        timescale.tt_jd(2_000_000)
        , timescale.tt_jd(2_001_000)
        , skyfield.almanac.seasons(planets)
        , epsilon = 1/24/2  # Half-hour precision is good enough
    )
    solstices = (p for p in zip(t, s) if p[1] == 3)
    print("\n".join((str(p[0].tt) for p in solstices)))


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

    options = parser.parse_args(args[1:])
    options.process(options)


if __name__ == "__main__":
    main(sys.argv)
