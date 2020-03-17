import argparse
import os.path
import skyfield
import sys


# Julian date = days since 1 Jan -4713

timescale = skyfield.api.load.timescale(builtin=True)
planets = skyfield.api.load("de406.bsp")  # goes from -3000 to +3000


def find_epoch(options):
    skyfield.almanac.find_discrete(
        timescale.tt_jd(1_000_000)  # Roughly -2000
        , timescale.tt_jd(2_500_000)  # Roughly +2000
        , skyfield.almanac.seasons(planets)
        , epsilon = 1/24/2  # Half-hour precision is good enough
    )


def main(args):
    parser = argparse.ArgumentParser(
        prog=os.path.basename(args[0]),
        description="Diary processing tool",
    )
    subcommands = parser.add_subparsers()
    command_find_epoch = subcommands.add_parser(
        "find-epoch"
    )
    command_find_epoch.set_defaults(process=find_epoch)

    options = parser.parse_args(args[1:])
    options.process(options)


if __name__ == "__main__":
    main(sys.argv)
