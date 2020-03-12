import argparse
import os.path
import sys


def find_epoch(options):
    pass


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
