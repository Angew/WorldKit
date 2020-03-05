import re
import skyfield


class Time:
    _textParser = re.compile(r"""
        ^\s*
        # Day+month or Special day
        (\d+)\. \s* (\d+)\. # Day+month
        |
        (II?)\. # Special day

        \s*
        # Year
        (-? \s* \d+)

        \s*
        # Time
        (?:
            (\d+) : (\d+)
        )?

        \s*$
    """)

    @classmethod
    def parseText(cls, text):
        pass

    @classmethod
    def story(cls, text):
        pass
