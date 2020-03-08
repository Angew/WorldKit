import re

import skyfield


class Time:
    _text_parser = re.compile(r"""
        ^\s*
        # Day+month or Special day
        (?P<day>\d+)\. \s* (?P<month>\d+)\.
        |
        (?P<special_day>II?)\.

        \s*
        (?P<year>-? \s* \d+)

        \s*
        (?:
            (?P<hour>\d+) : (?P<minute>\d+)
        )?

        \s*$
    """, flags=re.VERBOSE)

    def __init__(self, year=None, month=None, day=None, hour=None, minute=None, special_day=None, ordinal_day=None, base_year=None):
        if (month is None) != (day is None):
            raise ValueError("Must specify either month and day, or neither")
        if year is None and base_year is not None:
            raise ValueError("base_year without year makes no sense")
        if hour is None and minute is not None:
            raise ValueError("minute without hour makes no sense")

        self.__ordinal_day = None
        day_spec = 0
        if month is not None:
            day_spec += 1
            self.__ordinal_day = 1 + (month-1)*28 + (day-1)
        if special_day is not None:
            day_spec += 1
            if special_day in {1, "I"}:
                self.__ordinal_day = 0
            elif special_day in {2, "II"}:
                self.__ordinal_day = 365
            else:
                raise ValueError("Invalid value of special_day")
        if ordinal_day is not None:
            day_spec += 1
            self.__ordinal_day = ordinal_day
        if day_spec > 1:
            raise ValueError("Date defined multiple times")

        self.__year = year
        self.__base_year = base_year
        self.__hour = hour
        self.__minute = minute

    @classmethod
    def parse_text(cls, text, base_year=None):
        match = cls._text_parser.match(text)
        if not match:
            raise ValueError("Invalid time format")
        def convert(name):
            res = match.group(name)
            if res is not None and name != "special_day":
                res = int(res)
            return res
        kwargs = {name: convert(name) for name in ("year", "month", "day", "special_day", "hour", "minute")}
        kwargs["base_year"] = base_year
        return cls(**kwargs)

    @classmethod
    def story(cls, text):
        return cls.parse_text(text, base_year=0)

    @property
    def hour(self):
        return self.__hour or 0

    @property
    def minute(self):
        return self.__minute or 0

    def __str__(self):
        res = []
        d = self.__ordinal_day
        if d is not None:
            if d == 0:
                res += ["I."]
            elif d == 365:
                res += ["II."]
            else:
                res += [f"{(d-1)%28 + 1}.{(d-1)/28 + 1}."]
        if self.__year is not None:
            res += [str(self.__year)]
        if self.__hour is not None:
            res += ["f{self.hour}:{self.minute}"]
        return " ".join(res)



if __name__ == "__main__":
    print(Time.story("I.7"))
