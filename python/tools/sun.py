import math


# Copied from https://en.wikipedia.org/wiki/Position_of_the_Sun:
# declination = - arcsin( sin(-23.44d) * cos( 360d/365.24 * M + 360d/PI * 0.0167 * sin( 360d/365.24 * (M-12) ) ) )

declination_per_day = math.radians(360)/365.24
eccentricity = 0.0167
# perihelion_day = 12
perihelion_day = 0 # By DM fiat, perihelion occurs on solstice noon (it happened in 1246 AD, so what)
solstice_declination = math.radians(-23.44)

def declination(time):
    return math.degrees(math.asin(
        math.sin(solstice_declination) *
        math.cos(
            declination_per_day * time +
            math.radians(360)/math.pi * eccentricity * math.sin(declination_per_day * (time - perihelion_day))
        )
    ))


minutes_per_day = 24 * 60

def time_of_declination(target_declination, time_guess, interval_days=3):
    best_error = abs(declination(time_guess) - target_declination)
    best_time = time_guess
    minute = - interval_days * minutes_per_day
    while minute <= interval_days * minutes_per_day:
        time = time_guess + minute/minutes_per_day
        error = abs(declination(time) - target_declination)
        if error < best_error:
            best_error = error
            best_time = time
        minute += 1
    return best_time


if __name__ == '__main__':
    print("Winter solstice:", time_of_declination(-23.44, 0))
    print("Spring equinox:", time_of_declination(0, 91.25))
    print("Summer solstice:", time_of_declination(23.44, 182.5))
    print("Autumn equinox:", time_of_declination(0, 3*91.25))
    print("Winter solstice:", time_of_declination(-23.44, 365))
