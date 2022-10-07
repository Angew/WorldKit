import os
import skyfield.almanac
import skyfield.api
import sys

script_dir = sys.path[0]

timescale = skyfield.api.load.timescale(builtin=True)
loader = skyfield.api.Loader(os.path.join(script_dir, "skyfield_data", "loaded"))
eph = loader("de406.bsp")

my_time = timescale.tt_jd(2110762.5)

earth = eph["Earth"]
moon = eph["Moon"]
m = earth.at(my_time).observe(moon)
point = skyfield.api.wgs84.latlon_of(m)

print(point)
