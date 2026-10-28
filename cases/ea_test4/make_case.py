#!/usr/bin/env python3
"""EA 2D benchmark Test 4 (Neelz and Pender 2013, Environment Agency report SC120002): speed of flood propagation over
an extended floodplain. Writes the flat 1000 x 2000 m terrain and the uniform Manning n (0.05) on the 5 m grid the test
specifies; config.txt and bc_test4.csv (inflow read from the report's Figure A.8) are in this folder.

  python3 make_case.py && ../../bin/2D-DWE --config config.txt --out out      (run from this folder)

Needs numpy and rasterio. Output points (report Table A.4 / Test4Output.csv): (50,1000) (100,1000) (200,1000)
(300,1000) (400,1000) (300,1300) m, gauges P1-P6 in config.txt. See RESULT.md for the comparison.
"""
import numpy as np
import rasterio
from rasterio.transform import from_origin

NX, NY, DX = 200, 400, 5.0
TR = from_origin(0.0, 2000.0, DX, DX)
for name, val in (("dem.tif", 0.0), ("manning.tif", 0.05)):
    with rasterio.open(name, "w", driver="GTiff", width=NX, height=NY, count=1, dtype="float32",
                       crs="EPSG:27700", transform=TR) as d:
        d.write(np.full((NY, NX), val, np.float32), 1)
print("wrote dem.tif, manning.tif")
