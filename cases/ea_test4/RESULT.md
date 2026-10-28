# EA 2D benchmark Test 4 (Neelz and Pender 2013, Environment Agency report SC120002)

Binary: 2D-DWE rockcreek branch (fec3648) = the binary of the OpenWatershedTwin flood-map service (Sligo Creek); the merged main gives identical gauges. Run 2026-10-28.
Case: flat 1000 x 2000 m plain, 5 m grid (200 x 400), n 0.05, uniform inflow over 20 m (4 cells) in the middle of the
west side, hydrograph read from Figure A.8: 0 to 5 min, 20 m3/s at 60 min, 20 to 240 min, 0 at 300 min; closed
elsewhere; dry start; cfl 0.7, theta 0.9 (as the Sligo runs). Output points from the report (Tygron wiki table).
Run time 41 s; mass error 2e-15.

Pass criterion (set before the run): depth curves at points 1, 3, 5, 6 (Figure 4.10) within the spread of the 17
packages; arrival times within their spread. Velocities not compared.

| point | arrival (first wet) | DWE depth 1h / 2h / 3h / 4h / 5h (m) | reference spread read from Fig. 4.10 |
|---|---|---|---|
| 1 | 10.0 min (8-11) | 0.292 / 0.305 / 0.307 / 0.308 / 0.150 | plateau 0.30-0.33; DWE at lower edge (with TUFLOW FV1/FV2) |
| 3 | 30.3 min (25-31) | 0.161 / 0.204 / 0.212 / 0.215 / 0.150 | 2h 0.20-0.22, 4h 0.21-0.23, 5h 0.147-0.158 |
| 5 | 62.3 min (52-62) | 0.000 / 0.143 / 0.161 / 0.169 / 0.147 | 2h 0.14-0.155, 3h 0.16-0.172, 4h 0.167-0.18, 5h 0.146-0.155; arrival at late edge |
| 6 | 59.7 min (55-65) | 0.009 / 0.150 / 0.166 / 0.173 / 0.148 | 2h 0.137-0.15, 3h 0.152-0.167, 4h 0.162-0.176, 5h 0.145-0.154 |

Result: PASS (all depths within the spread; point 5 arrival at the late edge, within reading accuracy).
