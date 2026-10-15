**To:** eckard.specht@physik.uni-magdeburg.de
**Subject:** New packings for cci, csq, scu and ssp (Vidit Patankar)

Dear Dr. Specht,

I would like to submit new packings for four categories on Packomania: circles in a circle (cci),
circles in a square (csq), spheres in a cube (scu) and spheres in a sphere (ssp). Each one has a
larger radius than the packing currently listed for that N.

The attached zip follows your formatting hints: one `<category><N>.pck` file per packing, with the
radius on the first line, my name on the second, and the coordinates after that (container size as on
Packomania: unit square or cube centred at the origin, circle or sphere of radius 1). It also contains
`SUMMARY.tsv` with the published radius, the new radius and the relative gain for every entry.

I checked every packing in exact rational arithmetic: all objects inside the container, no two
overlapping, and a radius strictly larger than the published one. The checker is included in the zip
(`verify_packings.py`) and in https://github.com/vidigoat/optimization-records/tree/main/packomania.

Many of these entries had a smaller radius than a packing with more objects listed in the same
table. I started from those larger packings with objects removed and then re-optimised each one, so
most new radii are also above what removing objects alone gives.

Thank you for maintaining Packomania.

Best regards,
Vidit Patankar
