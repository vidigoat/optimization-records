**To:** Markus Grassl (contact address on http://www.codetables.de)
**Subject:** Four new linear codes over GF(5): [27,10,13], [64,10,39], [81,10,51], [82,10,52]

Dear Dr. Grassl,

I would like to submit four linear codes over GF(5) that improve the lower bounds in your tables:

| Code | Current lower bound | New |
|---|---|---|
| [27, 10] | 12 | 13 |
| [64, 10] | 38 | 39 |
| [81, 10] | 50 | 51 |
| [82, 10] | 51 | 52 |

The generator matrices are attached (one row per line, entries 0–4) and are also at
https://github.com/vidigoat/optimization-records/tree/main/codetables. Each minimum distance was
computed exactly by enumerating all 5^10 codewords, and confirmed with a second, independent
enumeration. The codes were found by Construction X applied to pairs of nested constacyclic codes.

Thank you for maintaining the tables.

Best regards,
Vidit Patankar
