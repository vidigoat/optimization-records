/* Exact minimum distance of a linear [n,k] code over GF(q), q in {2,3,4,5,7,8,9? no: 2,3,4,5,7,8}, by
   enumerating every message. Trusted checker: no heuristics, no early exit except when a codeword of
   weight below the target is found (then the code is certainly not good enough).

   stdin: q n k target, then k*n integers 0..q-1 (row-major generator matrix).
   GF(4) and GF(8) elements are bit vectors in the polynomial basis (alpha^2 = alpha + 1 for GF(4),
   alpha^3 = alpha + 1 for GF(8)), so addition is XOR; only additions are needed for enumeration once the
   rows alpha^j * g_i are precomputed.
   stdout: the minimum distance (or a value below target if one was found), and the codeword count. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int q, n, k, target;

static int gfmul(int a, int b) {             /* only used for q = 4, 8 (char 2) */
    int m = (q == 4) ? 2 : 3, poly = (q == 4) ? 0x7 : 0xB, r = 0;
    for (int i = 0; i < m; i++) if (b >> i & 1) r ^= a << i;
    for (int i = 2 * m - 2; i >= m; i--) if (r >> i & 1) r ^= poly << (i - m);
    return r;
}

int main(void) {
    if (scanf("%d %d %d %d", &q, &n, &k, &target) != 4) return 2;
    if (n < 1 || n > 256 || k < 1 || k > 40) return 2;
    int *G = malloc(sizeof(int) * k * n);
    for (int i = 0; i < k * n; i++) { if (scanf("%d", &G[i]) != 1 || G[i] < 0 || G[i] >= q) return 2; }
    int best = n + 1;
    unsigned long long count = 0;
    if (q == 2 || q == 4 || q == 8) {
        int m = (q == 2) ? 1 : (q == 4) ? 2 : 3, W = (n + 63) / 64, R = k * m;
        /* row (i, j) = alpha^j * g_i, stored as m bit-planes of n bits */
        uint64_t *rows = calloc((size_t)R * m * W, 8), *c = calloc((size_t)m * W, 8);
        for (int i = 0; i < k; i++) for (int j = 0; j < m; j++) {
            int r = i * m + j, aj = 1 << j;
            for (int col = 0; col < n; col++) {
                int v = (q == 2) ? G[i * n + col] : gfmul(G[i * n + col], aj);
                for (int b = 0; b < m; b++) if (v >> b & 1) rows[((size_t)r * m + b) * W + col / 64] |= 1ULL << (col % 64);
            }
        }
        unsigned long long total = 1ULL << R;
        for (unsigned long long s = 1; s < total; s++) {
            int r = __builtin_ctzll(s);
            for (int b = 0; b < m; b++) for (int w = 0; w < W; w++) c[b * W + w] ^= rows[((size_t)r * m + b) * W + w];
            int wt = 0;
            for (int w = 0; w < W; w++) { uint64_t o = 0; for (int b = 0; b < m; b++) o |= c[b * W + w]; wt += __builtin_popcountll(o); }
            count++;
            if (wt < best) { best = wt; if (best < target) break; }
        }
    } else {                                   /* prime q: reflected q-ary Gray code, one row added per step */
        int *digit = calloc(k, sizeof(int)), *dir = malloc(sizeof(int) * k), *c = calloc(n, sizeof(int));
        for (int i = 0; i < k; i++) dir[i] = 1;
        unsigned long long total = 1;
        for (int i = 0; i < k; i++) total *= q;
        for (unsigned long long s = 1; s < total; s++) {
            int i = 0;                           /* find the digit to move: lowest one not at its end */
            while ((dir[i] == 1 && digit[i] == q - 1) || (dir[i] == -1 && digit[i] == 0)) { dir[i] = -dir[i]; i++; }
            digit[i] += dir[i];
            int sgn = dir[i] == 1 ? 1 : q - 1, wt = 0;
            for (int col = 0; col < n; col++) { c[col] = (c[col] + sgn * G[i * n + col]) % q; wt += c[col] != 0; }
            count++;
            if (wt < best) { best = wt; if (best < target) break; }
        }
    }
    printf("%d %llu\n", best, count);
    return 0;
}
