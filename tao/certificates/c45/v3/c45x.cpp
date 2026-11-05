// c45x.cpp -- exact histogram of nu0 <= C for the Romanoff obstruction bound, by group factorisation.
//
// Same decomposition as c45g.cpp (see there): nu0(a) = sum_{s mod 4} sigma_s beta_s with
//   sigma_s = sum_{c = s (4)} prod_{G in B} alpha_G(c),   beta_s = prod_{G in A} alpha_G(s),
// independent over a. Counts are exact: per-group laws in 128 bits, the B side stored up to Aff(Z/12) and then
// sigma up to Aff(Z/4) (128-bit orbit weights), the A side in full (256-bit counts), and the histogram in 384 bits.
// Only bins nu0 <= C are produced; the rest is reported as one aggregate count (= M0 - sum of the bins), which the
// evaluator treats as the regime min(...) = 1/(2M) after checking that this is valid for C.
// A pair (sigma, beta) is skipped when sum(sigma) * min(beta) > C, since then nu0 > C.
//
//   ./c45x C threads B:3,5,7,13 ... A:29,43,113,127 ...  > hist.txt   (lines "nu0 count", then "# above count")
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
using namespace std;
typedef unsigned __int128 u128;
typedef uint64_t u64;

// ---- fixed-width unsigned integers, little-endian 64-bit limbs
template <int N> struct UInt {
    u64 w[N] = {};
    void add(const UInt& o) { u128 c = 0; for (int i = 0; i < N; i++) { c += (u128)w[i] + o.w[i]; w[i] = (u64)c; c >>= 64; } if (c) { fprintf(stderr, "overflow\n"); exit(2); } }
    template <int M> static UInt mul(const UInt<M>& a, u128 b) {   // a * b, b < 2^128
        UInt r; u64 bl[2] = {(u64)b, (u64)(b >> 64)};
        for (int j = 0; j < 2; j++) { u128 c = 0; for (int i = 0; i < M && i + j < N; i++) { c += (u128)a.w[i] * bl[j] + r.w[i + j]; r.w[i + j] = (u64)c; c >>= 64; }
            if (M + j < N) { u128 cc = c + r.w[M + j]; r.w[M + j] = (u64)cc; if (cc >> 64) { fprintf(stderr, "overflow\n"); exit(2); } } else if (c) { fprintf(stderr, "overflow\n"); exit(2); } }
        return r;
    }
    string str() const {   // decimal
        UInt t = *this; string s;
        bool zero = true; for (int i = 0; i < N; i++) if (t.w[i]) zero = false;
        if (zero) return "0";
        while (true) {
            bool z = true; for (int i = 0; i < N; i++) if (t.w[i]) z = false;
            if (z) break;
            u128 rem = 0; for (int i = N - 1; i >= 0; i--) { u128 cur = (rem << 64) | t.w[i]; t.w[i] = (u64)(cur / 10); rem = cur % 10; }
            s.push_back('0' + (int)rem);
        }
        reverse(s.begin(), s.end()); return s;
    }
};
typedef UInt<4> U256; typedef UInt<6> U384;
static U256 from128(u128 x) { U256 r; r.w[0] = (u64)x; r.w[1] = (u64)(x >> 64); return r; }

static long ord2(long q) { long x = 2 % q, d = 1; while (x != 1) { x = x * 2 % q; d++; } return d; }
struct Group { char side; vector<long> q; vector<int> d; int L, I; };
typedef vector<int> Vec;

// exact law of alpha for one group: map vector -> count (sum of counts = prod q)
static map<Vec, u128> group_law(const Group& g) {
    int L = g.L, I = g.I;
    vector<u128> cls(I, 0); u128 full = 0;
    for (int x = 0; x < L; x++) { cls[x % I] |= (u128)1 << x; full |= (u128)1 << x; }
    vector<int> idx(g.q.size()); iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return g.d[a] < g.d[b]; });
    auto km = [&](int i) { vector<u128> k(g.d[i], 0); for (int x = 0; x < L; x++) k[x % g.d[i]] |= (u128)1 << x; return k; };
    int ntail = min<int>(2, idx.size()), nhead = idx.size() - ntail;
    map<u128, u128> cur{{0, 1}};
    for (int h = 0; h < nhead; h++) {
        int i = idx[h]; u128 q = g.q[i], d = g.d[i]; auto K = km(i); map<u128, u128> nx;
        for (auto& [m, c] : cur) { nx[m] += c * (q - d); for (auto k : K) nx[m | k] += c; }
        cur.swap(nx);
    }
    vector<vector<pair<u128, u128>>> opts;
    for (int t = nhead; t < (int)idx.size(); t++) {
        int i = idx[t]; vector<pair<u128, u128>> o{{0, (u128)(g.q[i] - g.d[i])}};
        for (auto k : km(i)) o.push_back({k, 1}); opts.push_back(o);
    }
    map<Vec, u128> law; Vec a(I);
    auto emit = [&](u128 m, u128 c) {
        u128 alive = full & ~m;
        for (int s = 0; s < I; s++) { u128 v = alive & cls[s]; a[s] = __builtin_popcountll((u64)v) + __builtin_popcountll((u64)(v >> 64)); }
        law[a] += c;
    };
    for (auto& [m, c] : cur) {
        if (opts.empty()) emit(m, c);
        else if (opts.size() == 1) for (auto& [k1, c1] : opts[0]) emit(m | k1, c * c1);
        else for (auto& [k1, c1] : opts[0]) for (auto& [k2, c2] : opts[1]) emit(m | k1 | k2, c * c1 * c2);
    }
    return law;
}

static vector<vector<int>> affine(int n) {
    vector<vector<int>> G;
    for (int u = 1; u < n; u++) if (gcd(u, n) == 1) for (int t = 0; t < n; t++) { vector<int> p(n); for (int c = 0; c < n; c++) p[c] = (u * c + t) % n; G.push_back(p); }
    return G;
}
static const auto G12 = affine(12), G4 = affine(4);
static Vec canon(const Vec& v, const vector<vector<int>>& G) {
    Vec best = v, w(v.size());
    for (auto& p : G) { for (size_t c = 0; c < v.size(); c++) w[c] = v[p[c]]; if (w < best) best = w; }
    return best;
}
static u64 pack4(const array<int, 4>& s) { return ((u64)s[0] << 48) | ((u64)s[1] << 32) | ((u64)s[2] << 16) | (u64)s[3]; }
static u64 canon4(array<int, 4> s) {
    u64 best = ~0ull;
    for (auto& p : G4) { array<int, 4> w; for (int c = 0; c < 4; c++) w[c] = s[p[c]]; best = min(best, pack4(w)); }
    return best;
}

int main(int argc, char** argv) {
    long C = stol(argv[1]); int NT = stoi(argv[2]);
    vector<Group> gs;
    for (int i = 3; i < argc; i++) {
        string s = argv[i]; Group g; g.side = s[0]; g.I = (g.side == 'B') ? 12 : 4; g.L = g.I;
        stringstream ss(s.substr(2)); string tok;
        while (getline(ss, tok, ',')) { long q = stol(tok); g.q.push_back(q); g.d.push_back(ord2(q)); g.L = lcm(g.L, g.d.back()); }
        if (g.L > 128) { fprintf(stderr, "L too large\n"); return 1; }
        gs.push_back(g);
    }
    vector<Group> Bg, Ag; for (auto& g : gs) (g.side == 'B' ? Bg : Ag).push_back(g);
    // ---- B side: canonical product over all but the last B group, then the last group fused into sigma
    map<Vec, u128> B{{Vec(12, 1), 1}};
    for (size_t gi = 0; gi + 1 < Bg.size(); gi++) {
        auto law = group_law(Bg[gi]); map<Vec, u128> nx; Vec w(12);
        for (auto& [v, c] : B) for (auto& [u, cu] : law) { for (int k = 0; k < 12; k++) w[k] = v[k] * u[k]; nx[canon(w, G12)] += c * cu; }
        B.swap(nx); fprintf(stderr, "B product after group %zu: %zu reps\n", gi, B.size());
    }
    vector<pair<Vec, u128>> Bv(B.begin(), B.end()), lastv;
    { auto law = group_law(Bg.back()); lastv.assign(law.begin(), law.end()); fprintf(stderr, "last B law %zu\n", lastv.size()); }
    vector<unordered_map<u64, u128>> part(NT);
    {
        vector<thread> th;
        for (int t = 0; t < NT; t++) th.emplace_back([&, t]() {
            auto& M = part[t]; array<int, 4> s;
            for (size_t i = t; i < Bv.size(); i += NT) {
                auto& v = Bv[i].first; u128 c = Bv[i].second;
                for (auto& [u, cu] : lastv) {
                    for (int k = 0; k < 4; k++) s[k] = v[k] * u[k] + v[k + 4] * u[k + 4] + v[k + 8] * u[k + 8];
                    M[canon4(s)] += c * cu;
                }
            }
        });
        for (auto& x : th) x.join();
    }
    unordered_map<u64, u128> S4; for (auto& M : part) { for (auto& [k, c] : M) S4[k] += c; M.clear(); }
    fprintf(stderr, "sigma reps %zu\n", S4.size());
    // ---- A side: full product law, 256-bit counts
    map<array<int, 4>, U256> A; A[{1, 1, 1, 1}] = from128(1);
    for (auto& g : Ag) {
        auto law = group_law(g); map<array<int, 4>, U256> nx;
        for (auto& [v, c] : A) for (auto& [u, cu] : law) { array<int, 4> w{v[0] * u[0], v[1] * u[1], v[2] * u[2], v[3] * u[3]}; nx[w].add(U256::mul<4>(c, cu)); }
        A.swap(nx); fprintf(stderr, "A product: %zu vectors\n", A.size());
    }
    struct BE { array<int, 4> v; int mn; U256 c; };
    vector<BE> Av; for (auto& [v, c] : A) Av.push_back({v, *min_element(v.begin(), v.end()), c});
    sort(Av.begin(), Av.end(), [](const BE& a, const BE& b) { return a.mn < b.mn; });
    struct SE { array<int, 4> s; long sum; u128 w; };
    vector<SE> Sv;
    for (auto& [k, w] : S4) { array<int, 4> s{(int)(k >> 48), (int)((k >> 32) & 0xffff), (int)((k >> 16) & 0xffff), (int)(k & 0xffff)}; Sv.push_back({s, (long)s[0] + s[1] + s[2] + s[3], w}); }
    // ---- pairing with the cut
    vector<vector<U384>> hp(NT, vector<U384>(C + 1));
    {
        vector<thread> th;
        for (int t = 0; t < NT; t++) th.emplace_back([&, t]() {
            auto& h = hp[t];
            for (size_t i = t; i < Sv.size(); i += NT) {
                auto& S = Sv[i];
                for (auto& b : Av) {
                    if (S.sum > 0 && (long)b.mn * S.sum > C) break;
                    long n = (long)S.s[0] * b.v[0] + (long)S.s[1] * b.v[1] + (long)S.s[2] * b.v[2] + (long)S.s[3] * b.v[3];
                    if (n <= C) h[n].add(U384::mul<4>(b.c, S.w));
                }
            }
        });
        for (auto& x : th) x.join();
    }
    for (int t = 1; t < NT; t++) for (long n = 0; n <= C; n++) hp[0][n].add(hp[t][n]);
    for (long n = 0; n <= C; n++) { bool z = true; for (int i = 0; i < 6; i++) if (hp[0][n].w[i]) z = false; if (!z) printf("%ld %s\n", n, hp[0][n].str().c_str()); }
    fprintf(stderr, "done\n");
}
