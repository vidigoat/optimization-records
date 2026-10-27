// histB.cpp -- independent exact base-histogram computer ("checker B").
//
// Problem.  P = set of odd primes q with d_q = ord_q(2), T0 = lcm d_q, M0 = prod q.
// For a mod M0 let nu(a) = #{ r in Z_T0 : a mod q != 2^r mod q for all q }.
// Output delta(nu) = #{ a mod M0 : nu(a) = nu } exactly, as lines "nu delta".
//
// Reformulation.  For each q independently, a mod q is either a non-power of 2
// (q - d_q residues; removes nothing) or equals 2^s for exactly one s in Z_{d_q}
// (removes the class r = s mod d_q).  By CRT, delta is the weighted distribution
// of |Z_T0 minus union of chosen classes| over independent per-prime choices.
//
// Dynamic program ("cut" DP).  Process the primes in some order.  After k primes
// let Lp = lcm(processed orders), Lr = lcm(remaining orders), g = gcd(Lp, Lr).
// State: the vector c on Z_g,
//      c[x] = #{ r in Z_T0 : r = x mod g and r survives the processed primes }.
// For r uniform in Z_T0 and fixed r mod g, the residues r mod Lp and r mod Lr are
// independent (the pair runs over all compatible pairs equally often), hence c
// is a complete description of the past as far as the future is concerned.
// Transition with a prime of order d (d | Lr): G = lcm(g, d); expand
// e[X] = c[X mod g] / (G/g) on Z_G (exact); zero the class X = s mod d if a class
// is chosen; project to the next cut modulus g' (g' | G) by summing fibres.
// At the end g = 1 and c[0] = nu.
//
// Symmetry.  Let D* = {gcd(d, g) : d a remaining order}.  If sigma is a
// permutation of Z_g preserving every partition "x mod e" (e in D*), it lifts to a
// permutation of Z_Lr preserving every remaining class: write X in Z_Lr, per prime
// l, as X_l = x_l + l^{a_l} h_l (x_l in [0, l^{a_l}), a_l = v_l(g)); put
// sigma~(X) = (sigma(X mod g), same h).  For d remaining, X = Y mod d iff
// X = Y mod gcd(d, g) in Z_g and the h-digits agree up to level v_l(d) - a_l, and
// both conditions are preserved.  Since the future weight of a class choice is 1
// for every class, substituting X -> sigma~(X) shows that c and c o sigma have the
// same future nu-distribution.  Canonical forms used (all orbit invariants that
// return a member of the orbit): if D* is a divisibility chain, recursive sorting
// of the tree of nested blocks; otherwise a componentwise subgroup (per l-component
// a tree group with the occurring valuations as levels, the largest component
// sorted, the others enumerated) or the affine group x -> u x + t.  The last few
// primes are enumerated depth-first without merging (last two in closed form).
//
// Exact arithmetic: weights are fixed-width multi-limb unsigned integers (128,
// 192 or 256 bits, chosen from M0).  At the end the following are verified:
//   sum delta = M0, sum delta*nu = T0*phi(M0),
//   sum delta*nu^2 and sum delta*nu^3 against direct formulas.
//
// Build:  g++ -O3 -march=native -std=c++17 -pthread histB.cpp -o histB
// Usage:  histB [options] q1 q2 ...      (prints "nu delta" lines to stdout)
//   --threads=N        worker threads (default: hardware concurrency)
//   --order=auto|given|q,q,q,...   processing order (default auto)
//   --tail=K           number of final primes enumerated without merging (default 3)
//   --canon=tree|affine|none       symmetry used for merging (default tree)
//   --maxstates=N      abort (exit 3) if a stored layer exceeds N states
//   --brute            also run a direct enumeration over a mod M0 and compare
//   --bruteonly        only run the direct enumeration
//   --selftest         built-in brute-force comparisons on small sets
//   --plan             print the plan (order, cut moduli g_k) and exit
//   --quiet            do not print per-step progress
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include <memory>
#include <chrono>
#include <numeric>
#include <functional>
#include <sys/resource.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

// ---------------------------------------------------------------- big ints
template <int NL>
struct UBig {
    u64 l[NL];
    UBig() { for (int i = 0; i < NL; i++) l[i] = 0; }
    explicit UBig(u64 v) { l[0] = v; for (int i = 1; i < NL; i++) l[i] = 0; }
    bool iszero() const { for (int i = 0; i < NL; i++) if (l[i]) return false; return true; }
    u64 add(const UBig &o) {  // returns carry-out (nonzero = overflow)
        u64 c = 0;
        for (int i = 0; i < NL; i++) {
            u128 s = (u128)l[i] + o.l[i] + c;
            l[i] = (u64)s; c = (u64)(s >> 64);
        }
        return c;
    }
    u64 mul_small(u64 m) {  // *= m, returns overflow
        u64 c = 0;
        for (int i = 0; i < NL; i++) {
            u128 p = (u128)l[i] * m + c;
            l[i] = (u64)p; c = (u64)(p >> 64);
        }
        return c;
    }
    u64 divmod_small(u64 m) {  // /= m, returns remainder
        u128 r = 0;
        for (int i = NL - 1; i >= 0; i--) {
            u128 cur = (r << 64) | l[i];
            l[i] = (u64)(cur / m); r = cur % m;
        }
        return (u64)r;
    }
    bool operator==(const UBig &o) const { for (int i = 0; i < NL; i++) if (l[i] != o.l[i]) return false; return true; }
    std::string str() const {
        UBig t = *this;
        if (t.iszero()) return "0";
        std::vector<u64> parts;
        while (!t.iszero()) parts.push_back(t.divmod_small(10000000000000000000ULL));
        std::string s = std::to_string(parts.back());
        char buf[32];
        for (int i = (int)parts.size() - 2; i >= 0; i--) { snprintf(buf, sizeof buf, "%019llu", (unsigned long long)parts[i]); s += buf; }
        return s;
    }
    template <int NL2> UBig<NL2> widen() const {
        UBig<NL2> r;
        for (int i = 0; i < NL; i++) { if (i < NL2) r.l[i] = l[i]; else if (l[i]) { fprintf(stderr, "widen overflow\n"); exit(2); } }
        return r;
    }
};
typedef UBig<6> WC;  // check accumulators (moments up to M0*T0^3), and output histogram

static void die(const char *msg) { fprintf(stderr, "FATAL: %s\n", msg); exit(2); }

// ---------------------------------------------------------------- number theory
static u64 gcdu(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
static u64 lcmu(u64 a, u64 b) { return a / gcdu(a, b) * b; }
static bool isprime(u64 n) { if (n < 2) return false; for (u64 p = 2; p * p <= n; p++) if (n % p == 0) return false; return true; }
static u64 ord2(u64 q) {  // multiplicative order of 2 mod q, by direct iteration
    u64 x = 2 % q, k = 1;
    while (x != 1) { x = x * 2 % q; k++; if (k > q) die("order loop"); }
    return k;
}
static std::vector<std::pair<u64, int>> factor(u64 n) {
    std::vector<std::pair<u64, int>> f;
    for (u64 p = 2; p * p <= n; p++) if (n % p == 0) { int a = 0; while (n % p == 0) { n /= p; a++; } f.push_back({p, a}); }
    if (n > 1) f.push_back({n, 1});
    return f;
}
static int val(u64 n, u64 p) { int a = 0; while (n % p == 0) { n /= p; a++; } return a; }
static u64 ipow(u64 b, int e) { u64 r = 1; while (e--) r *= b; return r; }

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static double peak_rss_gb() {
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return ru.ru_maxrss / 1e9;          // bytes on macOS
#else
    return ru.ru_maxrss * 1024.0 / 1e9;  // kilobytes on Linux
#endif
}
static std::string join(const std::vector<u64> &v) { std::string s; for (size_t i = 0; i < v.size(); i++) { if (i) s += ","; s += std::to_string(v[i]); } return s; }

// ---------------------------------------------------------------- problem setup
struct Problem {
    std::vector<u64> q, d;  // in processing order
    u64 T0 = 1;
    int n = 0;
    std::vector<u64> g;     // g[k], k = 0..n
};

static std::vector<u64> compute_g(const std::vector<u64> &d) {
    int n = d.size();
    std::vector<u64> pre(n + 1, 1), suf(n + 1, 1), g(n + 1);
    for (int k = 0; k < n; k++) pre[k + 1] = lcmu(pre[k], d[k]);
    for (int k = n - 1; k >= 0; k--) suf[k] = lcmu(suf[k + 1], d[k]);
    for (int k = 0; k <= n; k++) g[k] = gcdu(pre[k], suf[k]);
    return g;
}

// Automatic processing order (affects speed/memory only, never correctness).
// Heuristic found by experiment: (1) take the primes whose order carries the
// highest power of 2 dividing T0, then the next power of 2 (largest first);
// (2) then, for each odd prime l >= 5 dividing T0, the primes whose order is
// divisible by l, as a contiguous block, larger l first; (3) then the remaining
// primes whose order is not divisible by 9 (largest order first); (4) finally
// the primes whose order is divisible by 9, smallest order first.
static std::vector<int> auto_order(const std::vector<u64> &d0) {
    int n = d0.size();
    u64 T0 = 1; for (u64 d : d0) T0 = lcmu(T0, d);
    std::vector<int> used(n, 0), ord;
    auto take = [&](std::vector<int> idx, bool asc) {
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return asc ? (d0[a] < d0[b] || (d0[a] == d0[b] && a < b)) : (d0[a] > d0[b] || (d0[a] == d0[b] && a < b)); });
        for (int i : idx) if (!used[i]) { used[i] = 1; ord.push_back(i); }
    };
    int e2 = val(T0, 2);
    for (int e = e2; e >= 3; e--) {
        std::vector<int> idx;
        for (int i = 0; i < n; i++) if (!used[i] && val(d0[i], 2) == e) idx.push_back(i);
        take(idx, false);
    }
    auto f = factor(T0);
    std::vector<u64> ls;
    for (auto &pa : f) if (pa.first >= 5) ls.push_back(pa.first);
    std::sort(ls.rbegin(), ls.rend());
    for (u64 l : ls) {
        std::vector<int> idx;
        for (int i = 0; i < n; i++) if (!used[i] && d0[i] % l == 0) idx.push_back(i);
        take(idx, true);
    }
    {
        std::vector<int> idx;
        for (int i = 0; i < n; i++) if (!used[i] && d0[i] % 9 != 0) idx.push_back(i);
        take(idx, false);
    }
    {
        std::vector<int> idx;
        for (int i = 0; i < n; i++) if (!used[i]) idx.push_back(i);
        take(idx, true);
    }
    return ord;
}

// ---------------------------------------------------------------- canonicaliser
// Component Z_{l^a} of Z_g with tree levels lev = {0 = j_0 < j_1 < ... < j_L = a}.
struct Comp {
    u64 l = 1; int a = 0; int N = 1;
    std::vector<int> lev;
    double logsize() const {  // log |tree automorphism group|
        double s = 0;
        for (size_t k = 1; k < lev.size(); k++) {
            double nodes = (double)ipow(l, lev[k - 1]);
            u64 b = ipow(l, lev[k] - lev[k - 1]);
            s += nodes * std::lgamma((double)b + 1);
        }
        return s;
    }
};

// All elements of the tree automorphism group of a component (as permutations of Z_N).
// Leaf x <-> digits (i_1..i_L), x = sum_k i_k l^{j_{k-1}}, i_k in [0, b_k).
// An automorphism applies a permutation pi_node to i_k at every internal node,
// the node being identified by the (source) prefix residue.
static std::vector<std::vector<uint16_t>> tree_elems(const Comp &c) {
    int L = c.lev.size() - 1;
    std::vector<u64> b(L + 1), w(L + 1);
    for (int k = 1; k <= L; k++) { b[k] = ipow(c.l, c.lev[k] - c.lev[k - 1]); w[k] = ipow(c.l, c.lev[k - 1]); }
    std::vector<std::pair<int, u64>> nodes;  // (k, y): internal node at depth k-1 with residue y mod l^{j_{k-1}}
    for (int k = 1; k <= L; k++) for (u64 y = 0; y < w[k]; y++) nodes.push_back({k, y});
    std::vector<std::vector<std::vector<int>>> perms(L + 1);
    for (int k = 1; k <= L; k++) {
        std::vector<int> p(b[k]); std::iota(p.begin(), p.end(), 0);
        do perms[k].push_back(p); while (std::next_permutation(p.begin(), p.end()));
    }
    std::vector<std::vector<uint16_t>> out;
    std::vector<size_t> choice(nodes.size(), 0);
    while (true) {
        std::vector<uint16_t> pm(c.N);
        for (int x = 0; x < c.N; x++) {
            u64 y = 0, yy = 0;
            for (int k = 1; k <= L; k++) {
                u64 ik = (x / w[k]) % b[k];
                size_t ni = 0;
                for (int kk = 1; kk < k; kk++) ni += w[kk];
                ni += y;
                u64 jk = perms[k][choice[ni]][ik];
                y += ik * w[k]; yy += jk * w[k];
            }
            pm[x] = (uint16_t)yy;
        }
        out.push_back(pm);
        size_t i = 0;
        while (i < nodes.size()) { if (++choice[i] < perms[nodes[i].first].size()) break; choice[i] = 0; i++; }
        if (i == nodes.size()) break;
    }
    return out;
}
static std::vector<std::vector<uint16_t>> affine_elems(const Comp &c) {
    std::vector<std::vector<uint16_t>> out;
    for (int u = 1; u <= c.N; u++) if (gcdu(u, c.N) == 1)
        for (int t = 0; t < c.N; t++) {
            std::vector<uint16_t> pm(c.N);
            for (int x = 0; x < c.N; x++) pm[x] = (uint16_t)(((u64)u * x + t) % c.N);
            out.push_back(pm);
        }
    return out;
}

struct CanonScratch {
    std::vector<uint16_t> G, S, best, cc;
    std::vector<int> order, tmp;
};

struct Canon {
    int mode = 0;  // 0 none, 1 affine (whole Z_g), 2 tree
    int g = 1;
    std::vector<std::vector<uint16_t>> mul;  // affine mode
    int N = 1, R = 1, ne = 1, L = 0;          // tree mode
    std::vector<uint16_t> Tab;      // ne * g : source index for (e, x*R + r)
    std::vector<uint16_t> outIdx;   // g : standard index for serial position p*R + r
    std::vector<u64> b, w;          // tree branching b[k], w[k] = l^{j_{k-1}}
    std::vector<int> nodesAt;       // number of nodes at depth k (k=0..L): l^{j_k}
    std::vector<int> unit;          // serial size of a depth-k subtree
    std::string desc;

    void init(int g_, const std::vector<u64> &rem, int mode_) {
        g = g_; mode = mode_;
        if (g == 1) { mode = 0; desc = "trivial"; return; }
        if (mode == 0) { desc = "none"; return; }
        if (mode == 1) {
            mul.clear();
            for (int u = 1; u <= g; u++) if (gcdu(u, g) == 1) {
                std::vector<uint16_t> m(g);
                for (int x = 0; x < g; x++) m[x] = (uint16_t)((u64)u * x % g);
                mul.push_back(m);
            }
            desc = "affine(" + std::to_string(g * mul.size()) + ")";
            return;
        }
        // D* = {gcd(d, g) : d remaining} U {g}: the symmetry group is the group of all
        // permutations of Z_g preserving the partitions "mod e", e in D* (see README).
        std::vector<u64> Ds;
        for (u64 d : rem) { u64 e = gcdu(d, g); if (e > 1) Ds.push_back(e); }
        Ds.push_back(g);
        std::sort(Ds.begin(), Ds.end()); Ds.erase(std::unique(Ds.begin(), Ds.end()), Ds.end());
        bool chain = true;
        for (size_t i = 1; i < Ds.size(); i++) if (Ds[i] % Ds[i - 1]) chain = false;
        if (chain) {
            // Aut = iterated wreath product (automorphisms of the tree of nested blocks)
            std::vector<u64> E(1, 1); for (u64 e : Ds) E.push_back(e);
            N = g; R = 1; ne = 1;
            Tab.resize(g); for (int p = 0; p < g; p++) Tab[p] = (uint16_t)p;
            std::vector<int> idx(g); std::iota(idx.begin(), idx.end(), 0);
            setup_tree(E, idx);
            desc = "chain[";
            for (size_t i = 0; i < E.size(); i++) { if (i) desc += "|"; desc += std::to_string(E[i]); }
            desc += "]";
            return;
        }
        std::vector<Comp> comps;
        for (auto &pa : factor(g)) {
            Comp c; c.l = pa.first; c.a = pa.second; c.N = (int)ipow(c.l, c.a);
            c.lev.push_back(0);
            for (u64 d : rem) { int v = val(d, c.l); if (v > 0 && v < c.a) c.lev.push_back(v); }
            c.lev.push_back(c.a);
            std::sort(c.lev.begin(), c.lev.end()); c.lev.erase(std::unique(c.lev.begin(), c.lev.end()), c.lev.end());
            comps.push_back(c);
        }
        int star = 0;
        for (size_t i = 0; i < comps.size(); i++) if (comps[i].logsize() > comps[star].logsize()) star = i;
        const double LIMIT = 3000;
        std::vector<int> kind(comps.size(), 2);  // 2 full tree group, 1 affine, 0 trivial
        auto gsize = [&](size_t i) -> double {
            if (kind[i] == 2) return std::exp(comps[i].logsize());
            if (kind[i] == 1) { int nu = 0; for (int u = 1; u <= comps[i].N; u++) if (gcdu(u, comps[i].N) == 1) nu++; return (double)nu * comps[i].N; }
            return 1;
        };
        while (true) {
            double prod = 1; int worst = -1; double ws = 0;
            for (size_t i = 0; i < comps.size(); i++) if ((int)i != star) { double s = gsize(i); prod *= s; if (s > ws) { ws = s; worst = i; } }
            if (prod <= LIMIT || worst < 0) break;
            if (kind[worst] == 2) { double full = gsize(worst); kind[worst] = 1; if (gsize(worst) >= full) kind[worst] = 0; }
            else kind[worst] = 0;
        }
        {
            double lg = 0; int nunits = 0;
            for (size_t i = 0; i < comps.size(); i++) lg += comps[i].logsize();
            for (int u = 1; u <= g; u++) if (gcdu(u, g) == 1) nunits++;
            if (std::fabs(lg - std::log((double)g * nunits)) < 1e-9) {  // the group is Aff(Z_g) itself
                init(g, rem, 1);
                return;
            }
        }
        desc = "";
        for (size_t i = 0; i < comps.size(); i++) {
            char buf[160];
            std::string lv; for (size_t k = 0; k < comps[i].lev.size(); k++) { if (k) lv += "|"; lv += std::to_string(ipow(comps[i].l, comps[i].lev[k])); }
            snprintf(buf, sizeof buf, "%sZ%d[%s]:%s", i ? "x" : "", comps[i].N, lv.c_str(), (int)i == star ? "sort" : kind[i] == 2 ? "tree" : kind[i] == 1 ? "aff" : "id");
            desc += buf;
        }
        std::vector<std::vector<std::vector<uint16_t>>> el(comps.size());
        for (size_t i = 0; i < comps.size(); i++) {
            if ((int)i == star) continue;
            if (kind[i] == 2) el[i] = tree_elems(comps[i]);
            else if (kind[i] == 1) el[i] = affine_elems(comps[i]);
            else { std::vector<uint16_t> id(comps[i].N); std::iota(id.begin(), id.end(), 0); el[i].push_back(id); }
        }
        const Comp &S = comps[star];
        N = S.N; R = g / N;
        // component tuple <-> X
        std::vector<int> lookup;
        size_t tot = 1; for (auto &c : comps) tot *= c.N;
        lookup.assign(tot, -1);
        auto keyof = [&](const std::vector<u64> &xs) { size_t key = 0; for (size_t i = 0; i < comps.size(); i++) key = key * comps[i].N + xs[i]; return key; };
        std::vector<u64> xs(comps.size());
        for (int X = 0; X < g; X++) { for (size_t i = 0; i < comps.size(); i++) xs[i] = X % comps[i].N; lookup[keyof(xs)] = X; }
        // idx[x*][r]: X with X = x* mod N and X = r mod R
        std::vector<int> idx(g);
        for (int X = 0; X < g; X++) idx[(X % N) * R + (X % R)] = X;
        ne = 1;
        for (size_t i = 0; i < comps.size(); i++) if ((int)i != star) ne *= el[i].size();
        Tab.assign((size_t)ne * g, 0);
        std::vector<size_t> ch(comps.size(), 0);
        for (int e = 0; e < ne; e++) {
            for (int p = 0; p < g; p++) {
                int X = idx[p];
                for (size_t i = 0; i < comps.size(); i++) { xs[i] = X % comps[i].N; if ((int)i != star) xs[i] = el[i][ch[i]][xs[i]]; }
                Tab[(size_t)e * g + p] = (uint16_t)lookup[keyof(xs)];
            }
            for (size_t i = 0; i < comps.size(); i++) {
                if ((int)i == star) continue;
                if (++ch[i] < el[i].size()) break;
                ch[i] = 0;
            }
        }
        std::vector<u64> E;
        for (int j : S.lev) E.push_back(ipow(S.l, j));
        setup_tree(E, idx);
    }

    // Tree on Z_N with nested block moduli E = {1 = e_0 | e_1 | ... | e_L = N}:
    // node at depth k = residue y mod e_k, its children y + e_k t, t < e_{k+1}/e_k.
    // Serial (DFS) position p, digits (i_1..i_L) with i_1 most significant, is the
    // leaf x = sum_k i_k e_{k-1}.  Tab is re-indexed so that gathering produces the
    // serial layout directly; unit[k] = size (in entries) of a depth-k subtree.
    void setup_tree(const std::vector<u64> &E, const std::vector<int> &idx) {
        L = E.size() - 1;
        b.assign(L + 1, 1); w.assign(L + 1, 1); nodesAt.assign(L + 1, 1);
        for (int k = 1; k <= L; k++) { b[k] = E[k] / E[k - 1]; w[k] = E[k - 1]; }
        for (int k = 0; k <= L; k++) nodesAt[k] = (int)E[k];
        std::vector<int> lay(N);
        outIdx.assign(g, 0);
        for (int p = 0; p < N; p++) {
            u64 rem2 = p, x = 0;
            for (int k = L; k >= 1; k--) { u64 ik = rem2 % b[k]; rem2 /= b[k]; x += ik * w[k]; }
            lay[p] = (int)x;
            for (int r = 0; r < R; r++) outIdx[p * R + r] = (uint16_t)idx[x * R + r];
        }
        std::vector<uint16_t> T2(Tab.size());
        for (int e = 0; e < ne; e++)
            for (int p = 0; p < N; p++)
                for (int r = 0; r < R; r++) T2[(size_t)e * g + p * R + r] = Tab[(size_t)e * g + lay[p] * R + r];
        Tab.swap(T2);
        unit.assign(L + 1, R);
        for (int k = L - 1; k >= 0; k--) unit[k] = unit[k + 1] * (int)b[k + 1];
    }

    void prep(CanonScratch &s) const {
        s.G.resize(g); s.S.resize(g); s.best.resize(g); s.cc.resize(2 * g);
        s.order.resize(std::max(N, 1)); s.tmp.resize(std::max(N, 1));
    }

    // out = canonical representative of the orbit of c (a fixed function of the orbit,
    // and itself a member of the orbit).  Tree mode: for every enumerated element e of
    // the non-sorted components, gather c into the serial layout, then sort sibling
    // subtrees bottom-up (any fixed total order on their serialisations works), and
    // keep the smallest result over e.
    void apply(const uint16_t *c, uint16_t *out, CanonScratch &s) const {
        if (mode == 0) { memcpy(out, c, 2 * g); return; }
        if (mode == 1) { apply_affine(c, out, s); return; }
        bool have = false;
        uint16_t *Sv = s.S.data(), *best = s.best.data(), *tmp = s.G.data();
        int *ord = s.order.data();
        for (int e = 0; e < ne; e++) {
            const uint16_t *tb = &Tab[(size_t)e * g];
            for (int p = 0; p < g; p++) Sv[p] = c[tb[p]];
            for (int k = L; k >= 1; k--) {
                const int u = unit[k], nb = (int)b[k];
                const int blk = u * nb, nblk = g / blk;
                if (u == 1) {  // scalar units: plain sort of each block
                    for (int B = 0; B < nblk; B++) {
                        uint16_t *v = Sv + (size_t)B * blk;
                        if (nb <= 32) { for (int i = 1; i < nb; i++) { uint16_t x = v[i]; int j = i - 1; while (j >= 0 && v[j] > x) { v[j + 1] = v[j]; j--; } v[j + 1] = x; } }
                        else std::sort(v, v + nb);
                    }
                    continue;
                }
                const size_t ub = (size_t)u * 2;
                for (int B = 0; B < nblk; B++) {
                    uint16_t *v = Sv + (size_t)B * blk;
                    if (nb == 2) {
                        if (memcmp(v, v + u, ub) > 0) { memcpy(tmp, v, ub); memcpy(v, v + u, ub); memcpy(v + u, tmp, ub); }
                        continue;
                    }
                    for (int i = 0; i < nb; i++) ord[i] = i;
                    bool sorted = true;
                    for (int i = 1; i < nb && sorted; i++) if (memcmp(v + (size_t)(i - 1) * u, v + (size_t)i * u, ub) > 0) sorted = false;
                    if (sorted) continue;
                    if (nb <= 32) { for (int i = 1; i < nb; i++) { int x = ord[i]; int j = i - 1; while (j >= 0 && memcmp(v + (size_t)ord[j] * u, v + (size_t)x * u, ub) > 0) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; } }
                    else std::sort(ord, ord + nb, [&](int A, int Bb) { return memcmp(v + (size_t)A * u, v + (size_t)Bb * u, ub) < 0; });
                    for (int i = 0; i < nb; i++) memcpy(tmp + (size_t)i * u, v + (size_t)ord[i] * u, ub);
                    memcpy(v, tmp, ub * nb);
                }
            }
            if (!have || memcmp(Sv, best, 2 * g) < 0) { memcpy(best, Sv, 2 * g); have = true; }
        }
        for (int p = 0; p < g; p++) out[outIdx[p]] = best[p];
    }

  private:
    void apply_affine(const uint16_t *c, uint16_t *out, CanonScratch &s) const {
        uint16_t *cc = s.cc.data();
        memcpy(cc, c, g * 2); memcpy(cc + g, c, g * 2);
        uint16_t mn = 0xffff;
        for (int x = 0; x < g; x++) mn = std::min(mn, c[x]);
        bool have = false;
        for (int t = 0; t < g; t++) {
            if (c[t] != mn) continue;
            const uint16_t *ct = cc + t;
            for (const auto &m : mul) {
                const uint16_t *mm = m.data();
                if (!have) { for (int x = 0; x < g; x++) out[x] = ct[mm[x]]; have = true; continue; }
                for (int x = 1; x < g; x++) {
                    uint16_t v = ct[mm[x]];
                    if (v < out[x]) { for (int y = x; y < g; y++) out[y] = ct[mm[y]]; break; }
                    if (v > out[x]) break;
                }
            }
        }
    }
};

static inline u64 hash_bytes(const uint8_t *k, int nb) {
    u64 h = 0x9E3779B97F4A7C15ULL ^ (u64)nb;
    int i = 0;
    for (; i + 8 <= nb; i += 8) {
        u64 w; memcpy(&w, k + i, 8);
        h ^= w; h *= 0xff51afd7ed558ccdULL; h ^= h >> 32;
    }
    if (i < nb) { u64 w = 0; memcpy(&w, k + i, nb - i); h ^= w; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 29; }
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
    return h;
}

// ---------------------------------------------------------------- state store
// Records: [weight NL*8 bytes][key kb bytes] padded to 8 bytes, stored in blocks.
template <int NL>
struct Arena {
    typedef UBig<NL> W;
    int kb = 1; size_t stride = 0, n = 0;
    static const size_t BLK = 1 << 10;
    std::vector<std::unique_ptr<char[]>> blocks;
    void init(int kb_) { kb = kb_; stride = (sizeof(W) + kb + 7) / 8 * 8; n = 0; blocks.clear(); }
    char *rec(size_t i) const { return blocks[i / BLK].get() + (i % BLK) * stride; }
    W *wt(size_t i) const { return (W *)rec(i); }
    uint8_t *key(size_t i) const { return (uint8_t *)(rec(i) + sizeof(W)); }
    size_t push() { if (n % BLK == 0) blocks.emplace_back(new char[BLK * stride]); return n++; }
    size_t bytes() const { return blocks.size() * BLK * stride; }
};

struct SpinLock {
    std::atomic_flag f = ATOMIC_FLAG_INIT;
    void lock() { while (f.test_and_set(std::memory_order_acquire)) {} }
    void unlock() { f.clear(std::memory_order_release); }
};

template <int NL>
struct Shard {
    typedef UBig<NL> W;
    SpinLock lk;
    Arena<NL> ar;
    std::vector<u64> tab;  // (tag<<32)|(idx+1), 0 = empty
    u64 mask = 0;
    char pad[64];
    void init(int kb) { ar.init(kb); tab.assign(16, 0); mask = 15; }
    void grow() {
        std::vector<u64> nt(tab.size() * 2, 0); u64 nm = nt.size() - 1;
        for (u64 e : tab) if (e) {
            size_t idx = (u32)e - 1;
            u64 h = hash_bytes(ar.key(idx), ar.kb);
            u64 p = (h >> 20) & nm;
            while (nt[p]) p = (p + 1) & nm;
            nt[p] = e;
        }
        tab.swap(nt); mask = nm;
    }
    void insert(const uint8_t *k, u64 h, const W &w) {
        int kb = ar.kb;
        u64 tag = h >> 32;
        lk.lock();
        u64 p = (h >> 20) & mask;
        while (true) {
            u64 e = tab[p];
            if (!e) break;
            if ((e >> 32) == tag) {
                size_t idx = (u32)e - 1;
                if (memcmp(ar.key(idx), k, kb) == 0) {
                    if (ar.wt(idx)->add(w)) die("weight overflow");
                    lk.unlock(); return;
                }
            }
            p = (p + 1) & mask;
        }
        size_t idx = ar.push();
        if (idx >= 0xfffffff0ULL) die("shard too large");
        *ar.wt(idx) = w; memcpy(ar.key(idx), k, kb);
        tab[p] = (tag << 32) | (u64)(idx + 1);
        if ((ar.n + 1) * 10 > tab.size() * 7) grow();
        lk.unlock();
    }
};

// ---------------------------------------------------------------- transitions
struct StepInfo {
    u64 q, d, g, g2, G, f;
    int kw_in, kw_out;               // bytes per entry of the stored keys
    std::vector<uint16_t> projG;     // X -> X mod g2 for X in Z_G
    std::vector<uint16_t> modg;      // X -> X mod g for X in Z_G
};
static void expand(const StepInfo &st, const uint16_t *c, uint32_t *e) {
    for (u64 X = 0; X < st.G; X++) {
        uint16_t v = c[st.modg[X]];
        if (v % st.f) die("non-exact expansion");
        e[X] = v / st.f;
    }
}
static inline void pack(const uint16_t *v, int g, int kw, uint8_t *out) {
    if (kw == 1) for (int i = 0; i < g; i++) out[i] = (uint8_t)v[i];
    else memcpy(out, v, 2 * g);
}
static inline void unpack(const uint8_t *k, int g, int kw, uint16_t *v) {
    if (kw == 1) for (int i = 0; i < g; i++) v[i] = k[i];
    else memcpy(v, k, 2 * g);
}

struct Options {
    int nthreads = 1, tail = 3, canon = 2;
    size_t maxstates = (size_t)-1;
    bool verbose = true;
};

struct StepLog { u64 q, d, g, g2; size_t states; double t; double rss; std::string canon; };

template <int NL>
static std::vector<WC> run_dp(const Problem &P, const Options &O, std::vector<StepLog> &logv) {
    typedef UBig<NL> W;
    int n = P.n;
    u64 T0 = P.T0;
    std::vector<WC> out(T0 + 1);
    if (n == 0) { out[T0] = WC(1); return out; }
    int tail = std::max(1, std::min(O.tail, n));
    int m = n - tail;
    std::vector<StepInfo> steps(n);
    for (int k = 0; k < n; k++) {
        StepInfo &s = steps[k];
        s.q = P.q[k]; s.d = P.d[k]; s.g = P.g[k]; s.g2 = P.g[k + 1];
        s.G = lcmu(s.g, s.d); s.f = s.G / s.g;
        if (s.G % s.g2) die("g2 does not divide G");
        if (s.G > 65535) die("G too large");
        s.projG.resize(s.G); s.modg.resize(s.G);
        for (u64 X = 0; X < s.G; X++) { s.projG[X] = (uint16_t)(X % s.g2); s.modg[X] = (uint16_t)(X % s.g); }
        s.kw_in = (T0 / s.g <= 255) ? 1 : 2;
        s.kw_out = (T0 / s.g2 <= 255) ? 1 : 2;
    }
    const int nsh = 1 << 12;
    typedef std::vector<Shard<NL>> Layer;
    Layer *cur = new Layer(1);
    (*cur)[0].init(steps[0].kw_in);
    {
        uint16_t k0 = (uint16_t)T0; uint8_t kb[2];
        pack(&k0, 1, steps[0].kw_in, kb);
        (*cur)[0].insert(kb, hash_bytes(kb, steps[0].kw_in), W(1));
    }
    double t0 = now();
    for (int k = 0; k < m; k++) {
        const StepInfo &st = steps[k];
        int g = st.g, g2 = st.g2, kwi = st.kw_in, kwo = st.kw_out;
        std::vector<u64> rem(P.d.begin() + k + 1, P.d.end());
        Canon cn; cn.init(g2, rem, O.canon);
        Layer *nxt = new Layer(nsh);
        for (auto &s : *nxt) s.init(g2 * kwo);
        std::vector<std::pair<int, size_t>> items;
        size_t nin = 0;
        for (size_t si = 0; si < cur->size(); si++) {
            auto &a = (*cur)[si].ar; nin += a.n;
            for (size_t bI = 0; bI * Arena<NL>::BLK < a.n; bI++) items.push_back({(int)si, bI});
        }
        std::atomic<size_t> next(0);
        auto worker = [&]() {
            std::vector<uint32_t> e(st.G), base(g2);
            std::vector<uint16_t> c(g), child(g2), outv(g2);
            std::vector<uint8_t> kb(g2 * kwo + 8);
            CanonScratch cs; cn.prep(cs);
            while (true) {
                size_t it = next.fetch_add(1);
                if (it >= items.size()) break;
                auto &a = (*cur)[items[it].first].ar;
                size_t lo = items[it].second * Arena<NL>::BLK, hi = std::min(a.n, lo + Arena<NL>::BLK);
                for (size_t i = lo; i < hi; i++) {
                    unpack(a.key(i), g, kwi, c.data());
                    const W &w = *a.wt(i);
                    expand(st, c.data(), e.data());
                    std::fill(base.begin(), base.end(), 0);
                    for (u64 X = 0; X < st.G; X++) base[st.projG[X]] += e[X];
                    auto emit = [&](const W &ww) {
                        cn.apply(child.data(), outv.data(), cs);
                        pack(outv.data(), g2, kwo, kb.data());
                        u64 h = hash_bytes(kb.data(), g2 * kwo);
                        (*nxt)[h & (nsh - 1)].insert(kb.data(), h, ww);
                    };
                    for (int y = 0; y < g2; y++) child[y] = (uint16_t)base[y];
                    W wn = w; if (wn.mul_small(st.q - st.d)) die("overflow");
                    emit(wn);
                    for (u64 s = 0; s < st.d; s++) {
                        for (int y = 0; y < g2; y++) child[y] = (uint16_t)base[y];
                        for (u64 X = s; X < st.G; X += st.d) child[st.projG[X]] -= (uint16_t)e[X];
                        emit(w);
                    }
                }
                a.blocks[items[it].second].reset();  // parent block consumed: release memory early
            }
        };
        std::vector<std::thread> th;
        for (int t = 0; t < O.nthreads; t++) th.emplace_back(worker);
        for (auto &t : th) t.join();
        delete cur;
        cur = nxt;
        size_t nout = 0, bytes = 0;
        for (auto &s : *cur) { nout += s.ar.n; bytes += s.ar.bytes() + s.tab.size() * 8; std::vector<u64>().swap(s.tab); }
        logv.push_back({st.q, st.d, (u64)g, (u64)g2, nout, now() - t0, peak_rss_gb(), cn.desc});
        if (O.verbose)
            fprintf(stderr, "step %2d q=%-6llu d=%-4llu g=%-4d -> g'=%-4d states %zu -> %zu  (%.0f MB) t=%.1fs peakrss=%.2fGB  sym %s\n", k,
                    (unsigned long long)st.q, (unsigned long long)st.d, g, g2, nin, nout, bytes / 1e6, now() - t0, peak_rss_gb(), cn.desc.c_str());
        if (nout > O.maxstates) { fprintf(stderr, "ABORT: layer exceeds --maxstates\n"); exit(3); }
    }
    // tail: DFS over steps m..n-1 for every stored state
    {
        {
            u128 pq = 1;
            for (int k = m; k < n; k++) { pq *= steps[k].q; if (pq >= ((u128)1 << 63)) die("tail too long for 64-bit local counts; reduce --tail"); }
        }
        std::vector<std::pair<int, size_t>> items;
        size_t nin = 0;
        for (size_t si = 0; si < cur->size(); si++) {
            auto &a = (*cur)[si].ar; nin += a.n;
            for (size_t bI = 0; bI * Arena<NL>::BLK < a.n; bI++) items.push_back({(int)si, bI});
        }
        std::atomic<size_t> next(0);
        std::mutex hm;
        std::vector<W> R(T0 + 1);
        auto worker = [&]() {
            // Per stored state, the tail is enumerated with small-integer multipliers
            // (products of q - d over "none" choices), accumulated into a local u64
            // histogram (total mass per state = prod of tail primes < 2^63, checked),
            // then flushed once per touched bin as W * count.
            std::vector<W> hist(T0 + 1);
            std::vector<u64> cnt(T0 + 1, 0);
            std::vector<u32> touched; touched.reserve(T0 + 1);
            int L = n - m;
            std::vector<std::vector<uint32_t>> E(L), B(L);
            std::vector<std::vector<uint16_t>> C(L + 1);
            C[0].resize(steps[m].g);
            for (int j = 0; j < L; j++) { E[j].resize(steps[m + j].G); B[j].resize(steps[m + j].g2); C[j + 1].resize(steps[m + j].g2); }
            auto bump = [&](u64 v, u64 add) { if (!cnt[v]) touched.push_back((u32)v); cnt[v] += add; };
            // Fused last two primes A, B (orders dA, dB): with Q = lcm(dA, dB) (a multiple
            // of the current cut modulus g), expand c to e on Z_Q; then for classes sA, sB
            // nu = tot - mA[sA] - mB[sB] + (sA = sB mod gcd ? e[CRT(sA, sB)] : 0).
            const bool fused = (L >= 2);
            u64 dA = 0, dB = 0, Q = 0, qA = 0, qB = 0, gF = 0;
            std::vector<int> crt;
            std::vector<u64> eQ, mA, mB;
            if (fused) {
                const StepInfo &sa = steps[n - 2], &sb = steps[n - 1];
                dA = sa.d; dB = sb.d; qA = sa.q; qB = sb.q; Q = lcmu(dA, dB); gF = sa.g;
                if (Q % gF) die("fused tail: g does not divide lcm");
                crt.assign(dA * dB, -1);
                for (u64 X = 0; X < Q; X++) crt[(X % dA) * dB + (X % dB)] = (int)X;
                eQ.resize(Q); mA.resize(dA); mB.resize(dB);
            }
            auto fused2 = [&](const uint16_t *c, u64 mu) {
                const u64 f = Q / gF;
                u64 tot = 0;
                for (u64 X = 0; X < Q; X++) {
                    uint16_t v = c[X % gF];
                    if (v % f) die("non-exact expansion (fused)");
                    eQ[X] = v / f; tot += eQ[X];
                }
                std::fill(mA.begin(), mA.end(), 0); std::fill(mB.begin(), mB.end(), 0);
                for (u64 X = 0; X < Q; X++) { mA[X % dA] += eQ[X]; mB[X % dB] += eQ[X]; }
                const u64 wA = qA - dA, wB = qB - dB;
                bump(tot, mu * wA * wB);
                for (u64 t = 0; t < dB; t++) bump(tot - mB[t], mu * wA);
                for (u64 sA = 0; sA < dA; sA++) {
                    u64 base = tot - mA[sA];
                    bump(base, mu * wB);
                    const int *cr = &crt[sA * dB];
                    for (u64 t = 0; t < dB; t++) {
                        u64 v = base - mB[t];
                        if (cr[t] >= 0) v += eQ[cr[t]];
                        bump(v, mu);
                    }
                }
            };
            std::function<void(int, const uint16_t *, u64)> rec = [&](int j, const uint16_t *c, u64 mu) {
                if (fused && j == L - 2) { fused2(c, mu); return; }
                const StepInfo &st = steps[m + j];
                uint32_t *e = E[j].data();
                expand(st, c, e);
                if (j == L - 1) {  // last prime: next modulus is 1, G == d
                    u64 tot = 0;
                    for (u64 X = 0; X < st.G; X++) tot += e[X];
                    bump(tot, mu * (st.q - st.d));
                    for (u64 s = 0; s < st.d; s++) {
                        u64 rm = 0;
                        for (u64 X = s; X < st.G; X += st.d) rm += e[X];
                        bump(tot - rm, mu);
                    }
                    return;
                }
                uint32_t *bb = B[j].data();
                int g2 = (int)st.g2;
                std::fill(bb, bb + g2, 0);
                for (u64 X = 0; X < st.G; X++) bb[st.projG[X]] += e[X];
                uint16_t *ch = C[j + 1].data();
                for (int y = 0; y < g2; y++) ch[y] = (uint16_t)bb[y];
                rec(j + 1, ch, mu * (st.q - st.d));
                for (u64 s = 0; s < st.d; s++) {
                    for (int y = 0; y < g2; y++) ch[y] = (uint16_t)bb[y];
                    for (u64 X = s; X < st.G; X += st.d) ch[st.projG[X]] -= (uint16_t)e[X];
                    rec(j + 1, ch, mu);
                }
            };
            while (true) {
                size_t it = next.fetch_add(1);
                if (it >= items.size()) break;
                if (O.verbose && items.size() >= 200 && it % (items.size() / 20) == 0)
                    fprintf(stderr, "  tail progress %3.0f%%  t=%.1fs\n", 100.0 * it / items.size(), now() - t0);
                auto &a = (*cur)[items[it].first].ar;
                size_t lo = items[it].second * Arena<NL>::BLK, hi = std::min(a.n, lo + Arena<NL>::BLK);
                for (size_t i = lo; i < hi; i++) {
                    unpack(a.key(i), steps[m].g, steps[m].kw_in, C[0].data());
                    rec(0, C[0].data(), 1);
                    const W &w = *a.wt(i);
                    for (u32 v : touched) {
                        W t = w;
                        if (t.mul_small(cnt[v])) die("overflow");
                        if (hist[v].add(t)) die("overflow");
                        cnt[v] = 0;
                    }
                    touched.clear();
                }
            }
            std::lock_guard<std::mutex> lg(hm);
            for (u64 v = 0; v <= T0; v++) if (R[v].add(hist[v])) die("overflow");
        };
        std::vector<std::thread> th;
        for (int t = 0; t < O.nthreads; t++) th.emplace_back(worker);
        for (auto &t : th) t.join();
        if (O.verbose) fprintf(stderr, "tail over %d primes from %zu states done t=%.1fs peakrss=%.2fGB\n", n - m, nin, now() - t0, peak_rss_gb());
        for (u64 v = 0; v <= T0; v++) out[v] = R[v].template widen<6>();
    }
    delete cur;
    return out;
}

// ---------------------------------------------------------------- brute force
// Literal definition: for every a in [0, M0) compute nu(a) = #{r in Z_T0 :
// a mod q != 2^r mod q for all q}.  For each prime and each residue a mod q the
// set {r in Z_T0 : 2^r = a mod q} is precomputed by testing every r in Z_T0.
static std::vector<u64> brute_force(const std::vector<u64> &qs, int nthreads) {
    u64 T0 = 1, M0 = 1;
    for (u64 q : qs) { T0 = lcmu(T0, ord2(q)); M0 *= q; if (M0 > (1ULL << 40)) die("brute: M0 too large"); }
    int nw = (T0 + 63) / 64;
    int np = qs.size();
    std::vector<std::vector<u64>> mask(np);
    for (int p = 0; p < np; p++) {
        u64 q = qs[p];
        mask[p].assign(q * nw, 0);
        u64 pw = 1 % q;
        for (u64 r = 0; r < T0; r++) { mask[p][pw * nw + r / 64] |= 1ULL << (r % 64); pw = pw * 2 % q; }
    }
    std::vector<std::vector<u64>> hs(nthreads, std::vector<u64>(T0 + 1, 0));
    std::atomic<u64> next(0);
    const u64 CH = 1 << 16;
    auto worker = [&](int tid) {
        std::vector<u64> acc(nw);
        while (true) {
            u64 lo = next.fetch_add(CH);
            if (lo >= M0) break;
            u64 hi = std::min(M0, lo + CH);
            for (u64 a = lo; a < hi; a++) {
                std::fill(acc.begin(), acc.end(), 0);
                for (int p = 0; p < np; p++) {
                    const u64 *mk = &mask[p][(a % qs[p]) * nw];
                    for (int w = 0; w < nw; w++) acc[w] |= mk[w];
                }
                u64 removed = 0;
                for (int w = 0; w < nw; w++) removed += __builtin_popcountll(acc[w]);
                hs[tid][T0 - removed]++;
            }
        }
    };
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; t++) th.emplace_back(worker, t);
    for (auto &t : th) t.join();
    std::vector<u64> h(T0 + 1, 0);
    for (auto &x : hs) for (u64 v = 0; v <= T0; v++) h[v] += x[v];
    return h;
}

// ---------------------------------------------------------------- moment checks
// sum_a nu(a)^j = sum_{r_1..r_j in Z_T0} #{a : all r_i survive}
//              = sum_{r_1..r_j} prod_q (q - #distinct{r_i mod d_q}),
// and by translation invariance = T0 * (same sum with r_1 = 0).
static bool moment_checks(const std::vector<u64> &qs, const std::vector<u64> &ds, u64 T0, const std::vector<WC> &H, int nthreads) {
    int np = qs.size();
    WC m0, m1, m2, m3;
    for (u64 v = 0; v <= T0; v++) {
        WC h = H[v];
        bool of = false;
        of |= m0.add(h) != 0;
        WC t = h; of |= t.mul_small(v) != 0; of |= m1.add(t) != 0;
        of |= t.mul_small(v) != 0; of |= m2.add(t) != 0;
        of |= t.mul_small(v) != 0; of |= m3.add(t) != 0;
        if (of) die("moment overflow");
    }
    WC M0(1), e1(T0), e2, e3;
    for (u64 q : qs) M0.mul_small(q);
    for (int p = 0; p < np; p++) e1.mul_small(qs[p] - 1);
    for (u64 x = 0; x < T0; x++) {
        WC t(1);
        for (int p = 0; p < np; p++) t.mul_small(qs[p] - 1 - (x % ds[p] != 0));
        e2.add(t);
    }
    e2.mul_small(T0);
    std::vector<WC> part(nthreads);
    std::atomic<u64> next(0);
    auto worker = [&](int tid) {
        while (true) {
            u64 x = next.fetch_add(1);
            if (x >= T0) break;
            for (u64 y = 0; y < T0; y++) {
                WC t(1);
                for (int p = 0; p < np; p++) {
                    u64 dd = ds[p], a = x % dd, b = y % dd;
                    u64 cnt = 1 + (a != 0) + (b != 0 && b != a);
                    t.mul_small(qs[p] - cnt);
                }
                part[tid].add(t);
            }
        }
    };
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; t++) th.emplace_back(worker, t);
    for (auto &t : th) t.join();
    for (auto &p : part) e3.add(p);
    e3.mul_small(T0);
    bool o0 = m0 == M0, o1 = m1 == e1, o2 = m2 == e2, o3 = m3 == e3;
    fprintf(stderr, "check sum delta        = M0          : %s (M0 = %s)\n", o0 ? "OK" : "FAIL", M0.str().c_str());
    fprintf(stderr, "check sum delta*nu     = T0*phi(M0)  : %s\n", o1 ? "OK" : "FAIL");
    fprintf(stderr, "check sum delta*nu^2   (direct form) : %s\n", o2 ? "OK" : "FAIL");
    fprintf(stderr, "check sum delta*nu^3   (direct form) : %s\n", o3 ? "OK" : "FAIL");
    return o0 && o1 && o2 && o3;
}

// ---------------------------------------------------------------- driver
static Problem make_problem(const std::vector<u64> &qs_in, const std::string &order_mode) {
    for (u64 q : qs_in)
        if (q < 3 || q % 2 == 0 || !isprime(q)) { fprintf(stderr, "not an odd prime: %llu\n", (unsigned long long)q); exit(2); }
    { auto s = qs_in; std::sort(s.begin(), s.end()); if (std::adjacent_find(s.begin(), s.end()) != s.end()) die("duplicate prime"); }
    std::vector<u64> d0;
    for (u64 q : qs_in) d0.push_back(ord2(q));
    std::vector<int> ord;
    if (order_mode == "given") { ord.resize(qs_in.size()); std::iota(ord.begin(), ord.end(), 0); }
    else if (order_mode == "auto") ord = auto_order(d0);
    else {
        std::vector<u64> lst; size_t p = 0;
        while (p < order_mode.size()) { size_t c = order_mode.find(',', p); if (c == std::string::npos) c = order_mode.size(); lst.push_back(std::stoull(order_mode.substr(p, c - p))); p = c + 1; }
        for (u64 q : lst) { auto it = std::find(qs_in.begin(), qs_in.end(), q); if (it == qs_in.end()) die("order list prime not in set"); ord.push_back(it - qs_in.begin()); }
        auto so = ord; std::sort(so.begin(), so.end());
        if (so.size() != qs_in.size() || std::adjacent_find(so.begin(), so.end()) != so.end()) die("order list is not a permutation");
    }
    Problem P; P.n = qs_in.size();
    for (int i : ord) { P.q.push_back(qs_in[i]); P.d.push_back(d0[i]); }
    for (u64 d : P.d) P.T0 = lcmu(P.T0, d);
    if (P.T0 > 65535) die("T0 too large for 16-bit state entries");
    P.g = compute_g(P.d);
    return P;
}

static int limbs_needed(const std::vector<u64> &qs) {
    UBig<8> m(1);
    for (u64 q : qs) if (m.mul_small(q)) die("M0 too large");
    int top = 0; for (int i = 0; i < 8; i++) if (m.l[i]) top = i;
    return top + 1;  // weights never exceed M0
}

static int run_one(const std::vector<u64> &qs, const std::string &order_mode, const Options &O, bool do_brute, bool print, bool plan_only) {
    Problem P = make_problem(qs, order_mode);
    int nl = limbs_needed(qs);
    if (O.verbose || plan_only) {
        fprintf(stderr, "primes (processing order): %s\n", join(P.q).c_str());
        fprintf(stderr, "orders                   : %s\n", join(P.d).c_str());
        fprintf(stderr, "T0 = %llu, cut moduli g_k: %s, weight limbs %d\n", (unsigned long long)P.T0, join(P.g).c_str(), nl);
    }
    if (plan_only) return 0;
    double t0 = now();
    std::vector<StepLog> lg;
    std::vector<WC> H;
    if (nl <= 2) H = run_dp<2>(P, O, lg);
    else if (nl == 3) H = run_dp<3>(P, O, lg);
    else if (nl == 4) H = run_dp<4>(P, O, lg);
    else die("M0 needs more than 256 bits");
    double t1 = now();
    bool ok = moment_checks(P.q, P.d, P.T0, H, O.nthreads);
    size_t ne = 0;
    for (u64 v = 0; v <= P.T0; v++) if (!H[v].iszero()) ne++;
    size_t mx = 0; for (auto &l : lg) mx = std::max(mx, l.states);
    fprintf(stderr, "histogram entries: %zu, max stored layer %zu states, DP time %.2fs, total %.2fs, peak RSS %.3f GB\n", ne, mx, t1 - t0, now() - t0, peak_rss_gb());
    if (do_brute) {
        std::vector<u64> b = brute_force(qs, O.nthreads);
        bool same = true;
        for (u64 v = 0; v <= P.T0; v++) if (!(H[v] == WC(b[v]))) same = false;
        fprintf(stderr, "brute force comparison: %s\n", same ? "MATCH" : "MISMATCH");
        ok = ok && same;
    }
    if (print) {
        std::string out;
        for (u64 v = 0; v <= P.T0; v++) if (!H[v].iszero()) { out += std::to_string(v); out += ' '; out += H[v].str(); out += '\n'; }
        fwrite(out.data(), 1, out.size(), stdout);
        fflush(stdout);
    }
    if (!ok) fprintf(stderr, "*** CHECK FAILED ***\n");
    return ok ? 0 : 1;
}

int main(int argc, char **argv) {
    Options O;
    O.nthreads = std::max(1u, std::thread::hardware_concurrency());
    std::string order_mode = "auto";
    bool brute = false, bruteonly = false, selftest = false, plan = false;
    std::vector<u64> qs;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a.rfind("--threads=", 0) == 0) O.nthreads = atoi(a.c_str() + 10);
        else if (a.rfind("--order=", 0) == 0) order_mode = a.substr(8);
        else if (a.rfind("--tail=", 0) == 0) O.tail = atoi(a.c_str() + 7);
        else if (a.rfind("--maxstates=", 0) == 0) O.maxstates = std::stoull(a.substr(12));
        else if (a.rfind("--canon=", 0) == 0) { std::string m = a.substr(8); O.canon = m == "none" ? 0 : m == "affine" ? 1 : m == "tree" ? 2 : -1; if (O.canon < 0) die("bad --canon"); }
        else if (a == "--brute") brute = true;
        else if (a == "--bruteonly") bruteonly = true;
        else if (a == "--selftest") selftest = true;
        else if (a == "--plan") plan = true;
        else if (a == "--quiet") O.verbose = false;
        else if (a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
        else {
            for (char &ch : a) if (ch == ',' || ch == ' ' || ch == '\t' || ch == '[' || ch == ']') ch = ' ';
            size_t p = 0;
            while (p < a.size()) { size_t c = a.find(' ', p); if (c == std::string::npos) c = a.size(); if (c > p) qs.push_back(std::stoull(a.substr(p, c - p))); p = c + 1; }
        }
    }
    if (O.nthreads < 1) O.nthreads = 1;
    if (selftest) {
        std::vector<std::vector<u64>> sets = {
            {3}, {7}, {3, 5}, {3, 5, 7, 13}, {3, 5, 7, 13, 17}, {3, 7, 5, 13, 17, 11, 43}, {5, 31, 11, 3, 41},
            {3, 5, 7, 11, 13, 17}, {73, 19, 7, 3, 37}, {257, 17, 5, 3, 7}, {3, 5, 7, 11, 13, 17, 31}, {127, 43, 29, 3, 5},
            {73, 19, 37, 3, 5}, {19, 37, 73, 13, 7}, {109, 37, 19, 5}, {29, 113, 43, 127, 3}, {241, 97, 17, 7}, {257, 97, 3, 5}};
        int bad = 0, runs = 0;
        Options Q = O; Q.verbose = false;
        for (auto &s : sets) {
            std::vector<u64> b = brute_force(s, O.nthreads);
            for (std::string om : {"auto", "given"})
                for (int tl : {1, 2, 3})
                    for (int cn : {2, 1, 0}) {
                        Q.tail = tl; Q.canon = cn;
                        Problem P = make_problem(s, om);
                        std::vector<StepLog> lg;
                        std::vector<WC> H = run_dp<2>(P, Q, lg);
                        bool same = true;
                        for (u64 v = 0; v <= P.T0; v++) if (!(H[v] == WC(b[v]))) same = false;
                        runs++;
                        if (!same) { bad++; fprintf(stderr, "MISMATCH: set %s order=%s tail=%d canon=%d\n", join(s).c_str(), om.c_str(), tl, cn); }
                    }
            fprintf(stderr, "selftest set %-28s compared with brute-force enumeration of a mod M0 in 18 configurations\n", join(s).c_str());
        }
        fprintf(stderr, "SELFTEST %s (%d runs, %d failures)\n", bad ? "FAILED" : "PASSED", runs, bad);
        return bad ? 1 : 0;
    }
    if (qs.empty()) { fprintf(stderr, "usage: histB [options] q1 q2 ...\n"); return 2; }
    if (bruteonly) {
        std::vector<u64> b = brute_force(qs, O.nthreads);
        for (size_t v = 0; v < b.size(); v++) if (b[v]) printf("%zu %llu\n", v, (unsigned long long)b[v]);
        return 0;
    }
    return run_one(qs, order_mode, O, brute, true, plan);
}
