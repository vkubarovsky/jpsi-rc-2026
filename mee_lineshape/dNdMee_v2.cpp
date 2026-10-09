// dNdMee_v2.cpp -- the M(e+e-) lineshape of V -> l+ l- (gamma), QED final-state
// radiation, optionally folded with a Gaussian detector resolution.
//
// Successor of VPK's dNdMee.cpp, same call  dNdMee(mee, M, sigma, br).  Three
// things are corrected (see mee_lineshape/report/Mee_lineshape.pdf):
//
//  1. The photon energy in the V rest frame is k = (M^2 - mee^2) / (2M), not
//     x = M - mee (they differ by up to 20 % at mee = 2 GeV), and a density in
//     mee carries the Jacobian dk/dmee = mee / M.
//  2. The spectrum is exponentiated.  The first-order Ehlotzky-Mitter spectrum
//     dNdk0 ~ 1/k is integrable only with a cut eps and puts too many photons
//     just above it.  Here the cumulative distribution is
//         F(K) = (K/E)^beta [1 + delta(K) - beta ln(K/E)] / [1 + delta(E)],
//     beta = (2 alpha/pi)(ln(M^2/m^2) - 1), delta(K) the first-order correction
//     with photons below K (Ehlotzky-Mitter, rmax(rmax-6)/16 term included,
//     ~/jpsi-rc-2026/rc_vector_mesons.py:delta_vm).  F(E) = 1 exactly: the
//     lineshape is a probability density, normalised by the measured
//     (radiation-inclusive) branching ratio, not by Gamma_0.
//  3. br is not a free parameter.  The fraction of decays with a photon above
//     eps is 1 - F(eps), fixed by QED; the default br < 0 uses it.  Passing
//     br in [0,1] overrides it (for studies only).
//
// The peak (k < eps) is a delta at mee = M of weight F(eps).  Keep eps well
// below the resolution sigma; the folded shape is then independent of eps
// (verified for eps = 0.1 ... 10 MeV against sigma = 30 MeV).
//
// Self-contained, no boost.  Build as a program (main below) or as a shared
// library:  g++ -O2 -shared -fPIC -DNO_MAIN dNdMee_v2.cpp -o libdNdMee.so
#include <cmath>
#include <cstdio>
#include <initializer_list>

namespace lineshape {

const double PI = 3.141592653589793;
const double ALPHA = 1.0 / 137.035999;
const double MASS_E = 0.51099895e-3;   // electron mass (M_E is a math.h macro)

inline double li2(double x) {
    // real dilogarithm Li2(x), 0 <= x <= 1 (series + reflection)
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return PI * PI / 6.0;           // avoid ln(1) ln(0) = 0 * inf
    if (x > 0.5) return PI * PI / 6.0 - std::log(x) * std::log(1.0 - x) - li2(1.0 - x);
    double s = 0.0, t = x;
    for (int n = 1; n < 200; ++n) { s += t / (double(n) * n); t *= x; if (t < 1e-18) break; }
    return s;
}

// delta(K): first-order correction to Gamma(V -> l l) with photons k < K,
// relative to the lowest-order width (Ehlotzky-Mitter; delta_vm in
// rc_vector_mesons.py).  Valid for 0 < K <= M/2.
inline double delta_vm(double K, double M, double m) {
    double E = M / 2.0;
    double r = K / E; if (r > 1.0) r = 1.0; if (r < 1e-300) r = 1e-300;
    double L = std::log(2.0 * E / m);
    double lsq = 0.5 * std::log(std::fmax(1.0 - r, 1e-300));
    double R = (r - r * r / 4.0 - 0.75) * (L + lsq - 0.5);
    double dd = r * (r - 6.0) / 16.0;
    return -4.0 * ALPHA / PI * (dd + (0.5 - L) * (std::log(r) + 1.0 / 3.0)
                                + 0.5 * (li2(r) - PI * PI / 6.0) + 17.0 / 72.0 + R);
}

// first-order photon spectrum (1/Gamma0) dGamma/dk (Ehlotzky-Mitter) = d delta/dK
inline double dNdk0(double k, double M, double m) {
    double E = M / 2.0, x = k / E;
    if (k <= 0.0 || x >= 1.0) return 0.0;
    return (2.0 * ALPHA / PI) * (2.0 * std::log(2.0 * E / m * std::sqrt(1.0 - x)) - 1.0)
           / k * (1.0 - x + x * x / 2.0);
}

inline double beta_ir(double M, double m) {
    return 2.0 * ALPHA / PI * (std::log(M * M / (m * m)) - 1.0);
}

// exponentiated cumulative distribution of the photon energy, F(E) = 1
inline double F_cdf(double K, double M, double m) {
    double E = M / 2.0;
    if (K <= 0.0) return 0.0;
    if (K >= E) return 1.0;
    double b = beta_ir(M, m);
    return std::pow(K / E, b) * (1.0 + delta_vm(K, M, m) - b * std::log(K / E))
           / (1.0 + delta_vm(E, M, m));
}

// its density dF/dK = (K/E)^beta [dNdk0(K) + beta (delta(K) - beta ln(K/E)) / K] / (1+delta(E))
inline double dNdk_exp(double k, double M, double m) {
    double E = M / 2.0;
    if (k <= 0.0 || k >= E) return 0.0;
    double b = beta_ir(M, m);
    double lk = std::log(k / E);
    return std::pow(k / E, b) * (dNdk0(k, M, m) + b * (delta_vm(k, M, m) - b * lk) / k)
           / (1.0 + delta_vm(E, M, m));
}

inline double gauss(double z, double s) {
    return std::exp(-0.5 * z * z / (s * s)) / (std::sqrt(2.0 * PI) * s);
}

// tail density in mee, for photons above eps, NOT folded
inline double tail_true(double mee, double M, double m) {
    double k = (M * M - mee * mee) / (2.0 * M);
    return dNdk_exp(k, M, m) * mee / M;
}

}  // namespace lineshape

// dN/dmee, normalised to 1 over the full range (unit: 1/GeV).
//   mee, M, sigma in GeV;  sigma <= 0 gives the unfolded tail (the delta peak
//   at M is then omitted -- use peak_weight() for it);
//   br < 0: QED value 1 - F(eps);  eps default 1 MeV;  m = lepton mass.
extern "C" double dNdMee_full(double mee, double M, double sigma, double br,
                              double eps, double m) {
    using namespace lineshape;
    double kmax = (M * M - 4.0 * m * m) / (2.0 * M);
    double Feps = F_cdf(eps, M, m);
    double brq = 1.0 - Feps;
    double scale = (br < 0.0) ? 1.0 : br / brq;      // tail scaled if br is forced
    double pk = (br < 0.0) ? Feps : 1.0 - br;
    if (sigma <= 0.0) {
        double k = (M * M - mee * mee) / (2.0 * M);
        if (k < eps || k > kmax) return 0.0;
        return scale * tail_true(mee, M, m);
    }
    double res = pk * gauss(mee - M, sigma);
    // tail folded with the Gaussian: integrate over k in ln k (the tail ~ 1/k),
    // restricted to |mee - mee(k)| < 8 sigma
    double mlo = mee - 8.0 * sigma, mhi = mee + 8.0 * sigma;
    double klo = eps, khi = kmax;
    if (mhi < M) { double kk = (M * M - mhi * mhi) / (2.0 * M); if (kk > klo) klo = kk; }
    if (mlo > 2.0 * m) { double kk = (M * M - mlo * mlo) / (2.0 * M); if (kk < khi) khi = kk; }
    if (klo >= khi) return res;
    const int N = 400;                                  // Simpson, even N
    double a = std::log(klo), b = std::log(khi), h = (b - a) / N, s = 0.0;
    for (int i = 0; i <= N; ++i) {
        double k = std::exp(a + i * h);
        double mk = std::sqrt(std::fmax(M * M - 2.0 * M * k, 0.0));
        double f = dNdk_exp(k, M, m) * k * gauss(mee - mk, sigma);   // dk = k dlnk
        s += f * ((i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0));
    }
    return res + scale * s * h / 3.0;
}

// weight of the unresolved peak (photon below eps)
extern "C" double peak_weight(double M, double eps, double m) {
    return lineshape::F_cdf(eps, M, m);
}

// VPK's interface: dNdMee(mee, M, sigma, br) with eps = 1 MeV, electrons;
// br < 0 (recommended) takes the QED value.
extern "C" double dNdMee(double mee, double M, double sigma, double br) {
    return dNdMee_full(mee, M, sigma, br, 1.0e-3, lineshape::MASS_E);
}

#ifndef NO_MAIN
int main() {
    const double M = 3.0969;
    std::printf("beta = %.5f   1 - F(1 MeV) = %.5f   P(mee < 3.0) = %.5f\n",
                lineshape::beta_ir(M, lineshape::MASS_E),
                1.0 - lineshape::F_cdf(1e-3, M, lineshape::MASS_E),
                1.0 - lineshape::F_cdf((M * M - 9.0) / (2 * M), M, lineshape::MASS_E));
    // normalisation check with sigma = 30 MeV: integral over mee must be 1
    double s = 0.0, h = 0.0005;
    for (double x = 1.0; x < 3.4; x += h) s += dNdMee(x, M, 0.030, -1.0) * h;
    std::printf("integral of dNdMee(sigma = 30 MeV) over 1.0-3.4 GeV = %.5f\n", s);
    for (double mee : {2.5, 2.8, 2.9, 3.0, 3.05, 3.08, 3.097, 3.12})
        std::printf("  dNdMee(%.3f, M, 30 MeV) = %.5f /GeV\n", mee, dNdMee(mee, M, 0.030, -1.0));
    return 0;
}
#endif
