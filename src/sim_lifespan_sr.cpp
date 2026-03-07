#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Xoroshiro128+ fast PRNG (not cryptographic, fine for simulation)
struct Xoro128 {
  uint64_t s[2];

  void seed(uint64_t seed) {
    // SplitMix64 to seed from single value
    uint64_t z = seed;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    s[0] = z ^ (z >> 31);
    z = (seed + 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    s[1] = z ^ (z >> 31);
  }

  uint64_t next() {
    uint64_t s0 = s[0], s1 = s[1];
    uint64_t result = s0 + s1;
    s1 ^= s0;
    s[0] = ((s0 << 24) | (s0 >> 40)) ^ s1 ^ (s1 << 16);
    s[1] = (s1 << 37) | (s1 >> 27);
    return result;
  }

  // Uniform (0, 1)
  double runif() {
    return (next() >> 11) * 0x1.0p-53;
  }

  // Box-Muller normal
  double rnorm() {
    double u1 = runif(), u2 = runif();
    if (u1 < 1e-300) u1 = 1e-300;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
  }

  // Exponential(rate)
  double rexp(double rate) {
    double u = runif();
    if (u < 1e-300) u = 1e-300;
    return -std::log(u) / rate;
  }
};

// [[Rcpp::export]]
NumericVector sim_lifespan_sr_cpp(NumericVector Xc, NumericVector m_ex,
                                   double eta, double beta,
                                   double epsilon, double kappa,
                                   double dt, double t_max,
                                   int rng_seed = 0) {
  int n = Xc.size();
  NumericVector death_age(n, t_max);
  std::vector<double> x(n, 0.0);

  double noise_sd = std::sqrt(2.0 * epsilon * dt);
  int n_steps = (int)std::ceil(t_max / dt);

  // Use fast internal RNG if seed provided, otherwise fall back to R's RNG
  bool use_fast_rng = (rng_seed != 0);
  Xoro128 rng;
  if (use_fast_rng) {
    rng.seed((uint64_t)rng_seed);
  }

  // Compact alive index array
  std::vector<int> alive_idx(n);
  for (int i = 0; i < n; i++) alive_idx[i] = i;
  int n_alive = n;

  for (int k = 0; k < n_steps; k++) {
    if (n_alive == 0) break;
    double t_now = k * dt;

    int write = 0;
    for (int r = 0; r < n_alive; r++) {
      int i = alive_idx[r];

      double rate = m_ex[i];
      if (rate < 1e-20) rate = 1e-20;
      double t_ext = use_fast_rng ? rng.rexp(rate) : R::rexp(1.0 / rate);

      double x_old = x[i];
      double drift = (eta * t_now - beta * x_old / (x_old + kappa)) * dt;
      double noise = use_fast_rng ? rng.rnorm() : R::rnorm(0.0, 1.0);
      double x_new = x_old + drift + noise_sd * noise;
      if (x_new < 0.0) x_new = 0.0;

      double t_int = dt + 1.0;
      if (x_new > Xc[i]) {
        double denom = x_new - x_old;
        if (std::abs(denom) < 1e-12) denom = 1e-12;
        double frac = (Xc[i] - x_old) / denom;
        if (frac < 0.0) frac = 0.0;
        if (frac > 1.0) frac = 1.0;
        t_int = frac * dt;
      }

      double t_death = (t_ext < t_int) ? t_ext : t_int;
      if (t_death < dt) {
        death_age[i] = t_now + t_death;
      } else {
        x[i] = x_new;
        alive_idx[write] = i;
        write++;
      }
    }
    n_alive = write;
  }

  return death_age;
}
