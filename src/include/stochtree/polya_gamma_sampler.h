/*!
 * Pólya-Gamma random variates PG(b, z), b > 0, z real.
 *
 * Written for the nashtree fork from the published algorithms only; no
 * third-party implementation was consulted or copied. References:
 *
 *  [D09]  Devroye, L. (2009). On exact simulation algorithms for some
 *         distributions related to Jacobi theta functions. Statistics &
 *         Probability Letters 79(21), 2251-2259.
 *  [PSW13] Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian
 *         inference for logistic models using Pólya-Gamma latent variables.
 *         Journal of the American Statistical Association 108(504),
 *         1339-1349. Definition (infinite convolution of gammas), the moment
 *         formulas, and Section 4 (Devroye's alternating-series sampler for
 *         J*(1, c), with PG(1, z) = J*(1, z / 2) / 4).
 *  [WPS14] Windle, J., Polson, N. G. and Scott, J. G. (2014). Sampling
 *         Pólya-Gamma random variates: alternate and approximate techniques.
 *         arXiv:1405.0506. (Truncated gamma-series and normal approximations
 *         for general b.)
 *  [MSH76] Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating
 *         random variates using transformations with multiple roots. The
 *         American Statistician 30(2), 88-90. (Inverse-Gaussian variates.)
 *
 * Dispatch (see PolyaGammaSampler::Sample):
 *
 *  - b integer and b <= kExactIntegerMaxB (50): exact. Sum of b independent
 *    PG(1, z) draws from Devroye's sampler [D09, PSW13 §4].
 *  - otherwise, b < kNormalMinB (10,000): truncated gamma series. PG(b, z) is
 *    (1 / (2 pi^2)) sum_k g_k / ((k - 1/2)^2 + z^2 / (4 pi^2)), g_k iid
 *    Gamma(b, 1) [PSW13 eq. (1)]. The first kSeriesTerms (100) terms are
 *    drawn; the remainder is replaced by its exact expectation, computed as
 *    the closed-form mean minus the mean of the drawn terms. The mean is
 *    therefore exact, and the variance omitted by the truncation is
 *    b / (4 pi^4) sum_{k > K} (k - 1/2)^-4 < 2 / (pi^4 K^3) of the total, about
 *    2e-8 at K = 100.
 *  - b >= kNormalMinB: normal approximation with the exact mean and variance
 *    [WPS14]. The skewness of PG(b, z) is largest at z = 0, where it equals
 *    1.96 / sqrt(b); it is below 0.02 at the threshold.
 *
 * Randomness comes from the package's std::mt19937 through the
 * platform-independent helpers in distributions.h.
 */
#ifndef STOCHTREE_POLYA_GAMMA_SAMPLER_H_
#define STOCHTREE_POLYA_GAMMA_SAMPLER_H_

#include <stochtree/distributions.h>

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace StochTree {

class PolyaGammaSampler {
 public:
  static constexpr double kPi = 3.141592653589793238462643383279502884;
  /*! \brief Devroye's truncation point for the J*(1, c) mixture proposal [D09] */
  static constexpr double kTruncation = 0.64;
  /*! \brief Largest integer b drawn exactly as a sum of PG(1, z) variates */
  static constexpr double kExactIntegerMaxB = 50.0;
  /*! \brief Smallest b drawn from the normal approximation */
  static constexpr double kNormalMinB = 10000.0;
  /*! \brief Number of gamma terms drawn in the truncated series */
  static constexpr int kSeriesTerms = 100;

  PolyaGammaSampler() {}
  ~PolyaGammaSampler() {}

  /*! \brief Draw PG(b, z). b = 0 returns 0 (the degenerate limit); b < 0 or non-finite inputs throw. */
  double Sample(double b, double z, std::mt19937& gen) {
    if (!std::isfinite(b) || !std::isfinite(z) || b < 0.0) {
      throw std::invalid_argument("PolyaGammaSampler: b must be finite and non-negative, z finite");
    }
    if (b == 0.0) return 0.0;
    if (b == std::floor(b) && b <= kExactIntegerMaxB) {
      return SampleExactSum(static_cast<int>(b), z, gen);
    }
    if (b < kNormalMinB) {
      return SampleGammaSeries(b, z, gen, kSeriesTerms);
    }
    return SampleNormalApproximation(b, z, gen);
  }

  /*! \brief Exact PG(1, z) via Devroye's J*(1, c) sampler, PG(1, z) = J*(1, |z| / 2) / 4 [PSW13 §4] */
  double SamplePG1(double z, std::mt19937& gen) {
    return 0.25 * SampleJStar1(0.5 * std::fabs(z), gen);
  }

  /*! \brief Exact PG(b, z) for integer b >= 1 as a sum of b independent PG(1, z) draws */
  double SampleExactSum(int b, double z, std::mt19937& gen) {
    double total = 0.0;
    for (int i = 0; i < b; i++) total += SamplePG1(z, gen);
    return total;
  }

  /*! \brief PG(b, z) by the gamma series truncated at `num_terms`, remainder replaced by its mean */
  double SampleGammaSeries(double b, double z, std::mt19937& gen, int num_terms) {
    const double c2 = (z * z) / (4.0 * kPi * kPi);
    const double inv_two_pi2 = 1.0 / (2.0 * kPi * kPi);
    double draw = 0.0;
    double drawn_mean_terms = 0.0;
    for (int k = 1; k <= num_terms; k++) {
      const double half = static_cast<double>(k) - 0.5;
      const double denom = half * half + c2;
      draw += sample_gamma(gen, b, 1.0) / denom;
      drawn_mean_terms += 1.0 / denom;
    }
    double remainder_mean = Mean(b, z) - b * inv_two_pi2 * drawn_mean_terms;
    if (remainder_mean < 0.0) remainder_mean = 0.0;
    return inv_two_pi2 * draw + remainder_mean;
  }

  /*! \brief Normal approximation N(E[PG(b, z)], Var[PG(b, z)]), redrawn if non-positive */
  double SampleNormalApproximation(double b, double z, std::mt19937& gen) {
    const double mean = Mean(b, z);
    const double sd = std::sqrt(Variance(b, z));
    double draw;
    do {
      draw = sample_standard_normal(mean, sd, gen);
    } while (draw <= 0.0);
    return draw;
  }

  /*! \brief E[PG(b, z)] = b tanh(z / 2) / (2 z), with limit b / 4 at z = 0 [PSW13] */
  static double Mean(double b, double z) {
    const double az = std::fabs(z);
    if (az < 1e-6) return b * (0.25 - az * az / 48.0);
    return b * std::tanh(0.5 * az) / (2.0 * az);
  }

  /*!
   * \brief Var[PG(b, z)] = b (sinh z - z) / (4 z^3 cosh^2(z / 2)), limit b / 24 at z = 0.
   *
   * For |z| < 1, sinh z - z is summed from its Taylor series to avoid
   * cancellation. For |z| >= 1 the identity sinh z = 2 sinh(z/2) cosh(z/2)
   * gives b (2 tanh(z / 2) - z sech^2(z / 2)) / (4 z^3), which cannot overflow.
   */
  static double Variance(double b, double z) {
    const double az = std::fabs(z);
    if (az < 1.0) {
      // sinh z - z = sum_{k >= 1} z^(2k + 1) / (2k + 1)!; divide by z^3 term by term
      const double z2 = az * az;
      double term = 1.0 / 6.0;
      double series = term;
      for (int k = 2; k < 12; k++) {
        term *= z2 / static_cast<double>((2 * k) * (2 * k + 1));
        series += term;
      }
      const double ch = std::cosh(0.5 * az);
      return b * series / (4.0 * ch * ch);
    }
    const double ch = std::cosh(0.5 * az);
    const double sech2 = std::isfinite(ch) ? 1.0 / (ch * ch) : 0.0;
    return b * (2.0 * std::tanh(0.5 * az) - az * sech2) / (4.0 * az * az * az);
  }

 private:
  /*! \brief Exp(1) variate; standard_uniform_draw_53bit lies in [0, 1), so 1 - u is in (0, 1] */
  static double Exponential(std::mt19937& gen) {
    return -std::log(1.0 - standard_uniform_draw_53bit(gen));
  }

  static double NormalCdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
  }

  /*!
   * \brief Coefficient a_n(x) of the alternating series for the J*(1, 0) density [D09, PSW13 §4].
   * The left representation is used on (0, t], the right one on (t, inf).
   */
  static double SeriesCoefficient(int n, double x) {
    const double k = static_cast<double>(n) + 0.5;
    if (x > kTruncation) {
      return kPi * k * std::exp(-0.5 * k * k * kPi * kPi * x);
    }
    return kPi * k * std::pow(2.0 / (kPi * x), 1.5) * std::exp(-2.0 * k * k / x);
  }

  /*!
   * \brief P(X <= t) for X ~ inverse Gaussian with mean 1 / c and shape 1.
   * Phi((t c - 1) / sqrt(t)) + exp(2 c) Phi(-(t c + 1) / sqrt(t)); the second
   * term is evaluated on the log scale so exp(2 c) cannot overflow. At c = 0 this
   * is the Lévy CDF 2 Phi(-1 / sqrt(t)).
   */
  static double InverseGaussianCdfAtTruncation(double c) {
    const double t = kTruncation;
    const double root_t = std::sqrt(t);
    const double first = NormalCdf((t * c - 1.0) / root_t);
    const double tail = NormalCdf(-(t * c + 1.0) / root_t);
    const double second = (tail > 0.0) ? std::exp(2.0 * c + std::log(tail)) : 0.0;
    return first + second;
  }

  /*!
   * \brief Inverse-Gaussian(mean 1 / c, shape 1) truncated to (0, t].
   *
   * When the mean 1 / c exceeds t (including c = 0), propose from the Lévy
   * law truncated to (0, t] and accept with probability exp(-c^2 x / 2). A
   * Lévy draw is 1 / Z^2 with Z standard normal; X <= t is Z beyond 1 / sqrt(t),
   * drawn by the exponential-proposal method for normal tails. Otherwise draw
   * the untruncated inverse Gaussian [MSH76] and reject values above t.
   */
  double SampleTruncatedInverseGaussian(double c, std::mt19937& gen) {
    const double t = kTruncation;
    if (c * t < 1.0) {
      while (true) {
        double e1, e2;
        do {
          e1 = Exponential(gen);
          e2 = Exponential(gen);
        } while (e1 * e1 > 2.0 * e2 / t);
        const double denom = 1.0 + t * e1;
        const double x = t / (denom * denom);
        if (standard_uniform_draw_53bit(gen) < std::exp(-0.5 * c * c * x)) return x;
      }
    }
    const double mu = 1.0 / c;
    while (true) {
      const double nu = sample_standard_normal(0.0, 1.0, gen);
      const double y = nu * nu;
      // The two roots of the MSH quadratic multiply to mu^2; compute the larger
      // one and divide to get the smaller one without cancellation.
      const double larger = mu + 0.5 * mu * mu * y + 0.5 * mu * std::sqrt(4.0 * mu * y + mu * mu * y * y);
      const double smaller = mu * mu / larger;
      const double x = (standard_uniform_draw_53bit(gen) <= mu / (mu + smaller)) ? smaller : larger;
      if (x <= t) return x;
    }
  }

  /*!
   * \brief Exact J*(1, c) draw by Devroye's alternating-series rejection sampler.
   *
   * The proposal mixes an exponential with rate K = pi^2 / 8 + c^2 / 2 on
   * (t, inf), mass p = pi / (2 K) exp(-K t), and the inverse Gaussian
   * IG(1 / c, 1) truncated to (0, t], mass q = 2 exp(-c) P(IG <= t). A proposal
   * x is accepted when U a_0(x) falls below the alternating partial sums of
   * sum_n (-1)^n a_n(x), which bracket the target density [D09].
   */
  double SampleJStar1(double c, std::mt19937& gen) {
    const double t = kTruncation;
    const double K = kPi * kPi / 8.0 + 0.5 * c * c;
    const double p = kPi / (2.0 * K) * std::exp(-K * t);
    const double q = 2.0 * std::exp(-c) * InverseGaussianCdfAtTruncation(c);
    const double prob_right = p / (p + q);
    while (true) {
      double x;
      if (standard_uniform_draw_53bit(gen) < prob_right) {
        x = t + Exponential(gen) / K;
      } else {
        x = SampleTruncatedInverseGaussian(c, gen);
      }
      double s = SeriesCoefficient(0, x);
      const double y = standard_uniform_draw_53bit(gen) * s;
      int n = 0;
      while (true) {
        n++;
        if (n % 2 == 1) {
          s -= SeriesCoefficient(n, x);
          if (y <= s) return x;
        } else {
          s += SeriesCoefficient(n, x);
          if (y > s) break;
        }
      }
    }
  }
};

}  // namespace StochTree

#endif  // STOCHTREE_POLYA_GAMMA_SAMPLER_H_
