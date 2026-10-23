#include <Rcpp.h>
#include <cmath>

namespace {

// log1p(n/a), without overflowing n/a for finite positive a and n.
double log1p_ratio(const double n, const double a) {
  if (n <= a) return std::log1p(n / a);
  const double ratio = n / a;
  if (R_FINITE(ratio)) return std::log1p(ratio);
  return std::log(n) - std::log(a) + std::log1p(a / n);
}

// lgamma(a+n) - lgamma(a) - n*log(a). Small integer counts use the
// normalized rising product. All other counts use an equivalent shifted
// Stirling difference; its work and storage do not grow with the count.
// Noninteger n is retained because prepared donor totals are doubles.
double rising_normalized(const double a, const double n) {
  if (!R_FINITE(a) || a <= 0 || !R_FINITE(n) || n < 0) return R_NaN;
  if (n == 0 || n == 1) return 0;
  if (n <= 32 && n == std::floor(n)) {
    double value = 0;
    for (int j = 1; j < static_cast<int>(n); ++j)
      value += log1p_ratio(static_cast<double>(j), a);
    return value;
  }

  // At z >= 32, retaining B_2 through B_16 makes the absolute remainder
  // in each log-gamma expansion < (43867/244188)/32^17 < 4.7e-27.
  const int shift = a < 32 ? static_cast<int>(std::ceil(32 - a)) : 0;
  const double z = a + shift;
  const double r = n / z;
  const double log_ratio = log1p_ratio(n, z);
  double mean_term;
  if (r <= 0.125) {
    // ((1+r)*log1p(r)-r)/r, summed without subtracting O(r) terms.
    // The omitted terms at r <= 1/8 are below double precision.
    double power = r;
    mean_term = 0;
    for (int j = 2; j <= 24; ++j) {
      mean_term += power / (static_cast<double>(j) * (j - 1));
      power *= -r;
    }
  } else {
    mean_term = (1 + z / n) * log_ratio - 1;
  }
  double value = n * mean_term - 0.5 * log_ratio;
  const double coefficients[] = {
    1.0 / 12, -1.0 / 360, 1.0 / 1260, -1.0 / 1680,
    1.0 / 1188, -691.0 / 360360, 1.0 / 156, -3617.0 / 122400
  };
  const double inverse = 1 / z;
  double inverse_power = inverse;
  for (int j = 0; j < 8; ++j) {
    const int power = 2 * j + 1;
    value += coefficients[j] * inverse_power * std::expm1(-power * log_ratio);
    inverse_power *= inverse * inverse;
  }
  if (shift) {
    value += n * log1p_ratio(static_cast<double>(shift), a);
    for (int j = 0; j < shift; ++j)
      value -= log1p_ratio(n, a + j);
  }
  return value;
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericVector ln_joint_rising_normalized_cpp(
    const double a, const Rcpp::NumericVector& counts) {
  Rcpp::NumericVector out(counts.size());
  for (R_xlen_t i = 0; i < counts.size(); ++i)
    out[i] = rising_normalized(a, counts[i]);
  return out;
}
