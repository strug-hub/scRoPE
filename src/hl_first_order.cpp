#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct ModeValue {
  bool ok;
  double h_eta;
  double h_etaeta;
};

inline bool finite_number(const double x) {
  return std::isfinite(x);
}

inline double logistic(const double x) {
  return R::plogis(x, 0.0, 1.0, 1, 0);
}

inline double logspace_add(const double a, const double b) {
  const double z = std::max(a, b);
  return z + std::log1p(std::exp(-std::abs(a - b)));
}

inline double d_minus_expm1(const double d) {
  if (std::abs(d) >= 1e-3) {
    return d - std::expm1(d);
  }
  double term = d * d / 2.0;
  double total = term;
  for (int kk = 3; kk <= 12; ++kk) {
    term *= d / static_cast<double>(kk);
    total += term;
  }
  return -total;
}

inline double gamma_log_normalizer(const double alpha,
                                   const double log_alpha) {
  if (alpha <= 1e5) {
    return alpha * log_alpha - R::lgammafn(alpha) - alpha;
  }
  const double inv = 1.0 / alpha;
  const double inv2 = inv * inv;
  const double inv3 = inv2 * inv;
  const double inv5 = inv3 * inv2;
  const double inv7 = inv5 * inv2;
  const double inv9 = inv7 * inv2;
  const double remainder = inv / 12.0 - inv3 / 360.0 +
    inv5 / 1260.0 - inv7 / 1680.0 + inv9 / 1188.0;
  return 0.5 * (log_alpha - std::log(2.0 * M_PI)) - remainder;
}

inline double gamma_prior_value(const double eta, const double v,
                                const double alpha,
                                const double log_alpha) {
  const double d = eta - v / 2.0;
  return alpha * d_minus_expm1(d) +
    gamma_log_normalizer(alpha, log_alpha);
}

inline double softplus(const double x) {
  return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x)));
}

inline double lgamma_increment(const double phi, const double y) {
  if (phi <= 1e5 || y != std::floor(y) || y > 1e5) {
    return R::lgammafn(phi + y) - R::lgammafn(phi);
  }
  if (y == 0.0) {
    return 0.0;
  }
  double correction = 0.0;
  const int count = static_cast<int>(y);
  for (int kk = 0; kk < count; ++kk) {
    correction += std::log1p(static_cast<double>(kk) / phi);
  }
  return y * std::log(phi) + correction;
}

ModeValue mode_value_subject(
    const Rcpp::NumericVector& y,
    const std::vector<double>& log_m_base,
    const std::vector<std::vector<int> >& rows,
    const int subject,
    const double eta,
    const double alpha,
    const double v,
    const double phi,
    const double xi) {
  const double d = eta - v / 2.0;
  const double exp_d = std::exp(d);
  const double lambda_z = alpha * exp_d;
  if (!finite_number(lambda_z)) {
    return ModeValue{false, NA_REAL, NA_REAL};
  }

  double score = 0.0;
  double curvature = 0.0;
  const std::vector<int>& subject_rows = rows[subject];
  for (std::size_t jj = 0; jj < subject_rows.size(); ++jj) {
    const int cell = subject_rows[jj];
    const double log_m = log_m_base[cell] + eta;
    const double r = logistic(log_m - xi);
    const double one_minus_r = logistic(xi - log_m);
    const double s = y[cell] * one_minus_r - phi * r;
    const double tt = -(y[cell] + phi) * r * one_minus_r;
    if (!finite_number(s) || !finite_number(tt)) {
      return ModeValue{false, NA_REAL, NA_REAL};
    }
    score += s;
    curvature += tt;
  }

  return ModeValue{
    true,
    score - alpha * std::expm1(d),
    curvature - lambda_z
  };
}

bool nb_cell_primitives(
    const double y,
    const double log_m,
    const double xi,
    const double phi,
    double& loglik,
    double& score,
    double& curvature) {
  const double r = logistic(log_m - xi);
  const double one_minus_r = logistic(xi - log_m);
  const double log_d = logspace_add(log_m, xi);
  const double inv_d = std::exp(-log_d);

  const double log_odds = log_m - xi;
  loglik = lgamma_increment(phi, y) - R::lgammafn(y + 1.0) -
    phi * softplus(log_odds) - y * softplus(-log_odds);
  score = y * one_minus_r - phi * r;
  curvature = -(y + phi) * r * one_minus_r;
  const double third = -(y + phi) * r * one_minus_r * (1.0 - 2.0 * r);
  const double fourth = -(y + phi) * r * one_minus_r *
    (1.0 - 6.0 * r + 6.0 * r * r);

  const double l_phi = R::digamma(y + phi) - R::digamma(phi) +
    xi + 1.0 - log_d - (y + phi) * inv_d;
  const double l_phiphi = R::trigamma(y + phi) - R::trigamma(phi) +
    1.0 / phi - inv_d + y * inv_d * inv_d - r * inv_d;
  const double s_phi = r * y * one_minus_r / phi - r * r;
  const double s_phiphi = -2.0 * s_phi * inv_d;
  const double t_phi = -r * one_minus_r +
    (y + phi) * r * one_minus_r * (1.0 - 2.0 * r) / phi;
  const double u_phi = -r * one_minus_r * (1.0 - 2.0 * r) +
    (y + phi) * r * one_minus_r *
      (1.0 - 6.0 * r + 6.0 * r * r) / phi;
  const double t_phiphi = -r * (
    2.0 * inv_d - 4.0 * (y + 2.0 * phi) * inv_d * inv_d +
      6.0 * phi * (y + phi) * inv_d * inv_d * inv_d
  );

  const double l_xi = phi * l_phi;
  const double l_xixi = phi * l_phi + phi * phi * l_phiphi;
  const double s_xi = phi * s_phi;
  const double s_xixi = phi * s_phi + phi * phi * s_phiphi;
  const double t_xi = phi * t_phi;
  const double t_xixi = phi * t_phi + phi * phi * t_phiphi;
  const double u_xi = phi * u_phi;

  const double values[] = {
    phi, r, loglik, score, curvature, third, fourth, l_phi, l_phiphi,
    s_phi, s_phiphi, t_phi, t_phiphi, u_phi, l_xi, l_xixi, s_xi,
    s_xixi, t_xi, t_xixi, u_xi
  };
  const std::size_t n_values = sizeof(values) / sizeof(values[0]);
  for (std::size_t ii = 0; ii < n_values; ++ii) {
    if (!finite_number(values[ii])) {
      return false;
    }
  }
  return true;
}

Rcpp::List mode_bracket_failure(const std::vector<bool>& unresolved) {
  int first_bad = -1;
  int n_bad = 0;
  for (std::size_t ii = 0; ii < unresolved.size(); ++ii) {
    if (unresolved[ii]) {
      if (first_bad < 0) {
        first_bad = static_cast<int>(ii);
      }
      ++n_bad;
    }
  }
  Rcpp::IntegerVector unresolved_subjects(n_bad);
  int at = 0;
  for (std::size_t ii = 0; ii < unresolved.size(); ++ii) {
    if (unresolved[ii]) {
      unresolved_subjects[at++] = static_cast<int>(ii) + 1;
    }
  }
  return Rcpp::List::create(
    Rcpp::_["success"] = false,
    Rcpp::_["status"] = "mode_failure",
    Rcpp::_["reason"] = "mode_bracket_failure",
    Rcpp::_["failed_subject"] = first_bad + 1,
    Rcpp::_["diagnostics"] = Rcpp::List::create(
      Rcpp::_["unresolved_subjects"] = unresolved_subjects
    )
  );
}

}  // namespace


// Evaluate the complete ordinary-h inner transaction in canonical
// (beta, tau, xi) coordinates.  subject_zero must use contiguous zero-based
// subject indices in first-appearance order.
//
// [[Rcpp::export]]
Rcpp::List hl_fo_inner_eval_cpp(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& X,
    const Rcpp::NumericVector& log_offset,
    const Rcpp::IntegerVector& subject_zero,
    const Rcpp::NumericVector& beta,
    const double tau,
    const double xi,
    const Rcpp::NumericVector& eta_start,
    const bool use_warm) {
  const int n = y.size();
  const int p = X.ncol();
  if (X.nrow() != n || log_offset.size() != n || subject_zero.size() != n) {
    Rcpp::stop("y, X, log_offset, and subject_zero must describe the same cells.");
  }
  if (beta.size() != p) {
    Rcpp::stop("beta must have one entry per design column.");
  }
  if (n < 1 || p < 1) {
    Rcpp::stop("The inner transaction requires at least one cell and one beta coordinate.");
  }

  int max_subject = -1;
  for (int ii = 0; ii < n; ++ii) {
    if (!finite_number(y[ii]) || y[ii] < 0.0 || y[ii] != std::round(y[ii])) {
      Rcpp::stop("Counts must be finite nonnegative integers.");
    }
    if (!finite_number(log_offset[ii]) ||
        Rcpp::IntegerVector::is_na(subject_zero[ii]) || subject_zero[ii] < 0) {
      Rcpp::stop("Offsets must be finite and subject_zero must be nonnegative.");
    }
    max_subject = std::max(max_subject, subject_zero[ii]);
    for (int aa = 0; aa < p; ++aa) {
      if (!finite_number(X(ii, aa))) {
        Rcpp::stop("Design entries must be finite.");
      }
    }
  }
  for (int aa = 0; aa < p; ++aa) {
    if (!finite_number(beta[aa])) {
      Rcpp::stop("beta must be finite.");
    }
  }

  const int m = max_subject + 1;
  std::vector<std::vector<int> > rows(m);
  std::vector<int> first_seen(m, -1);
  int next_subject = 0;
  for (int ii = 0; ii < n; ++ii) {
    const int subject = subject_zero[ii];
    if (first_seen[subject] < 0) {
      if (subject != next_subject) {
        Rcpp::stop("subject_zero must be contiguous and ordered by first appearance.");
      }
      first_seen[subject] = ii;
      ++next_subject;
    }
    rows[subject].push_back(ii);
  }
  for (int ii = 0; ii < m; ++ii) {
    if (rows[ii].empty()) {
      Rcpp::stop("Every subject must contribute at least one cell.");
    }
  }

  std::vector<double> log_m_base(n, 0.0);
  for (int ii = 0; ii < n; ++ii) {
    double value = log_offset[ii];
    for (int aa = 0; aa < p; ++aa) {
      value += X(ii, aa) * beta[aa];
    }
    log_m_base[ii] = value;
  }

  const double v = std::exp(tau);
  const double evm1 = std::expm1(v);
  const double alpha = 1.0 / evm1;
  const double log_alpha = -std::log(evm1);
  const double phi = std::exp(xi);

  bool warm = use_warm && eta_start.size() == m;
  if (warm) {
    for (int ii = 0; ii < m; ++ii) {
      if (!finite_number(eta_start[ii])) {
        warm = false;
        break;
      }
    }
  }
  std::vector<double> start(m, v / 2.0);
  if (warm) {
    for (int ii = 0; ii < m; ++ii) {
      start[ii] = eta_start[ii];
    }
  }

  std::vector<double> width0(m);
  std::vector<double> lo(m, NA_REAL);
  std::vector<double> hi(m, NA_REAL);
  std::vector<int> expansions(m, NA_INTEGER);
  std::vector<bool> unresolved(m, true);
  for (int ii = 0; ii < m; ++ii) {
    width0[ii] = std::max(1.0, 0.25 * (1.0 + std::abs(start[ii])));
  }

  for (int ee = 0; ee <= 60; ++ee) {
    const double multiplier = std::ldexp(1.0, ee);
    for (int ii = 0; ii < m; ++ii) {
      if (!unresolved[ii]) {
        continue;
      }
      const double width = width0[ii] * multiplier;
      const double candidate_lo = start[ii] - width;
      const double candidate_hi = start[ii] + width;
      const ModeValue value_lo = mode_value_subject(
        y, log_m_base, rows, ii, candidate_lo, alpha, v, phi, xi
      );
      const ModeValue value_hi = mode_value_subject(
        y, log_m_base, rows, ii, candidate_hi, alpha, v, phi, xi
      );
      if (value_lo.ok && value_hi.ok && value_lo.h_eta >= 0.0 &&
          value_hi.h_eta <= 0.0) {
        lo[ii] = candidate_lo;
        hi[ii] = candidate_hi;
        expansions[ii] = ee;
        unresolved[ii] = false;
      }
    }
    bool any_unresolved = false;
    for (int ii = 0; ii < m; ++ii) {
      any_unresolved = any_unresolved || unresolved[ii];
    }
    if (!any_unresolved) {
      break;
    }
  }
  for (int ii = 0; ii < m; ++ii) {
    if (unresolved[ii]) {
      return mode_bracket_failure(unresolved);
    }
  }

  std::vector<double> eta(m);
  std::vector<int> bisections(m, 0);
  std::vector<int> iterations(m, 0);
  std::vector<bool> converged(m, false);
  std::vector<ModeValue> current(m);
  std::vector<ModeValue> last(m);
  bool have_last = false;
  for (int ii = 0; ii < m; ++ii) {
    eta[ii] = std::min(std::max(start[ii], lo[ii]), hi[ii]);
  }

  for (int iter = 1; iter <= 100; ++iter) {
    for (int ii = 0; ii < m; ++ii) {
      current[ii] = mode_value_subject(
        y, log_m_base, rows, ii, eta[ii], alpha, v, phi, xi
      );
      if (!current[ii].ok) {
        return Rcpp::List::create(
          Rcpp::_["success"] = false,
          Rcpp::_["status"] = "mode_failure",
          Rcpp::_["reason"] = "nonfinite_mode_arithmetic"
        );
      }
    }
    for (int ii = 0; ii < m; ++ii) {
      if (!finite_number(current[ii].h_etaeta) ||
          current[ii].h_etaeta >= 0.0) {
        return Rcpp::List::create(
          Rcpp::_["success"] = false,
          Rcpp::_["status"] = "mode_failure",
          Rcpp::_["reason"] = "mode_nonnegative_curvature",
          Rcpp::_["failed_subject"] = ii + 1
        );
      }
    }

    bool all_converged = true;
    std::vector<double> residual(m);
    std::vector<double> step_residual(m);
    for (int ii = 0; ii < m; ++ii) {
      residual[ii] = std::abs(current[ii].h_eta) /
        std::sqrt(-current[ii].h_etaeta);
      step_residual[ii] = std::abs(current[ii].h_eta / current[ii].h_etaeta) /
        (1.0 + std::abs(eta[ii]));
      if (!converged[ii] && residual[ii] <= 1e-10 &&
          step_residual[ii] <= 1e-10) {
        converged[ii] = true;
        iterations[ii] = iter;
      }
      all_converged = all_converged && converged[ii];
    }
    if (all_converged) {
      Rcpp::NumericVector eta_out(m);
      Rcpp::List diagnostics(m);
      for (int ii = 0; ii < m; ++ii) {
        eta_out[ii] = eta[ii];
        diagnostics[ii] = Rcpp::List::create(
          Rcpp::_["success"] = true,
          Rcpp::_["status"] = "ok",
          Rcpp::_["reason"] = R_NilValue,
          Rcpp::_["eta"] = eta[ii],
          Rcpp::_["h_eta"] = current[ii].h_eta,
          Rcpp::_["h_etaeta"] = current[ii].h_etaeta,
          Rcpp::_["residual_scaled"] = residual[ii],
          Rcpp::_["step_scaled"] = step_residual[ii],
          Rcpp::_["iterations"] = iterations[ii],
          Rcpp::_["bracket_expansions"] = expansions[ii],
          Rcpp::_["bisections"] = bisections[ii],
          Rcpp::_["bracket"] = Rcpp::NumericVector::create(lo[ii], hi[ii]),
          Rcpp::_["start_type"] = warm ? "warm" : "cold_prior_mode"
        );
      }

      std::vector<double> nb_loglik(n);
      std::vector<double> nb_score(n);
      std::vector<double> nb_curvature(n);
      for (int cell = 0; cell < n; ++cell) {
        const int subject = subject_zero[cell];
        const double log_m = log_m_base[cell] + eta[subject];
        if (!finite_number(phi) || phi <= 0.0 ||
            !nb_cell_primitives(
              y[cell], log_m, xi, phi, nb_loglik[cell], nb_score[cell],
              nb_curvature[cell]
            )) {
          return Rcpp::List::create(
            Rcpp::_["success"] = false,
            Rcpp::_["status"] = "mode_failure",
            Rcpp::_["reason"] = !finite_number(phi) || phi <= 0.0 ?
              "nonfinite_nb_domain" : "nonfinite_nb_arithmetic"
          );
        }
      }

      std::vector<double> h_subject(m, 0.0);
      std::vector<double> a_subject(m, 0.0);
      std::vector<double> lambda_z(m);
      for (int ii = 0; ii < m; ++ii) {
        lambda_z[ii] = alpha * std::exp(eta[ii] - v / 2.0);
        const double prior_value = gamma_prior_value(
          eta[ii], v, alpha, log_alpha
        );
        h_subject[ii] = prior_value;
        a_subject[ii] = -lambda_z[ii];
      }

      Rcpp::NumericVector score_beta(p);
      Rcpp::NumericMatrix h_beta_subject(m, p);
      Rcpp::NumericMatrix h_bb(p, p);
      Rcpp::NumericMatrix d_subject(m, p);
      for (int cell = 0; cell < n; ++cell) {
        const int subject = subject_zero[cell];
        h_subject[subject] += nb_loglik[cell];
        a_subject[subject] += nb_curvature[cell];
        for (int aa = 0; aa < p; ++aa) {
          const double x_a = X(cell, aa);
          const double score_term = x_a * nb_score[cell];
          score_beta[aa] += score_term;
          h_beta_subject(subject, aa) += score_term;
          d_subject(subject, aa) += x_a * nb_curvature[cell];
          for (int bb = 0; bb < p; ++bb) {
            h_bb(aa, bb) += x_a * X(cell, bb) * nb_curvature[cell];
          }
        }
      }
      for (int ii = 0; ii < m; ++ii) {
        if (!finite_number(a_subject[ii]) || a_subject[ii] >= 0.0) {
          return Rcpp::List::create(
            Rcpp::_["success"] = false,
            Rcpp::_["status"] = "mode_failure",
            Rcpp::_["reason"] = "mode_nonnegative_curvature",
            Rcpp::_["failed_subject"] = ii + 1
          );
        }
      }

      Rcpp::NumericMatrix B_beta(p, p);
      for (int aa = 0; aa < p; ++aa) {
        for (int bb = 0; bb < p; ++bb) {
          double value = h_bb(aa, bb);
          for (int ii = 0; ii < m; ++ii) {
            value -= d_subject(ii, aa) * d_subject(ii, bb) / a_subject[ii];
          }
          B_beta(aa, bb) = value;
        }
      }
      const SEXP x_dimnames_sexp = X.attr("dimnames");
      if (!Rf_isNull(x_dimnames_sexp)) {
        const Rcpp::List x_dimnames(x_dimnames_sexp);
        if (x_dimnames.size() == 2 && !Rf_isNull(x_dimnames[1])) {
          B_beta.attr("dimnames") = Rcpp::List::create(
            x_dimnames[1],
            x_dimnames[1]
          );
        }
      }

      double P = 0.0;
      Rcpp::List subject_output(m);
      for (int ii = 0; ii < m; ++ii) {
        P += h_subject[ii];
        Rcpp::NumericVector h_beta(p);
        for (int aa = 0; aa < p; ++aa) {
          h_beta[aa] = h_beta_subject(ii, aa);
        }
        subject_output[ii] = Rcpp::List::create(
          Rcpp::_["p"] = h_subject[ii],
          Rcpp::_["h_beta"] = h_beta,
          Rcpp::_["p_beta_beta"] = Rcpp::LogicalVector::create(NA_LOGICAL),
          Rcpp::_["a"] = a_subject[ii]
        );
      }

      return Rcpp::List::create(
        Rcpp::_["success"] = true,
        Rcpp::_["status"] = "ok",
        Rcpp::_["beta"] = Rcpp::clone(beta),
        Rcpp::_["eta"] = eta_out,
        Rcpp::_["mode_diagnostics"] = diagnostics,
        Rcpp::_["P"] = P,
        Rcpp::_["S_beta"] = score_beta,
        Rcpp::_["B_beta"] = B_beta,
        Rcpp::_["subject"] = subject_output
      );
    }

    for (int ii = 0; ii < m; ++ii) {
      if (converged[ii]) {
        continue;
      }
      if (current[ii].h_eta > 0.0) {
        lo[ii] = eta[ii];
      } else {
        hi[ii] = eta[ii];
      }
      double proposal = eta[ii] - current[ii].h_eta / current[ii].h_etaeta;
      if (!finite_number(proposal) || proposal <= lo[ii] || proposal >= hi[ii]) {
        proposal = (lo[ii] + hi[ii]) / 2.0;
        ++bisections[ii];
      }
      eta[ii] = proposal;
    }
    last = current;
    have_last = true;
  }

  int bad = 0;
  while (bad < m && converged[bad]) {
    ++bad;
  }
  return Rcpp::List::create(
    Rcpp::_["success"] = false,
    Rcpp::_["status"] = "mode_failure",
    Rcpp::_["reason"] = "mode_iteration_exhausted",
    Rcpp::_["failed_subject"] = bad + 1,
    Rcpp::_["eta"] = Rcpp::wrap(eta),
    Rcpp::_["h_eta"] = have_last ? last[bad].h_eta : NA_REAL,
    Rcpp::_["h_etaeta"] = have_last ? last[bad].h_etaeta : NA_REAL
  );
}
