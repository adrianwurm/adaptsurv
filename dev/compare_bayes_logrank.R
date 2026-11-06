# Compare the Bayesian posterior probability P(HR < 1 | data) with
# 1 - one-sided log-rank p-value on a single data cut (weak prior:
# both numbers are nearly identical).

devtools::load_all()
set.seed(42)
rc <- c(0.10, 0.06)

trial <- simulate_trial(n = 400, accrual_time = 12,
                        rates_c = rc,       breaks_c = 6,
                        rates_t = 0.85 * rc, breaks_t = 6)
x <- cut_at_events(trial, 240)

print(analyse_bayes(x$data, 1, 12/log(2))[["prob"]])
print(1-(analyse_fixed(x$data)[["p_one_sided"]]))
