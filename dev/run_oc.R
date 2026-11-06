# Operating characteristics of both designs for all scenarios, stored as
# the package dataset oc_results (10,000 trials per scenario and calibration;
# each run takes up to about 15 minutes).

devtools::load_all()
lambda_c <- log(2) / 12
a0 <- 1; b0 <- 12 / log(2)
c_eff <- 0.9873009
D <- 247

# each scenario: hazard rates and change points for both arms
scenarios <- list(
  h0       = list(rates_c = lambda_c, breaks_c = numeric(0),
                  rates_t = lambda_c, breaks_t = numeric(0)),
  h0_break = list(rates_c = c(1, 0.6) * lambda_c , breaks_c = 6,    # both arms identical, with change point
                  rates_t = c(1, 0.6) * lambda_c, breaks_t = 6),
  ph_07    = list(rates_c = lambda_c, breaks_c = numeric(0),
                  rates_t = 0.7 * lambda_c, breaks_t = numeric(0)),
  delayed  = list(rates_c = lambda_c , breaks_c = numeric(0),
                   rates_t = c(1, 0.6) * lambda_c, breaks_t = 6),
  diminish = list(rates_c = lambda_c , breaks_c = numeric(0),
                  rates_t = c(0.5, 1) * lambda_c, breaks_t = 6),
  crossing = list(rates_c = lambda_c , breaks_c = numeric(0),
                  rates_t = c(1.3, 0.5) * lambda_c, breaks_t = 6)
)

# 1) Constant calibration (c_eff from dev/calibrate.R, step 1).
#    Computed once; the result is stored in dev/oc_const.rds.
# set.seed(2031)
# oc_const <- operating_characteristics(scenarios, n_sim = 10000,
#                                       c_eff = c_eff, a0 = a0, b0 = b0)
# saveRDS(oc_const, "dev/oc_const.rds")
oc_const <- readRDS("dev/oc_const.rds")

# 2) Robust calibration (c_robust from dev/calibrate.R, step 3)
c_robust <- 0.9918179602
set.seed(2033)
oc_robust <- operating_characteristics(scenarios, n_sim = 10000,
                                       c_eff = c_robust, a0 = a0, b0 = b0)
saveRDS(oc_robust, "dev/oc_robust.rds")

# 3) Combine both and save as the package dataset
oc_results <- rbind(cbind(calibration = "constant", oc_const),
                    cbind(calibration = "robust",   oc_robust))
usethis::use_data(oc_results, overwrite = TRUE)

