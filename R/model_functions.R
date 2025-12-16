
#' Fit age-dependent copper model
#'
#' Fits a log-linear model of serum copper as a function of age.
#'
#' @param data Data frame containing copper and age
#' @param cu_var Name of copper variable (µg/L)
#' @param age_var Name of age variable (months)
#'
#' @return An object of class 'cu_age_model'
#' @export
fit_age_model <- function(data, cu_var = "Cu", age_var = "age.months") {

  names(data)[names(data) == cu_var] <- "Cu"
  names(data)[names(data) == age_var] <- "age.months"

  formula <- log(Cu) ~ sqrt(age.months) + I(sqrt(age.months)^2)

  model <- lm(formula, data = data)

  out <- list(
    model = model,
    formula = formula,
    coefficients = coef(model),
    residual_sd = sd(residuals(model)),
    n = nrow(data)
  )

  class(out) <- "cu_age_model"
  return(out)
}


#' Fit age-adjusted reference intervals for serum copper
#'
#' Estimates lower and upper reference limits using bootstrap residuals.
#'
#' @param age_model Object returned by fit_age_model()
#' @param data Data used to estimate reference intervals
#' @param n_boot Number of bootstrap iterations
#'
#' @return An object of class 'cu_reference_intervals'
#' @export
fit_reference_intervals <- function(age_model, data, n_boot = 1000) {

  stopifnot(inherits(age_model, "cu_age_model"))

  model <- age_model$model
  kl_vals <- numeric(n_boot)
  ku_vals <- numeric(n_boot)

  set.seed(123)

  for (i in seq_len(n_boot)) {
    boot_data <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
    boot_lm <- boot_lm <- lm(formula = age_model$formula, data = boot_data)
    res <- residuals(boot_lm)

    kl_vals[i] <- quantile(res, 0.025, na.rm = TRUE)
    ku_vals[i] <- quantile(res, 0.975, na.rm = TRUE)
  }

  out <- list(
    Kl = mean(kl_vals),
    Ku = mean(ku_vals),
    n_boot = n_boot
  )

  class(out) <- "cu_reference_intervals"
  return(out)
}

#' Flag for systemic acute inflammation
#'
#' Computes a logical flag for systemic inflammation based on CRP, ESR and fibrinogen.
#'
#' @param data Data frame containing the marker columns
#' @param crp_col Name of CRP column (default "CRP")
#' @param esr_col Name of ESR column (default "ESR")
#' @param fibri_col Name of fibrinogen column (default "FIBRI")
#' @param crp_cutoff Numeric, threshold for CRP (default 10 mg/L)
#' @param esr_cutoff Numeric, threshold for ESR (default 30 mm/h)
#' @param fibri_cutoff Numeric, threshold for fibrinogen (default 4 g/L)
#'
#' @return Logical vector, TRUE if any marker is above threshold
#' @export
inflammation_flag <- function(data,
                              crp_col = "CRP",
                              esr_col = "ESR",
                              fibri_col = "fibrinogen",
                              crp_cutoff = 10,
                              esr_cutoff = 30,
                              fibri_cutoff = 4) {

  stopifnot(all(c(crp_col, esr_col, fibri_col) %in% colnames(data)))

  with(data,
       ( .subset2(data, crp_col) > crp_cutoff ) |
         ( .subset2(data, esr_col) > esr_cutoff ) |
         ( .subset2(data, fibri_col) > fibri_cutoff )
  )
}


#' Fit inflammation score using PLS
#'
#' @param age_model Object returned by fit_age_model()
#' @param data Data frame with inflammation markers and copper residuals
#' @param fibrinogen_col Name of fibrinogen column
#' @param crp_col Name of CRP column
#' @param esr_col Name of ESR column
#' @param cu_var Name copper variable
#' @param age_var Name of age variable (months)
#'
#' @return Object of class 'cu_inflammation_score'
#' @export
fit_inflammation_score <- function(age_model,
                                   data,
                                   fibrinogen_col = 'fibrinogen',
                                   crp_col = 'CRP',
                                   esr_col = 'ESR',
                                   cu_var = "Cu",
                                   age_var = "age.months") {

  # Renombrar internamente las columnas de inflamación
  data_renamed <- data
  names(data_renamed)[names(data_renamed) == fibrinogen_col] <- "fibrinogen"
  names(data_renamed)[names(data_renamed) == crp_col] <- "CRP"
  names(data_renamed)[names(data_renamed) == esr_col] <- "ESR"

  # Mantener el orden consistente
  marker_names <- c("fibrinogen", "CRP", "ESR")
  X <- data_renamed[, marker_names, drop = FALSE]

  X_scaled <- scale(X, center = FALSE, scale = TRUE)
  pls_df <- as.data.frame(X_scaled)

  log_cu_pred <- age_model$coefficients["(Intercept)"] +
    age_model$coefficients["sqrt(age.months)"] * sqrt(data[[age_var]]) +
    age_model$coefficients["I(sqrt(age.months)^2)"] * (sqrt(data[[age_var]])^2)

  # Calcular residuo
  cu_residual <- log(data[[cu_var]]) - log_cu_pred
  pls_df$Y <- cu_residual


  # Calcular residuo de cobre ajustado por edad
  log_cu_pred <- age_model$coefficients["(Intercept)"] +
    age_model$coefficients["sqrt(age.months)"] * sqrt(data_renamed[[age_var]]) +
    age_model$coefficients["I(sqrt(age.months)^2)"] * (sqrt(data_renamed[[age_var]])^2)

  cu_residual <- log(data_renamed[[cu_var]]) - log_cu_pred
  pls_df$Y <- cu_residual

  # Ajustar PLS
  pls_model <- pls::plsr(Y ~ ., ncomp = 2, data = pls_df, center = FALSE, verbose = FALSE)

  # score de la primera componente
  infl_score <- as.vector(pls::scores(pls_model)[, 1])

  # Ajuste lineal para obtener B3
  df_lm <- data.frame(
    Y = log(data_renamed[[cu_var]]),
    age_sqrt = sqrt(data_renamed[[age_var]]),
    age_sqrt2 = sqrt(data_renamed[[age_var]])^2,
    infl_score = infl_score
  )

  cu_model_adj_score <- lm(Y ~ age_sqrt + age_sqrt2 + infl_score, data = df_lm)
  B3 <- coef(cu_model_adj_score)["infl_score"]

  # Salida
  out <- list(
    weights = pls_model$loading.weights[, 1],
    sds = attr(X_scaled, "scaled:scale"),
    markers = marker_names,
    B3 = B3
  )

  class(out) <- "cu_inflammation_score"
  return(out)
}


#' Build a complete copper model
#'
#' Combines age model, reference intervals and inflammation score
#' into a single copper_model object.
#'
#' @param age_model Object from fit_age_model()
#' @param reference_intervals Object from fit_reference_intervals()
#' @param inflammation_score Object from fit_inflammation_score()
#' @param metadata Optional list with clinical metadata
#'
#' @return An object of class 'copper_model'
#' @export
build_copper_model <- function(age_model,
                               reference_intervals,
                               inflammation_score = NULL,
                               metadata = list()) {

  stopifnot(inherits(age_model, "cu_age_model"))
  stopifnot(inherits(reference_intervals, "cu_reference_intervals"))

  if (!is.null(inflammation_score)) {
    stopifnot(inherits(inflammation_score, "cu_inflammation_score"))
  }

  model <- list(
    age_model = age_model,
    reference_intervals = reference_intervals,
    inflammation_score = inflammation_score,
    metadata = c(
      list(
        created = format(Sys.Date(), "%d/%m/%Y")
      ),
      metadata
    )
  )

  class(model) <- "copper_model"
  return(model)
}

# Defino el print para que sea autoexplicativo
#' @export
print.copper_model <- function(x, ...) {
  cat("Copper clinical model\n")
  cat("----------------------\n")
  cat("Age model: fitted on", x$age_model$n, "subjects\n")

  if (!is.null(x$inflammation_score)) {
    cat("Inflammation adjustment: YES\n")
    cat("Markers:",
        paste(x$inflammation_score$markers, collapse = ", "),
        "\n")
  } else {
    cat("Inflammation adjustment: NO\n")
  }

  if (length(x$metadata) > 0) {
    cat("\nMetadata:\n")
    for (nm in names(x$metadata)) {
      cat(" -", nm, ":", x$metadata[[nm]], "\n")
    }
  }
}

#' Interpret copper model
#'
#' Provides a clinical interpretation of pediatric serum copper levels by adjusting
#' for age and inflammation using a copper_model object. Returns the observed copper,
#' age-adjusted value, inflammation-adjusted residual, reference intervals, and a final interpretation.
#'
#' @param cu Numeric. Serum copper concentration (µg/L) for the patient.
#' @param age_months Numeric. Age of the patient in months.
#' @param crp Numeric, optional. C-reactive protein (mg/L) for inflammation adjustment.
#' @param esr Numeric, optional. Erythrocyte sedimentation rate (mm/h) for inflammation adjustment.
#' @param fibrinogen Numeric, optional. Fibrinogen (g/L) for inflammation adjustment.
#' @param aat Numeric, optional. Alpha-1 antitrypsin (mg/L) for inflammation adjustment (if used in the model).
#' @param model Object of class 'copper_model' returned by build_copper_model().
#'
#' @return An object of class 'cu_interpretation' containing:
#' \itemize{
#'   \item cu: Observed serum copper value
#'   \item age_months: Patient age in months
#'   \item cu_age_adjusted: Age-adjusted copper
#'   \item cu_residual: Residual from age model
#'   \item cu_residual_adj: Residual adjusted for inflammation
#'   \item inflammation_score: Calculated inflammation score (if applicable)
#'   \item inflammation_flag: TRUE/FALSE if inflammation is present
#'   \item lwr, upr: Lower and upper reference limits for age
#'   \item lwr_adj, upr_adj: Reference limits adjusted for inflammation (if present)
#'   \item interpretation: "LOW", "NORMAL" or "HIGH"
#' }
#' @export
interpret_copper <- function(
                        cu,                  # cobre sérico del paciente (µg/L)
                        age_months,          # edad en meses
                        crp = NULL,          # marcador de inflamación CRP
                        esr = NULL,          # marcador ESR
                        fibrinogen = NULL,   # marcador fibrinógeno
                        model                # objeto copper_model
                    ) {
    stopifnot(inherits(model, "copper_model"))
    stopifnot(!is.null(cu) & !is.null(age_months))


    # AJUSTE POR EDAD
    age_mod <- model$age_model

    log_cu_pred <- age_mod$coefficients["(Intercept)"] +
      age_mod$coefficients["sqrt(age.months)"] * sqrt(age_months) +
      age_mod$coefficients["I(sqrt(age.months)^2)"] * (sqrt(age_months)^2)

    # Calcular residuo
    cu_residual <- log(cu) - log_cu_pred

    # Intervalos de referencia ajustados por edad
    ri <- model$reference_intervals
    lwr <- exp(log_cu_pred + ri$Kl)
    upr <- exp(log_cu_pred + ri$Ku)


    # AJUSTE POR INFLAMACIÓN
    if (!is.null(model$inflammation_score)) {
      stopifnot(!is.null(crp) & !is.null(esr) & !is.null(fibrinogen))

      # Escalar y calcular score
      X <- c(fibrinogen, crp, esr)
      names(X) <- model$inflammation_score$markers
      X_scaled <- X / model$inflammation_score$sds
      infl_score <- sum(X_scaled * model$inflammation_score$weights)

      # Determinar flag de inflamación
      infl_flag <- (crp > 10) | (esr > 30) | (fibrinogen > 4)

      # Ajuste del residuo
      cu_inflam_adj <- if (infl_flag) exp(log_cu_pred + cu_residual - model$inflammation_score$B3 * infl_score) else exp(log_cu_pred + cu_residual)


    } else {
      infl_score <- NA
      infl_flag <- FALSE
      cu_inflam_adj <- cu
    }

    # CLASIFICACIÓN POR LOS RI
    interpretation <- ifelse(cu_inflam_adj < ri$Kl, "LOW",
                             ifelse(cu_inflam_adj > ri$Ku, "HIGH", "NORMAL"))

    # DEVOLVER RESULTADO
    out <- list(
      cu = cu,
      age_months = age_months,
      cu_inflam_adj = unname(cu_inflam_adj),      # <-- quitar nombres
      lwr = unname(lwr),
      upr = unname(upr),
      inflammation_flag = infl_flag,
      inflammation_score = unname(infl_score),
      interpretation = unname(interpretation)
    )

    class(out) <- "cu_interpretation"
    return(out)

  }


#' Imprimir resultado
#' @export
print.cu_interpretation <- function(x, ...) {
  cat("Copper interpretation\n")
  cat("--------------------\n")
  cat("Observed Cu:", x$cu, "µg/L\n")
  if (!is.na(x$inflammation_score)) {
    cat("Inflammation?", ifelse(x$inflammation_flag, "YES", "NO"), "\n")
  } else {
    cat("Inflammation? NA\n")
  }
  cat("Cu (inflammation adjusted):", round(x$cu_inflam_adj, 3), "\n")
  cat("Reference values:", round(x$lwr, 1), "-", round(x$upr, 1), "\n")
  cat("Interpretation:", x$interpretation, "\n")
}
