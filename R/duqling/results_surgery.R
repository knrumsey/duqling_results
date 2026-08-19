load("data/DATA_TESTFUNCS.rda")
sim_study_testfuncs <- duqq

load("data/DATA_REALDATA.rda")
sim_study_realdata <- duqq

table(sim_study_testfuncs$df$method)
table(sim_study_realdata$df$method)

df <- sim_study_testfuncs$df

# Methods that are intentionally run on fewer cases
expensive_methods <- c("gp", "hetgp", "deepgp", "rgasp")

# Expected counts
expected_counts <- setNames(rep(3000L, length(unique(df$method))), unique(df$method))
expected_counts[intersect(names(expected_counts), expensive_methods)] <- 2400L

# Do not fill expensive methods for n_train = 5000
# So the target scenarios for expensive methods are only those with n_train != 5000
all_scenarios <- unique(df$sim_scenario)
baseline_df <- df[df$method == "baseline", ]

# Helpful lookup: scenario -> baseline row
# (should be unique already)
baseline_split <- split(baseline_df, baseline_df$sim_scenario)

rows_to_add <- list()

for (m in unique(df$method)) {
  if (m == "baseline") next

  # Which scenarios should this method have?
  if (m %in% expensive_methods) {
    target_scenarios <- unique(df$sim_scenario[df$method == "baseline" & df$n_train != 5000])
  } else {
    target_scenarios <- unique(df$sim_scenario[df$method == "baseline"])
  }

  have_scenarios <- unique(df$sim_scenario[df$method == m])
  missing_scenarios <- setdiff(target_scenarios, have_scenarios)

  if (length(missing_scenarios) > 0) {
    cat(sprintf("%s is missing %d scenarios\n", m, length(missing_scenarios)))
  }

  for (sc in missing_scenarios) {
    if (!sc %in% names(baseline_split)) {
      warning(sprintf("No baseline row found for scenario %s", sc))
      next
    }

    new_row <- baseline_split[[sc]]
    new_row$method <- m
    new_row$failure_type <- "time"

    rows_to_add[[length(rows_to_add) + 1]] <- new_row
  }
}

if (length(rows_to_add) > 0) {
  rows_to_add_df <- do.call(rbind, rows_to_add)
  df_fixed <- rbind(df, rows_to_add_df)
} else {
  df_fixed <- df
}

# Optional: sort nicely
df_fixed <- df_fixed[order(df_fixed$method, df_fixed$sim_scenario), ]

# Put back
df_testfuncs <- df_fixed



# SURGERY ON REAL DATA
df <- sim_study_realdata$df
df$n_train <- df$n + df$fold_size
df$fold_size <- df$n

front <- c(
  "id", "method", "sim_scenario", "replication",
  "n_train", "fold_size", "input_dim",
  "time_fit", "time_predict", "time"
)

rest <- setdiff(names(df), front)
df <- df[, c(front, rest), drop = FALSE]
df$fold <- NULL
df$n <- NULL

# -----------------------------
# Flags for future subsetting
# -----------------------------
expensive_methods <- c("gp", "hetgp", "deepgp", "rgasp")
df$big_data_flag <- df$n_train >= 1500
df$fast_method_flag <- !(df$method %in% expensive_methods)

# -----------------------------
# Helper objects
# -----------------------------
baseline_df <- df[df$method == "baseline", ]
baseline_key <- paste(baseline_df$id, baseline_df$replication, sep = "___")

all_ids <- sort(unique(df$id))
all_reps <- sort(unique(df$replication))
all_cases <- expand.grid(
  id = all_ids,
  replication = all_reps,
  stringsAsFactors = FALSE
)
all_cases$key <- paste(all_cases$id, all_cases$replication, sep = "___")

# -----------------------------
# 1. Fill missing RVM case with baseline, failure_type = "pred"
# -----------------------------
rvm_df <- df[df$method == "rvm", ]
rvm_key <- paste(rvm_df$id, rvm_df$replication, sep = "___")
missing_rvm <- all_cases[!(all_cases$key %in% rvm_key), , drop = FALSE]

if (nrow(missing_rvm) > 0) {
  add_rvm <- baseline_df[baseline_key %in% missing_rvm$key, , drop = FALSE]
  add_rvm$method <- "rvm"
  add_rvm$failure_type <- "pred"
  add_rvm$fast_method_flag <- TRUE
  add_rvm$big_data_flag <- add_rvm$n_train >= 1500
  df <- rbind(df, add_rvm)
}

# -----------------------------
# 2. Fill missing fast-method cases with baseline, failure_type = "time"
#    (excluding baseline and rvm, which was handled above)
# -----------------------------
fast_methods <- setdiff(sort(unique(df$method[df$fast_method_flag])), c("baseline", "rvm"))

rows_to_add <- list()

for (m in fast_methods) {
  mm <- df[df$method == m, ]
  mm_key <- paste(mm$id, mm$replication, sep = "___")
  missing_mm <- all_cases[!(all_cases$key %in% mm_key), , drop = FALSE]

  if (nrow(missing_mm) == 0) next

  add_mm <- baseline_df[baseline_key %in% missing_mm$key, , drop = FALSE]
  if (nrow(add_mm) == 0) next

  add_mm$method <- m
  add_mm$failure_type <- "time"
  add_mm$fast_method_flag <- TRUE
  add_mm$big_data_flag <- add_mm$n_train >= 1500

  rows_to_add[[m]] <- add_mm
}

if (length(rows_to_add) > 0) {
  df <- rbind(df, do.call(rbind, rows_to_add))
}

df <- df[order(df$method, df$id, df$replication), ]
df_realdata <- df

# SAVE RESULTS
sim_study_testfuncs <- process_sim_study(df_testfuncs, scale_CRPS=FALSE)$df
sim_study_realdata  <- process_sim_study(df_realdata,  scale_CRPS=FALSE)$df

save(sim_study_testfuncs, file="data/sim_study_testfuncs.rda")
save(sim_study_realdata, file="data/sim_study_realdata.rda")




