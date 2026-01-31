source("R/utils/helpers.R")

cfg <- read_config()
Sys.setenv(TZ = "UTC")
ensure_dir("data_intermediate")

project_id <- cfg$project_id
motus_dir  <- cfg$motus_dir

motus_file <- file.path(motus_dir, paste0("project-", project_id, ".motus"))
if (!file.exists(motus_file)) {
  stop("Missing .motus file at: ", motus_file,
       "\nRun scripts/01_download_update_motus_db.R or fix motus_dir in config.yml.")
}

sql_motus <- DBI::dbConnect(RSQLite::SQLite(), motus_file)

af <- cfg$activity_filter %||% list()

tbl_alltags <- motus::filterByActivity(
  sql_motus,
  minLen  = af$minLen  %||% 3,
  maxLen  = af$maxLen  %||% 4,
  maxRuns = af$maxRuns %||% 1000,
  ratio   = af$ratio   %||% 0.99,
  return  = "all"
)

start_date <- cfg$start_date_utc
if (!is.null(start_date) && nzchar(start_date)) {
  start_ts <- as.numeric(as.POSIXct(start_date, tz = "UTC"))
  tbl_alltags <- tbl_alltags %>% dplyr::filter(ts > start_ts)
}

df_alltags <- flatten_alltags(tbl_alltags)

df_keep <- dplyr::filter(df_alltags, probability == 1)
df_excl <- dplyr::filter(df_alltags, probability == 0) %>%
  dplyr::select(motusTagID, runID) %>%
  dplyr::distinct()

saveRDS(df_alltags, "data_intermediate/df_alltags_flat.rds")
saveRDS(df_keep,   "data_intermediate/df_alltags_keep.rds")
saveRDS(df_excl,   "data_intermediate/df_alltags_excl.rds")

DBI::dbDisconnect(sql_motus)

message("Saved. Next: scripts/03_blocking_rules.R (optional)")
