source("R/utils/helpers.R")

cfg <- read_config()
Sys.setenv(TZ = "UTC")

project_id <- cfg$project_id
motus_dir  <- cfg$motus_dir

ensure_dir(motus_dir)
ensure_dir("data_intermediate")

message("Updating Motus DB for project_id = ", project_id)
message("Directory: ", motus_dir)

sql_motus <- tagme(
  projRecv = project_id,
  new      = TRUE,
  update   = TRUE,
  dir      = motus_dir
)

saveRDS(list(project_id = project_id, motus_dir = motus_dir),
        file = "data_intermediate/motus_project_info.rds")

message("Done. Next: scripts/02_activity_filter_and_flatten.R")
