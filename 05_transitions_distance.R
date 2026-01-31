source("R/utils/helpers.R")

cfg <- read_config()
ensure_dir("outputs/transitions")

project_id <- cfg$project_id
motus_dir  <- cfg$motus_dir

motus_file <- file.path(motus_dir, paste0("project-", project_id, ".motus"))
sql_motus <- DBI::dbConnect(RSQLite::SQLite(), motus_file)

tbl_alltags <- dplyr::tbl(sql_motus, "alltags")

tag_keep <- cfg$tag_ids_keep
if (is.null(tag_keep) || length(tag_keep) == 0) {
  stop("Set tag_ids_keep in config.yml.")
}

for (tag_id in tag_keep) {
  trans <- motus::siteTrans(
    dplyr::filter(tbl_alltags, motusTagID == tag_id),
    latCoord = "recvDeployLat",
    lonCoord = "recvDeployLon"
  )

  trans <- trans %>%
    dplyr::filter(rate <= 6.8) %>%
    dplyr::select(motusTagID, ts.x, ts.y, lat.x, lon.x, lat.y, lon.y, dist, rate, bearing)

  out <- file.path("outputs/transitions", paste0("transitions_", tag_id, ".csv"))
  readr::write_csv(trans, out)
}

DBI::dbDisconnect(sql_motus)
