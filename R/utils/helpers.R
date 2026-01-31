source("R/utils/operators.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(DBI)
  library(RSQLite)
  library(motus)
  library(motusData)
  library(yaml)
})

read_config <- function(path = "config.yml") {
  if (!file.exists(path)) {
    stop("Missing config.yml. Copy config_example.yml to config.yml and edit paths.")
  }
  yaml::read_yaml(path)
}

ensure_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
  invisible(x)
}

as_utc_datetime <- function(x) {
  lubridate::as_datetime(x, tz = "UTC")
}

flatten_alltags <- function(tbl_alltags) {
  tbl_alltags %>%
    mutate(
      recvLat = if_else((is.na(gpsLat) | gpsLat %in% c(0, 999)), recvDeployLat, gpsLat),
      recvLon = if_else((is.na(gpsLon) | gpsLon %in% c(0, 999)), recvDeployLon, gpsLon),
      recvAlt = if_else(is.na(gpsAlt), recvDeployAlt, gpsAlt)
    ) %>%
    select(
      -noise, -slop, -burstSlop, -done, -bootnum, -codeSet, -mfg, -nomFreq,
      -markerNumber, -markerType, -deviceID, -recvDeployAlt, -speciesGroup, -recvAlt
    ) %>%
    collect() %>%
    as.data.frame() %>%
    mutate(
      ts = as_utc_datetime(ts),
      tagDeployStart = as_utc_datetime(tagDeployStart),
      tagDeployEnd   = as_utc_datetime(tagDeployEnd),
      recvLat = plyr::round_any(recvLat, 0.05),
      recvLon = plyr::round_any(recvLon, 0.05),
      recvDeployName = if_else(
        is.na(recvDeployName) | recvDeployName == "",
        paste(recvLat, recvLon, sep=":"),
        recvDeployName
      )
    )
}

get_path_by_run <- function(df, project_id) {
  df %>%
    filter(tagProjID == project_id, !is.na(recvLat), !(recvLat == 0)) %>%
    group_by(
      motusTagID, runID, recvDeployName, ambigID,
      tagDeployLon, tagDeployLat, recvLat, recvLon
    ) %>%
    summarize(
      max_runLen = max(runLen, na.rm = TRUE),
      ts_h = mean(as.numeric(ts), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(ts_h = as_utc_datetime(ts_h)) %>%
    arrange(motusTagID, ts_h)
}
