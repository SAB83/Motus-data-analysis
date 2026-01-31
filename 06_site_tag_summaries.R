source("R/utils/helpers.R")

cfg <- read_config()
ensure_dir("outputs/summaries")
ensure_dir("outputs/maps")

df_in <- if (file.exists("data_intermediate/df_alltags_blocked.rds")) {
  readRDS("data_intermediate/df_alltags_blocked.rds")
} else {
  readRDS("data_intermediate/df_alltags_keep.rds")
}

df_in <- df_in %>%
  dplyr::mutate(
    ts = lubridate::as_datetime(ts, tz = "UTC"),
    year = lubridate::year(ts),
    doy  = lubridate::yday(ts)
  ) %>%
  dplyr::filter(!is.na(recvLat))

readr::write_csv(motus::simSiteDet(df_in), "outputs/summaries/simSiteDet.csv")
readr::write_csv(motus::siteSum(df_in, units = "mins"), "outputs/summaries/siteSum.csv")
readr::write_csv(motus::siteSumDaily(df_in, units = "mins"), "outputs/summaries/siteSumDaily.csv")
readr::write_csv(motus::tagSum(df_in), "outputs/summaries/tagSum.csv")

tag_keep <- cfg$tag_ids_keep
if (!is.null(tag_keep) && length(tag_keep) > 0) {
  readr::write_csv(motus::tagSumSite(dplyr::filter(df_in, motusTagID %in% tag_keep)),
                   "outputs/summaries/tagSumSite_selected.csv")

  tag_id <- tag_keep[[1]]
  df_site <- dplyr::filter(df_in, motusTagID == tag_id)

  p_sig <- ggplot(df_site, aes(ts, sig,
                              colour = paste(recvDeployLat, recvDeployLon, recvDeployName, sep=": "))) +
    geom_point() +
    theme_bw() +
    labs(x = "Time (UTC)", y = "Signal strength", colour = "Site",
         title = paste("Signal strength over time | motusTagID", tag_id)) +
    facet_wrap(~ as.Date(ts), scales = "free_x")

  ggsave(paste0("outputs/maps/signal_by_day_", tag_id, ".png"),
         p_sig, width = 12, height = 8, dpi = 300)
}
