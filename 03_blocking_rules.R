source("R/utils/helpers.R")

cfg <- read_config()

block_cfg <- cfg$blocking
if (is.null(block_cfg) || !isTRUE(block_cfg$enabled)) {
  message("Blocking disabled. Skipping.")
  quit(save = "no")
}

df_keep <- readRDS("data_intermediate/df_alltags_keep.rds")
df_block_all <- readRDS("data_intermediate/df_alltags_excl.rds") %>%
  dplyr::mutate(probability = 0)

month_blocks <- block_cfg$month_blocks
if (!is.null(month_blocks) && length(month_blocks) > 0) {
  for (rule in month_blocks) {
    months <- rule$months
    tag_ids <- rule$tag_ids
    if (is.null(months) || is.null(tag_ids)) next

    tmp <- df_keep %>%
      dplyr::filter(lubridate::month(ts) %in% months, motusTagID %in% tag_ids) %>%
      dplyr::select(motusTagID, runID) %>%
      dplyr::distinct() %>%
      dplyr::mutate(probability = 0)

    df_block_all <- dplyr::bind_rows(df_block_all, tmp)
  }
}

date_blocks <- block_cfg$date_blocks
if (!is.null(date_blocks) && length(date_blocks) > 0) {
  for (rule in date_blocks) {
    start <- as.Date(rule$start)
    end   <- as.Date(rule$end)
    tag_ids <- rule$tag_ids
    if (is.na(start) || is.na(end) || is.null(tag_ids)) next

    tmp <- df_keep %>%
      dplyr::filter(as.Date(ts) >= start, as.Date(ts) <= end, motusTagID %in% tag_ids) %>%
      dplyr::select(motusTagID, runID) %>%
      dplyr::distinct() %>%
      dplyr::mutate(probability = 0)

    df_block_all <- dplyr::bind_rows(df_block_all, tmp)
  }
}

df_block_all <- df_block_all %>% dplyr::distinct(motusTagID, runID, .keep_all = TRUE)
df_blocked <- dplyr::anti_join(df_keep, df_block_all, by = c("motusTagID", "runID"))

saveRDS(df_blocked, "data_intermediate/df_alltags_blocked.rds")
message("Saved data_intermediate/df_alltags_blocked.rds")
