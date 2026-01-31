## Notes on customizing for your Motus project

### Blocking rules
Put your project-specific false-detection rules in `scripts/03_blocking_rules.R`
and control them via `config.yml`.

If you want blocks by `ambigID`, add:

```r
tmp <- df_keep %>%
  filter(ambigID %in% c(-359, -358)) %>%
  select(motusTagID, runID) %>%
  distinct() %>%
  mutate(probability = 0)

df_block_all <- bind_rows(df_block_all, tmp)
```

### Mapping
- Outline maps require no API keys.
- `ggmap` needs a Google API key stored in `config.yml` (git-ignored).
