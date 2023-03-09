# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2022-10-14
#
# Script Description: prepare data for inference with UnitEvents

# 1. Import libraries --------------------------------------------------
# Load functions from main folder
devtools::load_all()

packages <- c("here", "dplyr", "lubridate", "ggplot2")
base::lapply(packages, require)

# 2. Parameters --------------------------------------------------------
thr_obs <- 1 # Discard cameras that have less than thr_obs captures
thr_freq <- 1/30 # Discard cameras for which the mean frequency is less than thr_freq (day⁻¹)

scale <- 10e6

outputs_path <- here("outputs/05_example_real_data")

# 3. Read data ---------------------------------------------------------
dat <- read.csv(here("data/camtrap_data/data.csv"),
                stringsAsFactors = FALSE)

# 4. Format data and add columns ---------------------------------------

## 4.1 Add datetime ---------------------------------------------------------
dat <- dat %>%
  mutate(datetime = paste(eventDate, eventTime)) %>%
  select(-c(eventDate, eventTime))

## 4.2 Parse date time -----------------------------------------------------
fmt <- "%Y-%m-%d %H:%M:%S"
dat$datetime <- as.POSIXct(dat$datetime, format = fmt)
dat$datetime <- force_tz(dat$datetime, tz = "Etc/GMT-2")

## 4.3 Add stamps ----------------------------------------------------------
dat <- add_stamps(dat)
dat$row_ID <- 1:nrow(dat)

# 5. Data quality ----------------------------------------------------------
dat <- filter_inactive_cameras(df = dat, 
                               thr_obs = thr_obs, 
                               thr_freq = thr_freq)

# 6. Shift duplicates --------------------------------------------------

# Shift pictures taken at the same time
dat <- shift_duplicates(dat, precision = scale)

# 7. Data summary --------------------------------------------------
unique(dat$snapshotName)

# Number of occurrence events (should be the same)
nrow(dat)

dup <- sum(duplicated(dat[, c("cameraID","stamp")])) 
(nocc <- nrow(dat) - dup)

# Number of pictures (same datetime and same camera)
dup <- sum(duplicated(dat[, c("cameraID","datetime")])) 
(npic <- nrow(dat) - dup)

# Dates
range(dat$datetime)

# Camera activity
df_summ <- get_sampling_info(dat, return = TRUE)
min(df_summ$duration)
max(df_summ$duration)
sd(df_summ$duration)
nrow(df_summ)

# Species abundance
ggplot(dat) + geom_bar(aes(x = snapshotName))
dat %>% 
  group_by(snapshotName) %>%
  summarize(n = n())

# Points
ggplot(dat) +
  geom_point(aes(x = datetime, y = cameraID, 
                 col = snapshotName), show.legend = FALSE) +
  facet_wrap(facets = ~ locationID, scales = "free")

# 8. Write data --------------------------------------------------------
# Remove count (useless)
dat <- dat %>% select(-count)

write.csv(dat,
          file.path(outputs_path, "cleaned_data.csv"),
          row.names = FALSE)
