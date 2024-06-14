# baseline-views-forecasts.R
#
# Patrick T. Brandt
#
# 20240608 : Initial version, sets up provided baseline samples for 
#            forecast comparisons

library(arrow)

##### Read in the provided samples of forecasts ####

boot240 <- rbind(read_parquet("~/VIEWS2/shared_competition_data/benchmarks/boot_240/cm/window=Y2018/bm_boot_240_cm_2018.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/boot_240/cm/window=Y2019/bm_boot_240_cm_2019.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/boot_240/cm/window=Y2020/bm_boot_240_cm_2020.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/boot_240/cm/window=Y2021/bm_boot_240_cm_2021.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/boot_240/cm/window=Y2022/bm_boot_240_cm_2022.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/boot_240/cm/window=Y2023/bm_boot_240_cm_2023.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/boot_240/cm/window=Y2024/bm_boot_240_cm_2024.parquet")
)

# conflictology <- rbind(read_parquet("~/VIEWS2/shared_competition_data/benchmarks/conflictology/cm/window=Y2018/bm_conflictology_cm_2018.parquet"),
#                        read_parquet("~/VIEWS2/shared_competition_data/benchmarks/conflictology/cm/window=Y2019/bm_conflictology_cm_2019.parquet"),
#                        read_parquet("~/VIEWS2/shared_competition_data/benchmarks/conflictology/cm/window=Y2020/bm_conflictology_cm_2020.parquet"),
#                        read_parquet("~/VIEWS2/shared_competition_data/benchmarks/conflictology/cm/window=Y2021/bm_conflictology_cm_2021.parquet"),
#                        read_parquet("~/VIEWS2/shared_competition_data/benchmarks/conflictology/cm/window=Y2022/bm_conflictology_cm_2022.parquet"),
#                        read_parquet("~/VIEWS2/shared_competition_data/benchmarks/conflictology/cm/window=Y2023/bm_conflictology_cm_2023.parquet"),
#                        read_parquet("~/VIEWS2/shared_competition_data/benchmarks/conflictology/cm/window=Y2024/bm_conflictology_cm_2024.parquet")
# )

zero <- rbind(read_parquet("~/VIEWS2/shared_competition_data/benchmarks/zero/cm/window=Y2018/bm_zero_cm_2018.parquet"),
              read_parquet("~/VIEWS2/shared_competition_data/benchmarks/zero/cm/window=Y2019/bm_zero_cm_2019.parquet"),
              read_parquet("~/VIEWS2/shared_competition_data/benchmarks/zero/cm/window=Y2020/bm_zero_cm_2020.parquet"),
              read_parquet("~/VIEWS2/shared_competition_data/benchmarks/zero/cm/window=Y2021/bm_zero_cm_2021.parquet"),
              read_parquet("~/VIEWS2/shared_competition_data/benchmarks/zero/cm/window=Y2022/bm_zero_cm_2022.parquet"),
              read_parquet("~/VIEWS2/shared_competition_data/benchmarks/zero/cm/window=Y2023/bm_zero_cm_2023.parquet"),
              read_parquet("~/VIEWS2/shared_competition_data/benchmarks/zero/cm/window=Y2024/bm_zero_cm_2024.parquet")
)

lastobs <- rbind(read_parquet("~/VIEWS2/shared_competition_data/benchmarks/last/cm/window=Y2018/bm_last_cm_2018.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/last/cm/window=Y2019/bm_last_cm_2019.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/last/cm/window=Y2020/bm_last_cm_2020.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/last/cm/window=Y2021/bm_last_cm_2021.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/last/cm/window=Y2022/bm_last_cm_2022.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/last/cm/window=Y2023/bm_last_cm_2023.parquet"),
                 read_parquet("~/VIEWS2/shared_competition_data/benchmarks/last/cm/window=Y2024/bm_last_cm_2024.parquet")
)

# Read in the actuals for comparison
a18 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2018/cm_actuals_2018.parquet")
a19 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2019/cm_actuals_2019.parquet")
a20 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2020/cm_actuals_2020.parquet")
a21 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2021/cm_actuals_2021.parquet")
a22 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2022/cm_actuals_2022.parquet")
a23 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2023/cm_actuals_2023.parquet")
a24 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2024/cm_actuals_2024.parquet")

a18$year <- 2018; a19$year <- 2019; a20$year <- 2020; a21$year <- 2021; 
a22$year <- 2022; a23$year <- 2023; a24$year <- 2024;   

a18$month <- a18$month - min(a18$month) + 1
a19$month <- a19$month - min(a19$month) + 1
a20$month <- a20$month - min(a20$month) + 1
a21$month <- a21$month - min(a21$month) + 1
a22$month <- a22$month - min(a22$month) + 1
a23$month <- a23$month - min(a23$month) + 1
a24$month <- a24$month - min(a24$month) + 1

actuals <- rbind(a18,a19,a20,a21,a22,a23,a24)
actuals$date <- ISOdate(actuals$year, actuals$month, 1)

rm(a18,a19,a20,a21,a22,a23,a24)

# Label information for countries
country.labels <- read.csv("~/VIEWS2/matching_tables/countries.csv")[,c(1,2,5)]
country.labels$country_id <- country.labels$id

#### Combine forecasts together with the actuals ####

# Standardize column names for scoringutils and merging
colnames(actuals) <- c("observed", "month_id", "country_id", "year", "month", "date")

colnames(boot240) <- c("predicted", "month_id", "country_id", "sample_id")
boot240$model <- "Bootstrap"
boot240$sample_id <- as.integer(boot240$sample_id + 1)

colnames(lastobs) <- c("predicted", "month_id", "country_id", "sample_id")
lastobs$model <- "Constituent Poisson"
lastobs$sample_id <- as.integer(lastobs$sample_id + 1)

colnames(zero) <- c("predicted", "month_id", "country_id", "sample_id")
zero$model <- "Exactly Zero"
zero$sample_id <- as.integer(zero$sample_id + 1)

# Add country labels to the actuals
actuals <- merge(actuals, country.labels[,2:4])
colnames(actuals) <- c(colnames(actuals)[1:6], "location", "isoab")

# Now merge forecasts and format for scoringutils
boot240 <- merge(boot240, actuals, by.x = c("month_id", "country_id"),
                 by.y = c("month_id", "country_id"), all.x=TRUE)
lastobs <- merge(lastobs, actuals, by.x = c("month_id", "country_id"),
                 by.y = c("month_id", "country_id"), all.x=TRUE)
zero <- merge(zero, actuals, by.x = c("month_id", "country_id"),
                 by.y = c("month_id", "country_id"), all.x=TRUE)

rm(actuals)

# Now make the final object we need and save it

viewsforcs <- rbind(boot240,lastobs,zero)
rm(boot240,lastobs,zero)
gc()

# Reorder the columns so they match other data
viewsforcs <- viewsforcs[,c(2,1,3,4,5,6,7,8,10,11,9)]

save(viewsforcs, file="viewsforcs.RData")

