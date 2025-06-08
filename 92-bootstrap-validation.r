# 2024 05 bootstrap check of the standard normal null
# nboot = 1000 takes about 4 minutes on 4 cores

# Setup -----------------------------------------------------------------------
rm(list = ls())

# Set working directory to unbreakable-bh folder
if (basename(getwd()) != "unbreakable-bh") {
  # Try to find unbreakable-bh directory
  if (dir.exists("unbreakable-bh")) {
    setwd("unbreakable-bh")
  } else if (dir.exists("../unbreakable-bh")) {
    setwd("../unbreakable-bh")  
  } else {
    stop("Please run this script from the unbreakable-bh directory or its parent directory.")
  }
}

source("0-functions.r")
load("../data/emp_data.Rdata")

# load bootstrap of null model
load("../data/bootnull.Rdata")
bootnull[, tstat:=mean/vol*sqrt(nmonth)]

# summarize, add standard normal ---------------
histedge = c(-Inf, seq(-8,8,0.25), Inf)
boothist = bootnull %>%
    # histogram counts within each bootstrap
    mutate(bin = cut(abs(tstat), breaks = histedge, include.lowest = TRUE)) %>%
    group_by(booti, bin) %>%
    summarize(count = n(), .groups = 'drop')  %>% 
    # normalize by total signals
    left_join(bootnull %>% group_by(booti) %>% summarize(ntotal=n()), by = "booti") %>% 
    mutate(freq=count/ntotal)  %>% 
    # find bin midpoints 
    mutate(bin = str_remove_all(bin, "[\\(\\)\\[\\]]")
        , left = str_split(bin, ",")
        , left = sapply(left, function(x) as.numeric(x[1]))
        , right = str_split(bin, ",")
        , right = sapply(right, function(x) as.numeric(x[2]))
        , mid = (left + right) / 2
    )

# summarize
histdat <- boothist %>%
    group_by(left, mid, right) %>%
    summarize(
        f_mean = mean(freq), f_low = quantile(freq, 0.05), f_high = quantile(freq, 0.95),
        f_sd = sd(freq)
    ) %>%
    ungroup() 

# add normal
histdat <- histdat %>%
    mutate(fnorm =
        (mid > 0) * 2 * (pnorm(right) - pnorm(left))
    ) %>%
    setDT()

# plot simpler -----------------------------------------------------
lablist <- c("Bootstrapped Null", "Normal(0,1)")

p <- histdat[mid > 0] %>%
    select(mid, f_mean, fnorm) %>%
    pivot_longer(cols = c("f_mean", "fnorm")) %>%
    ggplot(aes(x = mid, y = value), size = 2) +
    geom_line(aes(color = name), size = 0.5) +
    geom_point(aes(color = name, shape = name), size = 3) +
    labs(x = "Absolute t-statistic", y = "Probability") +
    scale_x_continuous(breaks = seq(-8, 20, 1)) +
    coord_cartesian(xlim = c(0, 4)) +
    theme(
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.1, "cm"),
        legend.position = c(80, 80) / 100,
        legend.key.width = unit(1, "cm"),
        legend.spacing.y = unit(0.000001, "cm"),
        legend.background = element_rect(colour = "black", fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")
    ) +
    scale_color_manual(values = c(MATBLUE, MATRED), labels = lablist) +
    scale_shape_manual(values = c(1, 4), labels = lablist)

ggsave("../results/boot-null-tabs.pdf", p,
    width = 5, height = 2.5, scale = 1,
    device = cairo_pdf
)

