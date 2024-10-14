########################################################################
#Chapter 1
(n_wd <- getwd())  # Confirmation of working directory

# Table 1.1 Measurement Results (Days)
# Antigen Test Data
x <- c(
    11.5, 14.0, 15.0, 10.0, 14.0, 14.5, 12.5, 12.5, 12.5, 12.0,
    13.0, 12.0, 12.5, 13.5, 15.0, 13.5, 12.5, 12.0, 12.0, 17.0
)

# Table 1.2 Frequency Distribution of "Antigen Test" (Class Width: 1 Day)
(xs10 <- seq(9.5, 17.5, 1))  # Class intervals
(xt <- table(cut(x, xs10, right = F)))  # Frequency
(xp <- xt / length(x))  # Probability
cumsum(xt)  # Cumulative frequency
cumsum(xp)  # Cumulative probability

# Figure 1.1 Histogram - Left with Class Width 1 Day, Right with Class Width 0.5 Day
par(mfrow = c(1, 2))
xs05 <- seq(9.75, 17.25, 0.5)
hist(x, xs10, xlab = "Days", ylab = "", col = 4, right = F, main = "", cex.lab = 2, cex.axis = 2)
hist(x, xs05, xlab = "Days", ylab = "", col = 2, right = F, main = "", cex.lab = 2, cex.axis = 2)
par(mfrow = c(1, 1))
#dev.copy2pdf(file="../English_tex/fig/chap01/fig1_1.pdf")

# Summary Statistics for "Antigen Test" (Product Moment) 1.1.5  1.1.6
mean(x)  # Mean
round(s2 <- var(x), 3)  # Variance with denominator n
round(s <- sqrt(s2), 3)  # Standard Deviation

# Summary Statistics for "Antigen Test" (Quantiles) 1.1.7  1.1.8  1.1.9
sort(x)  # Sorted in ascending order
median(x)  # Median
quantile(x, 0.3, type = 1)  # 30th Percentile
quantile(x, 0.7, type = 1)  # 70th Percentile (Mode refers to the histogram)

# Probability Distribution Function 1.2.7
round(pnorm(12.5, 13.075, 1.502), 3)

# Probability of Observing Data in an Arbitrary Interval 1.2.9
round(pnorm(12.0, 13.075, 1.502), 3)
round(pnorm(10.0, 13.075, 1.502), 3)
round(pnorm(12.0, 13.075, 1.502) - pnorm(10.0, 13.075, 1.502), 3)

# 95% Prediction Interval  1.2.10
round(13.075 - 1.96 * 1.502, 3)
round(13.075 + 1.96 * 1.502, 3)

########################################################################
#Chapter 1 Practical Assignment
# Table 1.3 Measurement Results (Seconds) - Input Data for "Perceived Time"
y <- c(
    31.43, 31.09, 33.38, 30.49, 29.62, 35.40, 32.58, 28.96, 29.43, 28.52,
    25.39, 32.68, 30.51, 30.15, 32.33, 30.43, 32.50, 32.07, 32.35, 31.57
)

# Frequency Distribution Table
(ys50 <- seq(24.5, 35.5, 1))  # Class intervals
(yt <- table(cut(y, ys50, right = F)))  # Frequency
(yp <- yt / length(y))  # Probability
cumsum(yt)  # Cumulative frequency
cumsum(yp)  # Cumulative probability

# Histograms
par(mfrow = c(1, 2))
ys50 <- seq(24.5, 36.5, 1)
hist(y, ys50, xlab = "Time (Seconds)", col = 4, right = F)
ys25 <- seq(24, 36, 2)
hist(y, ys25, xlab = "Time (Seconds)", col = 2,right =F)
# "Perceived Time" Summary Statistics (Product Moment)
mean(y)  # Mean
round(ys2 <- var(y), 3)  # Variance with denominator n
round(ys <- sqrt(ys2), 3)  # Standard Deviation

# "Perceived Time" Summary Statistics (Quantiles)
sort(y)  # Sorted in ascending order
median(y)  # Median
quantile(y, 0.3, type = 1)  # 30th Percentile
quantile(y, 0.7, type = 1)  # 70th Percentile

########################################################################
# Chapter 1 Practical Problems

final_exam <- scan("dat/FinalExam.dat", sep = "")
income_men <- scan("dat/incomeMen.dat", sep = "")
geyser_wait_time <- scan("dat/Geyser.dat", sep = "")

par(mfrow = c(3, 1))
hist(final_exam, xlim = c(40, 100), breaks = 30, main = "", xlab = "Test Scores of 'Psychological Statistics' Final Exam", ylab = "", cex.axis = 2.0, cex.lab = 2.0)
hist(income_men, xlim = c(0, 2000), breaks = 30, main = "", xlab = "Annual Income of Men in the 2014", ylab = "", cex.axis = 2.0, cex.lab = 2.0)
hist(geyser_wait_time, breaks = 30, main = "", xlab = "Intermittent geyser Waiting Time (minutes)", ylab = "", cex.axis = 2.0, cex.lab = 2.0)
par(mfrow = c(1, 1))

# Self-Study
# 1. What are the mean and median of the final exam scores? The mode is 97 points.
# 2. What are the mean and median of men's annual income? Find the 90th, 95th, and 99th percentiles. The mode is 3.5 million yen.
# 3. What are the mean and median of the intermittent geyser waiting times? (The mode is 53 minutes, 77 minutes.)


########################################################################

