########################################################################
#Chapter 2
(getwd())  #Confirmation of working directory

#Figure 2.1: Posterior Distribution of the Probability of Pardon
set.seed(1234)                                  # Random seed
nod <- 1000000                                  # Number of random samples
BdAa <- runif(nod)                              # Setting a uniform distribution for f(Bd|Aa) in 2.4.5
AaBd <- BdAa / (BdAa + 1)                       # Probability that Prisoner A is pardoned (Equation 2.34)
mean(AaBd)                                      # Estimated Expected A Posteriori (EAP)
median(AaBd)                                    # Estimated Median (MED)
hist(AaBd, breaks = 50, freq = FALSE, col = 3)  # Posterior distribution, MAP estimate = 0.5

########################################################################
#Chapter 2, Section 2.6, Exercise
#(For data with mean 31.044 and standard deviation 2.068 as shown in Table 1.3)

##1
round(dnorm(28, 31.044, 2.068), 4)
round(dnorm(32, 31.044, 2.068), 4)

##2
1 - round(pnorm(33, 31.044, 2.068), 4)

##3
round(pnorm(34, 31.044, 2.068) - pnorm(28, 31.044, 2.068), 4)

##4
round(31.044 - 1.96 * 2.068, 2)
round(31.044 + 1.96 * 2.068, 2)

##5
round(qnorm(0.95, 31.044, 2.068), 4)
round(qnorm(0.025, 31.044, 2.068), 4)

##6
round(qnorm(0.25, 31.044, 2.068), 2)
round(qnorm(0.5, 31.044, 2.068), 2)
round(qnorm(0.75, 31.044, 2.068), 2)

########################################################################
#Chapter 2, Section 2.7, Practical assignments
(0.15 * 0.5) / ((0.15 * 0.5) + (0.02 * (1 - 0.5)))
(0.22 * 0.8823529) / ((0.22 * 0.8823529) + (0.01 * (1 - 0.8823529)))
(0.25 * 0.9939759) / ((0.25 * 0.9939759) + (0.01 * (1 - 0.9939759)))
