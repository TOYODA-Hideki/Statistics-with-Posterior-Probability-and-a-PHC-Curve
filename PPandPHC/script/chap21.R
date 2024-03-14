############ Chapter 21: Advanced Experimental Design
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
   library(cmdstanr)           #load and attach add-on package cmdstanr
   library(posterior)          #load and attach add-on package posterior
   library(rstan)              #load and attach add-on package rstan
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) # Defining probability points

################### 21.1 Inference of One-Way Matched Design

###### 21.1.1 Sensory Experiment of Two-Point Threshold

###### Table 21.1 Data of Two-Point Threshold Sensory Experiment
(TwoPointThreshold<-read.csv("./scrT/dat/touch_threshold.csv", header = TRUE, encoding="utf-8"))

###### Figure 21.1 Measurement Values of Two-Point Threshold (mm per body part)
boxplot(TwoPointThreshold$Length~TwoPointThreshold$BodyPart,xlab="Body Part",ylab="",cex.axis=2,cex.lab=2)
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_1.pdf")

###### Mean & sd of Two-Point Threshold Measurement Values per Body Part
tapply(TwoPointThreshold$Length,TwoPointThreshold$BodyPart,mean)
tapply(TwoPointThreshold$Length,TwoPointThreshold$BodyPart,sd)

###### Figure 21.2 Measurement Values of Two-Point Threshold (mm per subject)
boxplot(TwoPointThreshold$Length~TwoPointThreshold$Subject,xlab="Subject",cex.axis=1.4,cex.lab=1.5)
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_2.pdf")

###### Mean & sd of Two-Point Threshold Measurement Values per Subject, bottom, specific numbers not given
tapply(TwoPointThreshold$Length,TwoPointThreshold$Subject,mean)
tapply(TwoPointThreshold$Length,TwoPointThreshold$Subject,sd)

####################### Inference of One-Way Matched Design STAN Code
RBlockD <- '
data { 
  int<lower=0>  n;                     // Total number of data
  int<lower=0>  J;                     // Number of groups
  int<lower=0>  K;                     // Number of blocks
  vector[n]     y;                     // Criterion variable
  int<lower=0>  j[n];                  // Classification variable A
  int<lower=0>  k[n];                  // Block variable B
}
parameters {
  vector[J]   mu;                      // Effect of each group
  vector[K]   s;                       // Effect of block
  real<lower=0> s_e;                   // Error SD
  real<lower=0> s_s;                   // Block SD
}
model {
    y ~ normal(mu[j] + s[k], s_e);     // Equation (21.1), (21.2)
    s ~ normal(0, s_s);                // Equation (21.3)
}
generated quantities{
  real tmu;                            // Overall mean
  real s_mu;                           // Factor SD
  tmu = mean(mu);                      // Mean of μj
  s_mu = sqrt(variance(mu) * (J - 1) / J); // SD of μj
}
';
par<-c("mu","s_e","s_s","tmu","s_mu","s");                      # Parameters
dataStan <-list(n=length(TwoPointThreshold$Length),J=max(TwoPointThreshold$BodyPartIndex), # Input
                K=max(TwoPointThreshold$Subject),y=TwoPointThreshold$Length,j=TwoPointThreshold$BodyPartIndex,k=TwoPointThreshold$Subject)



########### Execution with cmdstanr
modfileRandom <- write_stan_file(RBlockD)             # Writing temporary file
modRandom <- cmdstan_model(modfileRandom)             # Compilation
csrfitRandom <- modRandom$sample(data = dataSetRandom, chains = 5, iter_sampling = 20000,
                                 iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
fitRandom <- rstan::read_stan_csv(csrfitRandom$output_files()) # Conversion to stan format

########### Execution with rstan
#fitRandom <- stan(model_code = RBlockD, data=dataSetRandom,
#  pars=par, seed=1234, chains=5, warmup=1000, iter=21000)
#save(fitRandom, file="./scrT/obje/stan_obje0801")
#load(file="./scrT/obje/stan_obje0801"); # Loading the object created beforehand with stan()

extRandom <- extract(fitRandom, par)

###### Table 21.2 Posterior distribution of parameters in Randomized Block Design
eta2 = (extRandom$s_mu^2)/((extRandom$s_mu^2)+(extRandom$s_e^2)); # Explanation rate (Equation 21.5)
gqcal(extRandom$mu[,1])
gqcal(extRandom$mu[,2])
gqcal(extRandom$mu[,3])
print('------------------------------------------')
gqcal(extRandom$s_e)
gqcal(extRandom$s_s)
print('------------------------------------------')
gqcal(extRandom$tmu)
gqcal(extRandom$s_mu)
print('------------------------------------------')
gqcal(extRandom$s[,4])
print('------------------------------------------')
gqcal(eta2)

###### Probability of difference between groups (verification of necessary condition for significance), bl.2
mean(extRandom$mu[,1] > extRandom$mu[,2])
mean(extRandom$mu[,2] > extRandom$mu[,3])
mean(extRandom$mu[,1] > extRandom$mu[,3])

###### Figure 21.3 Curve for difference in mean values using PHC
par(mfrow=c(2,1));                                    # Figure 21.3
phc01(seq01=seq(0,1.2,0.05), a=extRandom$mu[,1], b=extRandom$mu[,2],
      xlab="c<???��?���?��?��?���?��?���?��?���?��?���?��?��?���?��?���?��?��?���?��?���?��?��2")
phc01(seq01=seq(0,1.2,0.05), a=extRandom$mu[,2], b=extRandom$mu[,3],
      xlab="c<???��?���?��?��?���?��?���?��?���?��?���?��?��?���?��?���?��?��?���?��?���?��?��3")
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_3.pdf")

###### Probability of substantial difference  l.3
phc01(0.7, a=extRandom$mu[,1], b=extRandom$mu[,2], byoga="no")
phc01(0.5, a=extRandom$mu[,2], b=extRandom$mu[,3], byoga="no")

############# 21.2 Estimation of One Nested Factor
####### 21.2.1 Latent Learning in a Maze Task for Mice

####### Table 21.3 Data of Latent Learning in Maze Task for Mice
(LatentLearning <- read.csv("./scrT/dat/latent_learning.csv", header = TRUE))

###### Figure 21.4 Latent Learning in Maze Task for Mice (Per Condition)
boxplot(LatentLearning$Time ? LatentLearning$Condition, xlab="Condition", cex.lab=1.2, cex.axis=1.8)
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_4.pdf")

###### Figure 21.5 Latent Learning in Maze Task for Mice (Per Mouse)
boxplot(LatentLearning$Time ? LatentLearning$mouse, xlab="Mouse", ylab="")
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_5.pdf")

###### Mean and sd Per Condition (No specific numbers given in textbook)
tapply(LatentLearning$Time, LatentLearning$Condition, mean)
tapply(LatentLearning$Time, LatentLearning$Condition, sd)

############################ Estimation of One Nested Factor STAN Code
nestedD <- '
data { 
  int<lower=0>  n;                     // Total number of data
  int<lower=0>  J;                     // Number of groups
  int<lower=0>  K;                     // Number of blocks
  vector[n]     y;                     // Criterion variable
  int<lower=0>  j[n];                  // Classification variable A
  int<lower=0>  k[n];                  // Block variable B
}
parameters {
  vector[J]     mu;                    // Effect of each group
  vector[J*K]   s;                     // Effect of block
  real<lower=0> s_e;                   // Error SD
  real<lower=0> s_s;                   // Block SD
}
model {
    y ~ normal(mu[j] + s[k], s_e);     // Equation (21.6), (21.7)
    s ~ normal(0, s_s);                // Equation (21.8)
}
';

par <- c("mu", "s_e", "s_s")                                        # Parameters
dataSetNest <- list(n=length(LatentLearning$Time), J=max(LatentLearning$ConditionIndex), # Input
                    K=max(LatentLearning$Rat), y=LatentLearning$Time, j=LatentLearning$ConditionIndex, k=LatentLearning$Rat)

########### Execution with cmdstanr
modfileNest <- write_stan_file(nestedD)             # Writing temporary file
modNest <- cmdstan_model(modfileNest)               # Compilation
csrfitNest <- modNest$sample(data=dataSetNest, chains=5, iter_sampling=20000,
                             iter_warmup=1000, parallel_chains=5, seed=1234)  # MCMC
fitNest <- rstan::read_stan_csv(csrfitNest$output_files()) # Conversion to stan format

########### Execution with rstan
#fitNest <- stan(model_code=nestedD, data=dataSetNest,
#                pars=par, seed=1234, chains=5, warmup=1000, iter=21000)
#save(fitNest, file="./scrT/obje/stan_obje0802")
#load(file="./scrT/obje/stan_obje0802"); # Loading the object created beforehand with stan()

extNest <- extract(fitNest, par);               # Extracting random numbers

###### Table 21.4 Summary of Posterior Distribution of Parameters and Quantities for Branching Plan (First Part)
gqcal(extNest$mu[,1])
gqcal(extNest$mu[,2])
gqcal(extNest$mu[,3])
gqcal(extNest$s_e)
gqcal(extNest$s_s)

###### Table 21.4 (Second Part) Posterior Distribution of Mean Differences
def32 = extNest$mu[,3] - extNest$mu[,2];    # (Equation 21.10)
def13 = extNest$mu[,1] - extNest$mu[,3];    # (Equation 21.11)
gqcal(def32);
gqcal(def13);                               # End of Table 21.4


##################### 21.3 Estimation of Two-Factor Design with Correspondence
###### 21.3.1 Comparative Study of Persuasive Messages

###### Table 21.5 Data of Comparative Study of Persuasive Messages
(PersuasiveMessages <- read.csv("./scrT/dat/persuasive_messages.csv", header = TRUE))

###### Figure 21.6 Boxplot for Each Cell
boxplot(PersuasiveMessages$Rating ? PersuasiveMessages$Message + PersuasiveMessages$Font, xlab="Message/Font", ylab="Rating")
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_6.pdf")

###### Figure 21.7 Boxplot for Each Subject
boxplot(PersuasiveMessages$Rating ? PersuasiveMessages$Subject, xlab="Subject", ylab="Rating")
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_7.pdf")

####################################### Estimation of Two-Factor Design with Correspondence STAN Code
RBFD <- '
functions{
  vector zero_sum_vector(int a, vector m1A){// Main effect zero-sum
    vector[a] muA;            
    for(i in 1:(a-1)){ muA[i] = m1A[i];}   // Copy up to a-1 as it is
    muA[a] = -sum(m1A);                    // For the ath, insert -1 times the sum
    return(muA);
  }
  matrix zero_sum_matrix(int a, int b, matrix m1AB){// Interaction effect zero-sum
    vector [a-1] m1a;                       // Vector for row sum
    vector [b-1] m1b;                       // Vector for column sum
    matrix [a,b] muAB;                      // Matrix to be returned
    for(i in 1:(a-1)){ m1a[i] = 0.0;}      // Zero clear
    for(j in 1:(b-1)){ m1b[j] = 0.0;}      // Zero clear
    for(i in 1:(a-1)){ for(j in 1:(b-1)){muAB[i,j] = m1AB[i,j];}} // Copy as it is
    for(i in 1:(a-1)){ for(j in 1:(b-1)){m1a[i] += m1AB[i,j];}} // Create row sum
    for(j in 1:(b-1)){ for(i in 1:(a-1)){m1b[j] += m1AB[i,j];}} // Create column sum
    for(i in 1:(a-1)){muAB[i,b] = (-1)*m1a[i]; }    // Insert -1 times the row sum into bth column
    for(j in 1:(b-1)){muAB[a,j] = (-1)*m1b[j]; }    // Insert -1 times the column sum into ath row
    muAB[a,b] = sum(m1a);  // Insert the total into ath row, bth column, note: not minus
    return(muAB);
  }
}
data { 
  int<lower=0>  n;                             // Total number of data 
  int<lower=0>  a;                             // Number of levels for A
  int<lower=0>  b;                             // Number of levels for B
  int<lower=0>  s;                             // Number of blocks
  vector[n]     y;                             // Criterion variable   
  int<lower=0>  A[n];                          // Classification variable A
  int<lower=0>  B[n];                          // Classification variable B
  int<lower=0>  S[n];                          // Block variable S
}
parameters {
  real                           mu;           // Overall mean
  vector                   [a-1] m1A;          // One less mean for A
  vector                   [b-1] m1B;          // One less mean for B
  vector                   [s-1] m1S;          // One less mean for S
  matrix               [a-1,b-1] m1AB;         // One less interaction effect
  matrix               [a-1,s-1] m1AS;         // One less AS interaction effect
  matrix               [b-1,s-1] m1BS;         // One less BS interaction effect
  real<lower=0>             s_E;               // E standard deviation
  real<lower=0>             s_S;               // S standard deviation
  real<lower=0>             s_AS;              // AS standard deviation
  real<lower=0>             s_BS;              // BS standard deviation
}
transformed parameters {
  vector                   [a] muA;            // Mean for A
  vector                   [b] muB;            // Mean for B
  matrix                 [a,b] muAB;           // Interaction effect
  muA = zero_sum_vector(a, m1A);               // Equation (21.12), second line, first item
  muB = zero_sum_vector(b, m1B);               // Second item
  muAB = zero_sum_matrix(a, b, m1AB);          // Third and fourth items
}
model {
  vector                   [s] muS;            // Mean for S
  matrix                 [a,s] muAS;           // AS interaction effect
  matrix                 [b,s] muBS;           // BS interaction effect
  m1S ~ normal(0, s_S);                        // Equation (21.12), first line below, second item
  muS = zero_sum_vector(s, m1S);               // Equation (21.13), second line below, first item
  for(j in 1:(a-1)){ for(l in 1:(s-1)){m1AS[j,l] ~ normal(0, s_AS);}} // Third item
  for(k in 1:(b-1)){ for(l in 1:(s-1)){m1BS[k,l] ~ normal(0, s_BS);}} // Fourth item
  muAS = zero_sum_matrix(a, s, m1AS);          // Make the sum of as zero
  muBS = zero_sum_matrix(b, s, m1BS);          // Make the sum of bs zero
  for(i in 1:n){
     y[i] ~ normal(mu + muA[A[i]] + muB[B[i]] + muS[S[i]] +
            muAB[A[i],B[i]] + muAS[A[i],S[i]] + muBS[B[i],S[i]], s_E); // Equation (21.12)
}
generated quantities{
  real s_a;                                    // SD for factor A
  real s_b;                                    // SD for factor B
  real s_ab;                                   // SD for AB interaction effect
  s_a = sqrt(variance(muA)*(a-1)/a);           // SD for factor A
  s_b = sqrt(variance(muB)*(b-1)/b);           // SD for factor B
  s_ab = sqrt(variance(muAB)*((a*b)-1)/(a*b)); // SD for AB interaction effect
}';

par <- c("mu", "muA", "muB", "muAB", "s_S", "s_AS", "s_BS", "s_E", "s_a", "s_b", "s_ab")  # Parameters
dataSetRBFD <- list(n=length(PersuasiveMessages$Rating), a=max(PersuasiveMessages$Message),
                    b=max(PersuasiveMessages$Font), s=max(PersuasiveMessages$Subject), y=PersuasiveMessages$Rating,
                    A=PersuasiveMessages$Message, B=PersuasiveMessages$Font, S=PersuasiveMessages$Subject)  # Input

########### Execution with cmdstanr
modfileRBFD <- write_stan_file(RBFD)             # Writing temporary file
modRBFD <- cmdstan_model(modfileRBFD)            # Compilation
csrfitRBFD <- modRBFD$sample(data=dataSetRBFD, chains=5, iter_sampling=20000,
                             iter_warmup=1000, parallel_chains=5, seed=1234)  # MCMC
fitRBFD <- rstan::read_stan_csv(csrfitRBFD$output_files()) # Conversion to stan format

########### Execution with rstan
#fitRBFD <- stan(model_code=RBFD, data=dataSetRBFD,
#                pars=par, seed=1234, chains=5, warmup=1000, iter=21000)
#save(fitRBFD, file="./scrT/obje/stan_obje0803")
#load(file="./scrT/obje/stan_obje0803"); # Loading the object created beforehand with stan()

extRBFD <- extract(fitRBFD, par);  # Extracting random numbers

###### Table 21.6 Posterior Distribution of Parameters in Randomized Block Design
gqcal(extRBFD$mu       )
gqcal(extRBFD$muA[,1]   )
gqcal(extRBFD$muA[,2]   )
gqcal(extRBFD$muA[,3]   )
gqcal(extRBFD$muB[,1]   )
gqcal(extRBFD$muAB[,1,1])
gqcal(extRBFD$muAB[,2,1])
gqcal(extRBFD$muAB[,3,1])
gqcal(extRBFD$s_S      )
gqcal(extRBFD$s_AS     )
gqcal(extRBFD$s_BS     )
gqcal(extRBFD$s_E      )
print('------------------------------------------')
gqcal(extRBFD$s_a      )
gqcal(extRBFD$s_b      )
gqcal(extRBFD$s_ab     )

###### 21.3.4 PHC Probability of Difference Between Levels of 'Type'
phc02(0, extRBFD$muA)

# Probability of Difference Between Levels of 'Font'
mean(extRBFD$muB[,2] - extRBFD$muB[,1] > 0)

# Probability of Conjunction Proposition Being True 'Penalty' < 'Direct' and 'Direct' < 'Public Spirit'
round(mean((extRBFD$muA[,1] < extRBFD$muA[,2]) & (extRBFD$muA[,2] < extRBFD$muA[,3])), 3)

###### Figure 21.8 Curve of Probability of Difference Greater Than c in Level Effects
par(mfrow=c(3,1))
phc01(seq01=seq(0.5,1.4,0.05),a=extRBFD$muA[,2],b=extRBFD$muA[,1],
      xlab="c<a2-a1",cex.lab=2.5)
phc01(seq01=seq(1,2,0.05),a=extRBFD$muA[,3],b=extRBFD$muA[,2],
      xlab="c<a3-a2",cex.lab=2.5)
phc01(seq01=seq(1,2,0.05),a=extRBFD$muB[,2],b=extRBFD$muB[,1],
      xlab="c<b2-b1",cex.lab=2.5)
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_8.pdf")

###### Probability of a Substantial Difference
phc01(seq01=0.8, a=extRBFD$muA[,2], b=extRBFD$muA[,1], byoga="no")
phc01(seq01=1.3, a=extRBFD$muA[,3], b=extRBFD$muA[,2], byoga="no")
phc01(seq01=1.5, a=extRBFD$muB[,2], b=extRBFD$muB[,1], byoga="no")

################## 21.4 Estimation of Mixed Two-Factor Design
###### 21.4.1 Follow-Up Survey of Health Examinations

###### Table 21.7 Follow-Up Survey of Health Examinations
(DepressionFollowUp <- read.csv("./scrT/dat/depression_follow_up.csv", header = TRUE))

###### Figure 21.9 Scatter Plot of Follow-Up Survey of Depression Test Made Clearer by Adding Random Numbers
ma<-57;mi<-48
a1<-DepressionFollowUp[  1: 50,1]+runif(50)
a2<-DepressionFollowUp[ 51:100,1]+runif(50)
b1<-DepressionFollowUp[101:150,1]+runif(50)
b2<-DepressionFollowUp[151:200,1]+runif(50)
plot(a1,a2,pch=21,xlim=c(mi,ma),ylim=c(mi,ma),xlab="pre.test",ylab="post.test")
abline(0,1)
par(new=T)
plot(b1,b2,pch=19,xlim=c(mi,ma),ylim=c(mi,ma),xlab="",ylab="")
legend(49,56,c("Control","Intervention"),pch=c(21,19));  
#dev.copy2pdf(file="../English_tex/fig/chap21/fig21_9.pdf")

##################### Estimation of Mixed Two-Factor Design STAN Code
SPD <- '
functions{
  vector zero_sum_vector(int a, vector m1A){// Main effect zero-sum
    vector[a] muA;            
    for(i in 1:(a-1)){ muA[i] = m1A[i];}   // Copy up to a-1 as it is
    muA[a] = -sum(m1A);                    // For the ath, insert -1 times the sum
    return(muA);
  }
  matrix zero_sum_matrix(int a, int b, matrix m1AB){// Interaction effect zero-sum
    vector [a-1] m1a;                       // Vector for row sum
    vector [b-1] m1b;                       // Vector for column sum
    matrix [a,b] muAB;                      // Matrix to be returned
    for(i in 1:(a-1)){ m1a[i] = 0.0;}      // Zero clear
    for(j in 1:(b-1)){ m1b[j] = 0.0;}      // Zero clear
    for(i in 1:(a-1)){ for(j in 1:(b-1)){muAB[i,j] = m1AB[i,j];}} // Copy as it is
    for(i in 1:(a-1)){ for(j in 1:(b-1)){m1a[i] += m1AB[i,j];}} // Create row sum
    for(j in 1:(b-1)){ for(i in 1:(a-1)){m1b[j] += m1AB[i,j];}} // Create column sum
    for(i in 1:(a-1)){muAB[i,b] = (-1)*m1a[i]; }    // Insert -1 times the row sum into bth column
    for(j in 1:(b-1)){muAB[a,j] = (-1)*m1b[j]; }    // Insert -1 times the column sum into ath row
    muAB[a,b] = sum(m1a);  // Insert the total into ath row, bth column, note: not minus
    return(muAB);
  }
}
data { 
  int<lower=0>  n;                             // Total number of data 
  int<lower=0>  a;                             // Number of levels for A (between-subject factor)
  int<lower=0>  b;                             // Number of levels for B (within-subject factor)
  int<lower=0>  s;                             // Number of blocks
  vector[n]     y;                             // Criterion variable   
  int<lower=0>  A[n];                          // Classification variable A (between-subject factor)
  int<lower=0>  B[n];                          // Classification variable B (within-subject factor)
  int<lower=0>  S[n];                          // Block variable S
}
parameters {
  real                           mu;           // Overall mean
  vector                   [a-1] m1A;          // One less mean for A
  vector                   [b-1] m1B;          // One less mean for B
  vector                 [s*a-1] m1S;          // One less mean for S
  matrix               [a-1,b-1] m1AB;         // One less interaction effect
  real<lower=0>             s_E;               // E standard deviation
  real<lower=0>             s_S;               // S standard deviation
}
transformed parameters {
  vector                   [a] muA;            // Mean for A
  vector                   [b] muB;            // Mean for B
  matrix                 [a,b] muAB;           // Interaction effect
  muA = zero_sum_vector(a, m1A);               // Equation (21.13), second line below, first item
  muB = zero_sum_vector(b, m1B);               // Second item
  muAB = zero_sum_matrix(a, b, m1AB);          // Third and fourth items
}
model {
  vector                 [s*a] muS;            // Mean for S
  m1S ~ normal(0, s_S);                        // Equation  (21.13), first line below, second item
  muS = zero_sum_vector(s*a, m1S);             // Equation  (21.13), second line below, first item
  for(i in 1:n){
    y[i] ~ normal(mu + muA[A[i]] + muB[B[i]] + muS[S[i]] +
                    muAB[A[i],B[i]], s_E);}                // Equation (21.13)
}
generated quantities{
  real s_a;                                    // SD for factor A
  real s_b;                                    // SD for factor B
  real s_ab;                                   // SD for AB interaction effect
  s_a = sqrt(variance(muA)*(a-1)/a);           // SD for factor A
  s_b = sqrt(variance(muB)*(b-1)/b);           // SD for factor B
  s_ab = sqrt(variance(muAB)*((a*b)-1)/(a*b)); // SD for AB interaction effect
}
';

par <- c("mu", "muA", "muB", "muAB", "s_S", "s_E", "s_a", "s_b", "s_ab")  # Parameters
dataSetSPD <- list(n=length(DepressionFollowUp$Score), a=max(DepressionFollowUp$GroupIndex),
                   b=max(DepressionFollowUp$Time), s=max(DepressionFollowUp$Student), y=DepressionFollowUp$Score,
                   A=DepressionFollowUp$GroupIndex, B=DepressionFollowUp$Time, S=DepressionFollowUp$Student)  # Input

########### Execution with cmdstanr
modfileSPD <- write_stan_file(SPD)            
modSPD <- cmdstan_model(modfileSPD)                
csrfitSPD <- modSPD$sample(data = dataSetSPD,chains = 5,iter_sampling = 20000,
                           iter_warmup = 1000,parallel_chains = 5,seed=1234)  
fitSPD <- rstan::read_stan_csv(csrfitSPD$output_files()) 

extSPD <- extract(fitSPD, par);  # Extracting random numbers

###### Table 21.8 Summary of Posterior Distribution of Parameters in Mixed Design
gqcal(extSPD$mu        )
gqcal(extSPD$muA[,1]   )
gqcal(extSPD$muB[,1]   )
gqcal(extSPD$muAB[,1,1])
gqcal(extSPD$muAB[,1,2])
gqcal(extSPD$muAB[,2,1])
gqcal(extSPD$muAB[,2,2])
gqcal(extSPD$s_S       )
gqcal(extSPD$s_E       )
print('------------------------------------------')
gqcal(extSPD$s_a       )
gqcal(extSPD$s_b       )
gqcal(extSPD$s_ab      )

#### Posterior Distribution of Generated Quantities (Effect of Cell jk in Group j at Time k) (Equation 21.15)
ce11<-extSPD$muA[,1]+extSPD$muB[,1]+extSPD$muAB[,1,1]
ce12<-extSPD$muA[,1]+extSPD$muB[,2]+extSPD$muAB[,1,2]
ce21<-extSPD$muA[,2]+extSPD$muB[,1]+extSPD$muAB[,2,1]
ce22<-extSPD$muA[,2]+extSPD$muB[,2]+extSPD$muAB[,2,2]


###### Table 21.9 Posterior Distribution of Within-Group Differences Over Time
gqcal(ce11-ce12);             # equation (21.16)
gqcal(ce21-ce22);

###### Histogram of Posterior Distribution of Within-Group Differences Over Time (Not in the textbook)
hist(ce11-ce12,main="",ylab="")
hist(ce21-ce22,main="",ylab="")

################################# Practical assignment
#1
#def32
phc01(seq01=seq(-3,6,0.05),a=def32,cc="gtc",byoga="yes",xlab="d32",cex.lab=2.5)
phc01(seq01=seq(-3,6,0.5 ),a=def32,cc="gtc",byoga="no")

#def32,ROPE
phc01(seq01=seq(0,6,0.05),a=def32,cc="rope",byoga="yes",xlab="|d32|<c",
      cex.lab=2.5)
phc01(seq01=seq(0,6,0.5 ),a=def32,cc="rope",byoga="no")

#def13
phc01(seq01=seq(17,25,0.05),a=def13,cc="gtc",byoga="yes",xlab="d32",
      cex.lab=2.5)
phc01(seq01=seq(17,25,0.5),a=def13,cc="gtc",byoga="no")


#2

#
#     Correct answers are omitted
#

