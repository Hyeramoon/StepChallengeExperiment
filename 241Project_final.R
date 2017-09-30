setwd('~/Documents/Berkeley/Courses/w241/Project')

# 4 SECTIONS:
# SECTION 1: LINEAR REGRESSION MODELS - WITHOUT CLUSTERING AND COMPLIANCE
# SECTION 2: LINEAR REGRESSION MODELS - WITH COMPLIANCE (NO CLUSTERING)
# SECTION 3: LINEAR REGRESSION MODELS - WITH CLUSTERING AND COMPLIANCE
# SECTION 4: COVARIATES
# SECTION 5: ATTRITION

# NOTES:
# 1. pre experiment average was taken over 7 days - which means only one weekend. Thus to compare apples to apples, the average steps 
#    on the RHS of lm is also calculated on last 7 days of the Step Challenge ("step_avg3"). Note step_avg is taking into account all days
#    of the challenge and step_avg2 is the average from July 28 (first treatment email sent) to end 
# 2. Our model should be analyzing the difference in differences, i.e., increase in control vs increase in treatment groups: thus the
#    coefficients of the treatment variables is simply a shift in the regression line (ie how much more from control group in average)
#    but the coeffs of the interaction terms in the slope of the line (i.e for each step of pre-challenge, how many more step they took
#    from receiving the treatment)



###### Load Packages and Data ############
library(data.table)
library(foreign)
library(lmtest)
library(sandwich)
library(multiwayvcov)
library(stargazer)

d <- read.csv("241dataPrePost.csv")
head(d)
tail(d)


################## SECTION 1: LINEAR REGRESSION MODELS - WITHOUT CLUSTERING AND COMPLIANCE ###########################################

# Our model should be analyzing the difference in differences, i.e., increase in control vs increase in treatment groups

# 1. Difference in difference equation
#    Intercept: average step for control group in pre-treatment period
#    Meaning of coefficient of treat_positive and treat_negative: average steps people in the treatment group took more compared to control in pre-period
#    Meaning of coefficient of time_post: average steps the people in the control took more from pre-period to post-period
#    Interaction term: change in step during experiment from treatment (that's what we want to measure)
mod1 <- lm(step_avg_prepost ~ 1 + treat_positive + treat_negative + time_post + treat_positive*time_post + treat_negative*time_post, data=d)
summary(mod1)


# Note: num_read_emails will not be used since it may be conditional on treatment


############################## SECTION 2: COMPLIANCE - CACE ######################################################################

# Since we know who read (i.e. treatment email delivered successfully to subjects) and who didn't, not only for the treatment groups,
# but also for the control group, we have a placebo design. Thus, instead of using 2-stage linear models with "predict", we will drop
# all the non compliers (i.e. those you never read emails)

d_comp = subset(d, Compliance==1)  # Drop the non-compliers from dataset

mod3 <- lm(step_avg_prepost ~ 1 + treat_positive + treat_negative + time_post + treat_positive*time_post + treat_negative*time_post, data=d_comp)
summary(mod3)


# Print models side by side for easier comparison
stargazer(mod1, mod3, omit = "street", type = "text")


############################## SECTION 3: CLUSTERING ######################################################################

# Clustered standard errors
mod3$cluster.vcov <- cluster.vcov(mod3, ~ d_comp$Cluster_ID)
mod3$cluster.vcov
coeftest(mod3, mod3$cluster.vcov)


## Extract the SE from the VCOV (we need the square root of the diagonals of the VCOV matrix)
mod3$cluster.se <- sqrt(diag(mod3$cluster.vcov))
mod3$cluster.se


## Compute the OLS SE (i.e. not clustered SE)
# for model 3 with compliers only
mod3$ols.vcov <- vcovHC(mod3, "const")
mod3$ols.se   <- sqrt(diag(mod3$ols.vcov))
mod3$ols.se

# for model 1 with compliers and non-compliers
mod1$ols.vcov <- vcovHC(mod1, "const")
mod1$ols.se   <- sqrt(diag(mod1$ols.vcov))


stargazer(mod1,mod3, mod3,
          se = list(mod1$ols.se,
                    mod3$ols.se, 
                    mod3$cluster.se),
          type = "text")

# Notes: SEs from clustering are smaller than than SEs from no clustering


############################## SECTION 4: COVARIATES ######################################################################

#  Purpose: add covariates which are the most predictive of the outcomes (step average during challenge) in order to narrow the
#           the SEs of the coefficients of the treatment effects estimates (i.e. more precise estimates)

### 1. Covariates balance check
d_post = subset(d, time_post==1)  # for 84 observations. Note, the dataset d has double observations since pre and post periods covered

# 1. positive treatment group
check1 <- lm(treat_positive ~ Gender + used_ctr_before + activity_level_num 
             + age + weight, data=d_post)
summary(check1)

check1b <- lm(treat_positive ~ Gender + used_ctr_before + activity_level_num 
              + age + weight, data=d)
summary(check1b)
# Results: 
# no coeffs are significant (except weight) and also F-test p-value > 0.05 (thus mean of both population - positive and not positive are same)


# 2. negative treatment group
check2 <- lm(treat_negative ~ Gender + used_ctr_before + activity_level_num 
             + age + weight, data=d_post)
summary(check2)
check2b <- lm(treat_negative ~ Gender + used_ctr_before + activity_level_num 
              + age + weight, data=d)
summary(check2b)
# Results: 
# no coeffs are significant (except age)
# F-test p-value > 0.05 (thus mean of both population - negative and not negative are same) with 84 observations
# BUT F-test p-value < 0.05 (thus mean of both population - negative and not negative are not the same) with 168 observations - smaller SE

stargazer(check1, check1b, check2, check2b, omit = "street", type = "text")
stargazer(check1, check2, omit = "street", type = "text")


### 2. Adding covariates to the model
mod4 <- lm(step_avg_prepost ~ 1 + treat_positive + treat_negative + time_post + treat_positive*time_post + treat_negative*time_post 
           + Gender_num + activity_level_num
           + age + weight, data=d_comp)  # note: only compliers data

# Clustered standard errors
mod4$cluster.vcov <- cluster.vcov(mod4, ~ d_comp$Cluster_ID)
mod4$cluster.se <- sqrt(diag(mod4$cluster.vcov))


# Print models side by side for easier comparison
stargazer(mod1,mod3, mod3, mod4,
          se = list(mod1$ols.se,
                    mod3$ols.se, 
                    mod3$cluster.se,
                    mod4$cluster.se),
          type = "text")

# Observations: The covariates are not really helping here, in particular for coeffs interpretation become more difficult
#               (the covariate weight is a continuous measure and negative for each pound)
#              ** Thus, covariates not added in our final model**


###### ATTRITION CHECK #####
da <- read.csv("241dataZeroStep.csv")
head(da)

# num_step0 is the number of days with step count = 0
mod5 <- lm(num_step0 ~ treat_positive, data=da)
summary(mod5)

mod6 <- lm(num_step0 ~ treat_negative, data=da)
summary(mod6)

da$num_step0_bin <- as.numeric(da$num_step0!=0)  # dummy variable for if there days with steps count = 0
mod5a <- lm(num_step0_bin ~ treat_positive, data=da)
summary(mod5a)

mod6a <- lm(num_step0_bin ~ treat_negative, data=da)
summary(mod6a)

# Observations: except for model 5a (positive treatment group) with binary dummy variable num_step0_bin, no significant results
#               for number of days of step=0 to test for possible differential attrition

