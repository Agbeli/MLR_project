######  Project source code
#### libraries required 
library(ggplot2)
library(dplyr)
library(rlang)
library(lubridate)
library(gt)
library(corrplot)
library(tidyverse)
library(GGally)
library(tidyr)
library(DescTools)
library(leaps)
library(performance)
library(car)
library(MASS)

### set seed 
set.seed(2024)
### path to the datasets
filename = "DataSets/Huntington_MAS_pkg/Huntington_2019_calibration_data.csv"

### load the datsets
data_huntington <- read.csv(filename)

###  check the first five rows of the datasets
head(data_huntington)

### check the names of the columns: piping the process:
data_huntington %>% 
  names()

##################################################### Exploratory Data Analysis Pipeline ##############################################################
### check for missing values in dataframe  
data_huntington %>% summarise(across(everything(), ~ sum(is.na(.))))


### find the total and average Ecoli recorded over the years from 2005 to 2015
ecoli_summary <- function(data){
  result <- data %>%
    mutate(date_column = mdy(Date)) %>%  
    mutate(year = year(date_column)) %>%                                  
    group_by(year)%>%
    summarise(total_yearly_bateria = sum(EcoliAve_CFU), 
              average_yearly_bacteria = mean(EcoliAve_CFU))
  return(result)
}

### print out the results 
print(ecoli_summary(data = data_huntington))

#### summarize and plot the total and average of Ecoli concenntration 
df_ecoli_summary <- data_huntington %>%
  mutate(date_column = mdy(Date)) %>%  
  mutate(year = year(date_column)) %>%                                  
  group_by(year)%>%
  summarise(total_yearly_bateria = sum(EcoliAve_CFU), 
            average_yearly_bacteria = mean(EcoliAve_CFU)) %>%
  pivot_longer(cols = c(total_yearly_bateria, average_yearly_bacteria),  # Reshape for plotting
               names_to = "Type", 
               values_to = "Value")

### plot the results of the summarize table for the concentration
ggplot(df_ecoli_summary, aes(x = factor(year), y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Value, 1)),     # Combine Type and Value
            position = position_dodge(width = 0.8), 
            vjust = -0.3, 
            size = 3.5) + 
  labs(title = "Total and Average of Escherichia coli Concentration recorded by Year",
       x = "Year",
       y = "Value") +
  scale_fill_manual(values = c("skyblue", "darkblue"), labels = c("Average", "Total")) +
  theme_minimal()

### boxplot for raw data E.coli concentration by year 
data_huntington %>%
  mutate(date_column = mdy(Date)) %>%
  mutate(year = year(date_column)) %>%
  mutate(year = format(year))%>%
  ggplot(aes(x = year, y = EcoliAve_CFU)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Boxplot of EcoliAve_CFU by Year", x = "Year", y = "Value") +
  theme_minimal()


### Create a boxplot log transformation grouped by year
data_huntington %>%
  mutate(date_column = mdy(Date)) %>%
  mutate(year = year(date_column)) %>%
  mutate(year = format(year))%>%
  ggplot(aes(x = year, y = log1p(EcoliAve_CFU))) +
    geom_boxplot(fill = "lightblue", color = "black") +
    labs(title = "Boxplot of log(EcoliAve_CFU+1) by Year", x = "Year", y = "Value") +
    theme_minimal()

### histogram of the raw dataset: 
ggplot(data_huntington, aes(x=EcoliAve_CFU))+
  geom_histogram(color="darkblue", fill="lightblue",bins = 20,na.rm = TRUE) +
  ggtitle(label = "Histogram plot of Raw EcoliAve_CFU")

#### visualize the histogram distribution of the transformed target : Preprocessed data and perform winsorization on the few outlier
ggplot(data_huntington, aes(x=log1p(EcoliAve_CFU)))+
  geom_histogram(color="darkblue", fill="lightblue",bins = 30,na.rm = TRUE) +
  ggtitle(label = "Histogram plot of log(EcoliAve_CFU+1)")


#### summary statistic of each column
df_summary <- data_huntington %>%
  select(where(is.numeric)) %>%  # Select only numeric and integer columns
  summarise(
    across(
      everything(),
      list(
        mean = ~ mean(.),
        sd = ~ sd(.),
        min = ~ min(.),
        max = ~ max(.),
        median = ~ median(.)
      ),
      .names = "{col}_{fn}"
    )
  )

print(df_summary) ### output the summary statistics 

# Display in tabular format
#kable(df_summary, caption = "Summary Statistics") alternative apporach to summarise the datasets
df_summary %>%
  gt() %>%
  tab_header(
   title = "Summary Statistics for each features"
   )


#### correlation coefficient of the features  
corr_matrix <- data_huntington %>% select_if(is.numeric) %>% 
  cor(method="pearson", use="pairwise.complete.obs") %>%
  corrplot(method="ellipse",  
           type="lower", order="hclust", 
           insig = "p-value",
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation, 
           #title = "Correlation between EcoliAve_CFU and Predictors"
  ) 


### Select numeric columns
df_numeric <- data_huntington[sapply(data_huntington, is.numeric)]

### select all the predictors in the dataframe 
df_predictors <- df_numeric %>%
  select(c("Lake_Temp_C","Lake_Turb_NTRU","WaveHt_Ft","LL_PreDay","AirportRain48W_in"))

### Create pairplot using GGally's ggpairs
ggpairs(df_numeric,aes(color = "darkblue", alpha = 0.5),upper = list(continuous = "points"),
        title = "Scatterplot among target and predictors",
        diag = list(continuous = "densityDiag")) +
  theme_minimal()

### Reshape the data into long format
df_long <- df_predictors %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

### Create boxplot for each features 
ggplot(df_long, aes(x = Variable, y = Value)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  theme_minimal() +
  labs(title = "Boxplot for Each predictor", x = "Predictors", y = "Values")


################################ Data Preprocessing pipeline : Feature engineering the features. ######################################################
#######################################################################################################################################################
data_huntington_cp <- data.frame(data_huntington) ### make a copy of the original datasets 

### Features transformation pipeline 
data_huntington_cp<- data_huntington_cp %>% 
  mutate_at(c("EcoliAve_CFU","Lake_Turb_NTRU"), log1p)


### Define fourth function 
fourth_root <- function(x){
  return(x^(1/4))
}

### apply the forth root to the features 
data_huntington_cp<- data_huntington_cp %>% 
  mutate_at(c("WaveHt_Ft","AirportRain48W_in"), fourth_root)


### convert the date into numeric  and after that extract the day of recorded events #######
data_huntington_cp <- data_huntington_cp %>%
  mutate(date_column = mdy(Date)) %>%
  mutate(date_column_numeric = as.numeric(day(date_column)))


df_numeric_cp <- data_huntington_cp %>%
  select(c("EcoliAve_CFU","Lake_Temp_C","Lake_Turb_NTRU","WaveHt_Ft","LL_PreDay","AirportRain48W_in"))

### Create pairplot using GGally's ggpairs: Visualization pipeline
ggpairs(df_numeric_cp,aes(color = "darkblue", alpha = 0.5),upper = list(continuous = "points"),
        title = "Features Transformation Scatterplot among target and predictors",
        diag = list(continuous = "densityDiag")) +
  theme_minimal()

###  Time series plot of E.coli concentration 
p <- ggplot(data_huntington, aes(x=mdy(Date), y=EcoliAve_CFU)) +
  geom_line(color="turquoise4") +
  theme_minimal() + 
  labs(x="Daily records", y="E.coli", title="E.coli concentration from 2005 to 2015") +
  theme(plot.title = element_text(hjust=0.5, size=20, face="bold"))

p ### print out the graph 

### Boxplot after feature engineering 
### Reshape the data into long format
df_numeric_cpf <- data_huntington_cp %>%
  select(c("Lake_Temp_C","Lake_Turb_NTRU","WaveHt_Ft","LL_PreDay","AirportRain48W_in"))

#### Restructure the predictors to visualize the boxplot 
df_long_cpf <- df_numeric_cpf %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

### Create boxplot for each features 
ggplot(df_long_cpf, aes(x = Variable, y = Value)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  theme_minimal() +
  labs(title = "Boxplot for Each transformed predictor", x = "Predictors", y = "Values")


#### correlation transformation: Visualization pipeline 
corr_matrix_cp <- data_huntington_cp %>% select_if(is.numeric) %>% 
  cor(method="pearson", use="pairwise.complete.obs") %>%
  corrplot(method="ellipse",  
           type="lower", order="hclust", 
           insig = "p-value",
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation, 
           #title = "Correlation between EcoliAve_CFU and Predictors"
  ) 


##### Cap the outliers with the max values using winsorization.: Outlier manipulation.. 
winsorize2 <- function(x) {
  Min <- which(x == min(x))
  Max <- which(x == max(x))
  ord <- order(x)
  x[Min] <- x[ord][length(Min)+1]
  x[Max] <- x[ord][length(x)-length(Max)]
  x
}

###### applied the winsorization function to the features of interest
data_huntington_cp <- data_huntington_cp %>% 
  mutate_at(c('EcoliAve_CFU', 'Lake_Temp_C', 'Lake_Turb_NTRU', 'WaveHt_Ft','LL_PreDay',
              'AirportRain48W_in','date_column_numeric'), winsorize2)



################################# Model Fitting Pipeline after feature engineering #################################################################
###################################################################################################################################################
model_processed <- lm(EcoliAve_CFU ~ Lake_Temp_C + Lake_Turb_NTRU + WaveHt_Ft + LL_PreDay  + AirportRain48W_in, 
                     data = data_huntington_cp)

summary(model_processed)

########## 
model_performance(model_processed)


#### fit the model with time factor  
model_processed_with_date <- lm(EcoliAve_CFU ~ Lake_Temp_C + Lake_Turb_NTRU + WaveHt_Ft + LL_PreDay  + AirportRain48W_in + date_column_numeric, 
                      data = data_huntington_cp)

summary(model_processed_with_date)

######### General Linear F-Test 
anova(model_processed,model_processed_with_date)

################################################## Model Assessment pipeline ###############################################################################
huntington_cp_models <- regsubsets(EcoliAve_CFU ~ Lake_Temp_C + Lake_Turb_NTRU + WaveHt_Ft + LL_PreDay  + AirportRain48W_in + date_column_numeric, nbest=6, data=data_huntington_cp)

#### Visual comparison by adjusted R-squared, 
plot(huntington_cp_models, scale="adjr2",col = "orange")

### BIC 
plot(huntington_cp_models, scale="bic",col = "blue")

############################################### Model Diagnostic Pipeline ######################################################################
vif(model_processed)

# Extract fitted values and studentized residuals
fitted_values <- fitted(model)
studentized_residuals <- rstudent(model)

#### extract the fitted values 
fitted_values <- fitted(model_processed)

#### studentized residuals from the model without time factors 
studentized_residuals <- rstudent(model_processed)

# Create a data frame for ggplot
plot_data <- data.frame(Fitted = fitted_values, Residuals = studentized_residuals)

# Plot using ggplot2
ggplot(plot_data, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Studentized Residuals vs. Fitted Values",
       x = "Fitted Values",
       y = "Studentized Residuals") +
  theme_minimal()


#### Create a data frame for ggplot
qq_data <- data.frame(
  Theoretical = qqnorm(studentized_residuals, plot.it = FALSE)$x,
  Sample = qqnorm(studentized_residuals, plot.it = FALSE)$y
)

#### Plot the QQ plot with ggplot2
ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  labs(title = "QQ Plot of Studentized Residuals",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

