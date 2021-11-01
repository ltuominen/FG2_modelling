library('tidyverse') 

data <- read_csv('/Users/laurituominen/Documents/Research/FG2/DATA/SCRgen.csv')
data <- rename(data, M9=`CSminus`, ID=`Row`) 
data$subject_number <- seq(1,dim(data)[1])


# get x-axis
morphs <- read_csv('/Users/laurituominen/Documents/Research/FG2/morphs.csv', col_names=FALSE)
reor <- ifelse( (morphs$X1 > morphs$X8), TRUE, FALSE  )
morphs<- morphs[reor,order(colnames(morphs),decreasing=TRUE)] 
xaxis <- c(as.vector(apply(morphs, 2, mean)), 100)


#data$ok <- ifelse( (data$M1-data$M9 ) > 0, TRUE, FALSE )
#data <- data[data$ok, ]
long <- data %>%
  pivot_longer(cols = starts_with("M"), names_to = 'Event', names_prefix = 'M', values_to = 'SCR',
               names_transform = list(Event = as.integer)) 

#long$logSCR <- log(long$SCR)
long$SCRz <-( long$SCR - mean(long$SCR) ) / sd(long$SCR) 

a2b <- function(a, b){b[a]} 
long$morphs <- unlist(lapply(long$Event,a2b, b=xaxis)  )

write_csv(x=long, file = '/Users/laurituominen/Documents/Research/FG2/DATA/SCRlong.csv')
long_anon <- long %>%
  select(-c('ID', 'SCR'))

write_csv(x=long_anon, file= '/Users/laurituominen/Documents/Research/FG2/DATA/SCRanon.csv')

