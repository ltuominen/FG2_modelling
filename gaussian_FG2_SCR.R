# libraries
lapply(c('brms', 'tidyverse', 'tidybayes', 'reshape2', 'lme4'), library, character.only = TRUE, quietly = TRUE) 

# load data 
long <- read_csv( '/Users/laurituominen/Documents/Research/FG2/DATA/SCRanon.csv')

# plot data 
ggplot(data    = long,
       aes(x   = morphs,
           y   = SCRz,
           col = subject_number))+
  geom_point(size     = 1.2,
             alpha    = .6,
             position = "jitter",
             width = 1)+ 
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  theme_minimal()+
  ylab('Skin Conductance Response (z-score)') +
  xlab('Distance from Conditioned Stimulus') +
  theme(legend.position = "none")+
  labs(title    = "SCR vs. Events")

# test some reasonable priors 
e =seq(0,100, by=0.1)
h=0.614
gauss <- function(w) { 
  e =seq(0,100, by=0.1)
  h=0.614
  y= h * exp( -1 * ( e^2 / (2*w^2 )^2 ))
  return(y)
  } 
w1 <- gauss(w=1); w3 <- gauss(w=3); w45<-gauss(w=4.5); w6 <- gauss(w=6); w9 <- gauss(w=9)
df <- tibble(e=e, w1=w1, w3=w3, w45=w45, w6=w6, w9=w9)

# plot priors 
ggplot(data=df)+
  geom_line(aes(x=e, y=w1, col='red'))+
  geom_line(aes(x=e, y=w3, col='green'))+
  geom_line(aes(x=e, y=w45, col='cyan'))+
  geom_line(aes(x=e, y=w6, col='blue'))+
  geom_line(aes(x=e, y=w9, col='magenta'))+
  scale_color_discrete( name = "W",
                       breaks = c("red", "green", "cyan", "blue", "magenta"),
                       labels = c("1", "3", "4.5", "6", "9")) +
  ylab("response") +
  xlab("distance")

  
# write in the form 
nlform2 <- bf(SCRz ~ i + h * exp( -1 *  (morphs ^2 / (2*w^2 )^2 )), 
              i ~ (1 | subject_number),
              h ~ 1, 
              w ~ 1,
              nl = TRUE)

# choose maybe appropriate priors 
nlprior2 <- c(prior(normal(0.5,0.2), nlpar = "h"),
              prior(cauchy(4,1), nlpar = "w"))

# fit 
gaussian_fit <- brm(formula = nlform2, data = long, family = gaussian(),
                    prior = nlprior2, control = list(adapt_delta = 0.995, max_treedepth = 15),
                    iter = 20000, warmup = 2000, chains = 4, cores = 4)

prior_summary(gaussian_fit)
summary(gaussian_fit)

# plot fit and data
means <- long %>%
  select(-c('morphs')) %>%
  pivot_wider(values_from = SCRz, names_from = Event, names_prefix = 'M')  %>%   
  summarise(across(-c("subject_number"), ~ mean(.x)))

morphs <- unique(long$morphs)
df <- tibble(meanSCR = t(means), morphs=morphs)

p <- conditional_effects(gaussian_fit,
                    spaghetti = T,
                    nsamples = 200) 
  
ugly_plot <- plot(p, plot = FALSE)[[1]] 
ugly_plot + 
  geom_point(data=df, aes(y=meanSCR, x=morphs), alpha=3/4, inherit.aes=FALSE)+
  labs(title = "Gaussian fit to the SCR") + 
  ylab('Skin Conductance Response (z-score)') +
  xlab('Distance from Conditioned Stimulus') +
  theme(text=element_text(size =16)) 

