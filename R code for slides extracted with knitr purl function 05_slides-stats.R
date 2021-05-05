## ----setup, include=FALSE--------------------------------------------------------
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(
  fig.width=9, fig.height=3.5, fig.retina=3,
  out.width = "100%",
  cache = FALSE,
  echo = TRUE,
  message = FALSE, 
  warning = FALSE, 
  hiline = TRUE,
   comment = NA 
)
options(width = 70)

#Load data...
  #load(file = "Z:/Projects/Fridley_Brooke/2813_R_Programming_Course_2020/data/alldata.RData")
   load(file = "Z:/Projects/Fridley_Brooke/2813_R_Programming_Course_2020/data/Rclassdata.RData")
  
  library(easystats)
  library(tidyverse)
  library(kableExtra)
  library(ggplot2)
  library(janitor)
  library(gridExtra)
  library(broom)
   
   Rclassdata %>% filter(acronym %in% c("KIRC","HNSC")) -> tcga

   
  


## ----xaringan-extras, echo=FALSE, results="asis"---------------------------------
# remotes::install_github("gadenbuie/xaringanExtra")
xaringanExtra::use_xaringan_extra(c(
  "tile_view"
  # "editable",
  # "animate",
  # "panelset"
))


## ----xaringan-panelset, echo=FALSE-----------------------------------------------
xaringanExtra::use_panelset()


## /* Define title slide image or logo here */

## .talk-logo {

##   width: 400px;

##   height: 750px;

##   position: absolute;

##   top: 6%;

##   right: 7%;

##   /* background-image: url('img/r4ds-cover.png'); */

##   background-size: contain;

##   background-repeat: no-repeat;

##   background-position: contain;

## }


## ----comment=NA------------------------------------------------------------------

janitor::tabyl(tcga, radiation_therapy, vital_status )
   


## ----comment=NA------------------------------------------------------------------

janitor::tabyl(tcga, radiation_therapy, vital_status, 
               show_na = FALSE )



## ----comment=NA------------------------------------------------------------------

tcga %>% 
  tabyl(smoking, gender, show_na = FALSE )


## ----comment=NA------------------------------------------------------------------

tcga %>% 
  tabyl(smoking, gender, show_na = FALSE ) %>% 
  adorn_totals(where = c("row","col"))


## ----comment=NA------------------------------------------------------------------


tcga %>% 
  tabyl(smoking, gender, show_na = FALSE ) %>% 
  adorn_totals(where = c("row","col")) %>% 
  adorn_percentages(denominator = "col") 


## ----comment=NA------------------------------------------------------------------

tcga %>% 
  tabyl(smoking, gender, show_na = FALSE ) %>% 
  adorn_totals(where = c("row","col")) %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting(digits = 0)  



## ----comment=NA------------------------------------------------------------------

tcga %>% 
  tabyl(smoking, gender, show_na = FALSE ) %>% 
  adorn_totals(where = c("row","col")) %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting(digits = 0) %>% 
  adorn_ns(position = "front") 
 


## ----chisqtest, echo=FALSE ,  include=TRUE---------------------------------------
 tcga %>% 
  tabyl(smoking, gender, show_na = FALSE ) %>% 
  adorn_totals(where = c("row","col")) %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting(digits = 0)   %>% 
  adorn_ns(position = "front") 
 

 


## ----chisqtest1, eval=TRUE,  include=TRUE----------------------------------------
 
tcga %>% tabyl( smoking, gender , show_na = FALSE ) %>%
  chisq.test()



## ----geombar, echo = TRUE, fig.show = "hide"-------------------------------------
ggplot(tcga %>% filter(!is.na(smoking)),
       aes(x = gender, fill = smoking )) +  
  geom_bar(position = "fill") +
  labs(y = "proportion")



## ----ref.label = "geombar", echo = FALSE-----------------------------------------


## ----comment=NA------------------------------------------------------------------
 tcga %>% 
  tabyl( gender, vital_status , show_na = FALSE ) %>% 
  adorn_totals(where = c("row","col")) %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting(digits = 0) %>% 
  adorn_ns(position = "front") 



## ----comment=NA------------------------------------------------------------------
tcga  %>% tabyl( gender, vital_status , show_na = FALSE ) %>% 
  fisher.test()



## ----geombar2, echo = TRUE-------------------------------------------------------
ggplot(tcga %>% filter(!is.na(vital_status)), 
       aes(x = gender, fill = vital_status )) +  
       geom_bar(position = "fill") + labs(y = "proportion")



## ----comment=NA------------------------------------------------------------------

 
r1 <- tcga %>% correlation::cor_test("DUOXA1_exp", 
                                          "DUOX1_exp", 
                                          method = c("pearson") ) 
 
r2 <- tcga %>% correlation::cor_test("BUB1_exp",
                                          "C10orf32_exp",
                                          method = c("pearson") ) 
 
r3 <- tcga %>% correlation::cor_test("BRAF_exp", 
                                          "DTL_exp", 
                                          method = c("pearson") )
 


## ----comment=NA------------------------------------------------------------------
 knitr::kable(bind_rows(r1,r2,r3 ), format = 'html', digits = 3) %>%
  kable_styling(font_size = 12)



## ----geom, echo = TRUE, fig.show = "hide"----------------------------------------
ggplot(tcga) + 
  aes(DUOXA1_exp, DUOX1_exp) + 
  geom_point() + 
  annotate(geom = "text", x = 10, y = 25000, 
           label = paste("r = ", round(r1$r, 2), sep = ""),
           color = "red")



## ----ref.label = "geom", echo = FALSE--------------------------------------------


## ----geom2, echo = TRUE, fig.show = "hide"---------------------------------------

ggplot(tcga) + 
  aes(BUB1_exp, C10orf32_exp) + 
  geom_point() + 
  annotate(geom = "text", x = 1000, y = 4000, 
           label = paste("r = ", round(r2$r, 2), sep = ""),
           color = "red")
 



## ----ref.label = "geom2", echo = FALSE-------------------------------------------


## ----geom3, echo = TRUE, fig.show = "hide"---------------------------------------
ggplot(tcga) + 
  aes(BRAF_exp, DTL_exp) + 
  geom_point() + 
  annotate(geom = "text", x = 0, y = 1100, 
           label = paste("r = ", round(r3$r, 2), sep = ""),  
           color = "red")
 

 



## ----ref.label = "geom3", echo = FALSE-------------------------------------------


## ----comment=NA------------------------------------------------------------------




## ----comment=NA------------------------------------------------------------------
tcga %>% 
  select(DUOXA1_exp, DUOX1_exp, BUB1_exp, BRAF_exp, DTL_exp ) %>%
  cor() %>% round(.,2)



## ----echo=FALSE, comment=NA------------------------------------------------------
 

x <- seq(10,20,.25)
y <- (x-15)^3 + 100
pdata <-  data.frame(x,y)
 
ggplot(pdata, aes(x,y)) + geom_point()  



## ----comment=NA------------------------------------------------------------------

 
r1 <- tcga %>% correlation::cor_test("DUOXA1_exp", 
                                          "DUOX1_exp", 
                                          method = c("spearman") ) 
 
r2 <- tcga %>% correlation::cor_test("BUB1_exp",
                                          "C10orf32_exp",
                                          method = c("spearman") ) 
 
r3 <- tcga %>% correlation::cor_test("BRAF_exp", 
                                          "DTL_exp", 
                                          method = c("spearman") )
 


## ----comment=NA------------------------------------------------------------------
 knitr::kable(bind_rows(r1, r2, r3), format = 'html', digits = 3) %>%
  kable_styling(font_size = 12)



## ----ref.label = "wilcoxonhists", echo = FALSE-----------------------------------


## ----wilcoxonhists, echo = TRUE, fig.show = "hide"-------------------------------
mplot <- ggplot(tcga %>% filter(gender == "MALE"), 
                aes(x = BUB1_exp )) +
  geom_histogram(  colour = "black", position = "dodge") +
  ggtitle("Males") 

wplot <- ggplot(tcga  %>% filter(gender == "FEMALE"), 
                aes(x = BUB1_exp)) +
  geom_histogram(  colour= "black", position = "dodge")  +
  ggtitle("Females") 
  
grid.arrange(mplot,wplot, ncol=2) 
 



## ----wilcoxonboxplots, echo = TRUE, fig.show = "hide"----------------------------
  
ggplot(tcga , aes(gender, BUB1_exp, fill = gender))  +
  geom_boxplot() +
  scale_fill_manual( values = c("yellow", "blue")) + 
  theme(legend.position = "none") +
  labs(title = "BUB1 expression" , x = "Gender", y = "expression")  



## ----ref.label = "wilcoxonboxplots", echo = FALSE--------------------------------


## ----comment=NA------------------------------------------------------------------
 
 wilcox.test(BUB1_exp ~ gender, data = tcga,
                   alternative = c("two.sided"))
 
 


## ----t2boxplots, echo = TRUE, fig.show = "hide"----------------------------------
  
ggplot(tcga , aes(gender, PTEN_exp, fill = gender))  +
  geom_boxplot() +
  scale_fill_manual( values = c("yellow", "blue")) + 
  theme(legend.position = "none") +
  labs(title = "PTEN expression" , x = "Gender", y = "expression")  



## ----ref.label = "t2boxplots", echo = FALSE--------------------------------------


## ----comment=NA------------------------------------------------------------------
print(dim(tcga))
tcga <- tcga %>% filter(PTEN_exp < 4000)  
print(dim(tcga))


## ----t2hist, echo = TRUE, fig.show = "hide"--------------------------------------
  ggplot(tcga, aes(x=PTEN_exp, fill=gender)) +
  scale_fill_manual( values = c("red", "blue")) + 
  geom_histogram( alpha=0.6, position="identity")
 


## ----ref.label = "t2hist", echo = FALSE------------------------------------------


## ----comment=NA------------------------------------------------------------------


t.test( tcga$PTEN_exp ~ tcga$gender, 
        alternative = c("two.sided"))
 


## ----echo=FALSE, comment=NA------------------------------------------------------
 
ggplot(tcga) + 
  aes(DUOXA1_exp, DUOX1_exp ) + 
  geom_point() + 
  annotate(geom = "text", x = 10, y = 25000, 
           label = paste("r = ", round(r1$r,2), sep = ""),
           color = "red")


## ----echo=FALSE, comment=NA------------------------------------------------------

mod <- lm(DUOX1_exp ~  DUOXA1_exp, data = tcga)

print(summary(mod))



## ----comment=NA------------------------------------------------------------------

mod <- lm(DUOX1_exp ~  DUOXA1_exp, data = tcga)

tidy(mod)



## ----comment=NA------------------------------------------------------------------
 knitr::kable(tidy(mod), format = 'html', digits = 3) %>%
  kable_styling(font_size = 22)



## ----echo=FALSE, comment=NA------------------------------------------------------
 
ggplot(tcga, aes(DUOXA1_exp, DUOX1_exp )) + 
  geom_point() + 
  annotate(geom = "text", x = 10, y = 25000, 
           label = paste("r = ", round(r1$r,2), sep = ""),
           color = "red") + 
   geom_smooth(method = 'lm')



## ----echo=FALSE, comment=NA------------------------------------------------------
knitr::kable(tidy(mod), format = 'html', digits = 3) %>%   kable_styling(font_size = 22)

