### polynomial regression   
*just hit me that polynomial regression is another way to model non-linear relationships like this one so thought I'd poke around*    
*but of course this is wrong because I'm assuming data is normally distributed*
  
  ```{r}

m1_p=lm(proportion~poly(month, 3, raw=T), data=PB_gam_con)
summary(m1_p)
ggplot(PB_gam_con, aes(y=proportion, x=month)) +
  geom_point() +
  stat_smooth(method = 'gam', formula = y ~ poly(x, 3, raw=TRUE))+
  ggtitle("PB males (control)")
```

Notes:
  
  * polynomial regression line looks roughly the same as GAM. but I don't think polynomials allow you to determine when rate of change differed but is good to model non-linear relationship  
