## Geometrically Weighted Degree

There is ambiguity and confusion in the research community about how to interpret GWDegree estimates in exponential random graph models (ERGMs). This app aims to help.

There is a working version [online](michaellevy.shinyapps.io/gwdegree), but it has pretty conservative limits. If you want to do anything more than play with it, you can install and run it locally with the code below, which should launch the app on your machine. 

```
if (!"devtools" %in% installed.packages()[, "Package"]) install.packages("devtools")
devtools::install_github("michaellevy/gwdegree")
library(gwdegree)
gwdegree()
```

I presented a poster on this at Political Networks 2016. You can view that [here](https://figshare.com/articles/Interpretation_of_GW-Degree_Estimates_in_ERGMs/3465020). It got an honorable mention for methodological contribution at the conference.