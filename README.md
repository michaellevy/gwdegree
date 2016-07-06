## Geometrically Weighted Degree

![travis badge](https://travis-ci.org/michaellevy/gwdegree.svg?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/gwdegree)](https://cran.r-project.org/package=gwdegree)

There is ambiguity and confusion in the research community about how to interpret GWDegree estimates in exponential random graph models (ERGMs). This app aims to help by providing an interactive platform that demonstrates:

1. how the GWD statistic responds to adding edges to nodes of various degrees, contingent on the value of the shape parameter, $\theta_S$;

1. how the degree distribution of networks of various size and density are shaped by GWD parameter and $\theta_S$ values;

1. how GWD and GWESP -- an ERGM term used to model triadic closure -- interact to affect network centralization and clustering.

All three tabs aim to provide intuition on how GWD parameter values and shape parameter values relate to network structures. For the applied researcher trying to decide whether to estimate or fix the shape parameter value, how to choose or interpret the shape parameter value, and how to interpret the GWD parameter value, the second tab, "Parameter & Degree Distribution" may be particularly useful. Adjust the sliders to match the observed network's size and density. To choose a fixed decay parameter value, examine the possible degree distribution shapes for a given decay parameter value. Once an estimate of the GWD parameter is obtained, examine the implication of the parameter estimate on the degree distribution, *ceteris parabus*.

There is a working version [online](michaellevy.shinyapps.io/gwdegree), but it has pretty conservative limits. If you want to do anything more than play with the app, install it with `install.packages("gwdegree")` and launch the app locally with `library(gwdegree); gwdegree()`.

### Related Conference Poster

I presented a poster on this at Political Networks 2016. You can view that [here](https://figshare.com/articles/Interpretation_of_GW-Degree_Estimates_in_ERGMs/3465020). It got an honorable mention for methodological contribution at the conference.
