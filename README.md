# Geometrically Weighted Degree

![travis badge](https://travis-ci.org/michaellevy/gwdegree.svg?branch=master) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/gwdegree)](https://cran.r-project.org/package=gwdegree) [![status](http://joss.theoj.org/papers/0f4eda35180f77176ce495cd4d711075/status.svg)](http://joss.theoj.org/papers/0f4eda35180f77176ce495cd4d711075)

### Overview

There is ambiguity and confusion in the research community about how to interpret GWDegree estimates in exponential random graph models (ERGMs). This app aims to help by providing an interactive platform that demonstrates:

1. how the GWD statistic responds to adding edges to nodes of various degrees, contingent on the value of the shape parameter, $\theta_S$;

1. how the degree distribution of networks of various size and density are shaped by GWD parameter and $\theta_S$ values;

1. how GWD and GWESP -- an ERGM term used to model triadic closure -- interact to affect network centralization and clustering.

### Installation

The application is available as a web-app [online](https://michaellevy.shinyapps.io/gwdegree/), but with conservative simulation limits and limited bandwidth. If you want to do anything more than play with it, install the `gwdegree` package for R from CRAN via `install.packages("gwdegree")`.

### Use

To run the application from R, first attach the package with `library(gwdegree)`, then launch it with `gwdegree()`.

The three functionalities described above are split into separate tabs in the application. Researchers trying to decide whether to estimate or fix the shape parameter value, how to choose or interpret the shape parameter value, or how to interpret the GWD parameter value may find the second tab, "Parameter & Degree Distribution", especially useful by setting the network size and density to match the observed network and exploring the influence of the two parameters on the shape of the degree distribution.

### Feedback, questions, and contributions

If you would like to contribute to this app, please submit a pull request to the [GitHub repository](https://github.com/michaellevy/gwdegree). Feel free to report issues, provide feedback, suggest improvements, and ask questions by opening a [new issue](https://github.com/michaellevy/gwdegree/issues).


## Related Conference Poster

I presented a poster on this at the 2016 Political Networks conference. You can view that [here](https://figshare.com/articles/Interpretation_of_GW-Degree_Estimates_in_ERGMs/3465020). It won an honorable mention for methodological contribution at the conference.
