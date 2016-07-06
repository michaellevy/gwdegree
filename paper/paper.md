---
title: 'gwdegree: A Shiny app to aid interpretation of geometrically-weighted degree
  estimates in exponential random graph models'
bibliography: paper.bib
date: "05 July 2016"
output: pdf_document
tags:
- social network analysis
- ERGM
- R
authors:
- affiliation: University of California, Davis; Department of Environmental Science
    and Policy
  name: Michael A Levy
  orcid: 0000-0002-4188-2527
---

# Summary

Exponential random graph models (ERGMs) are maximum entropy statistical models that provide estimates on network tie formation of variables both exogenous (covariate) and endogenous (structural) to a network. Network centralization -- the tendency for edges to accrue among a small number of popular nodes -- is a key network variable in many fields, and in ERGMs it is primarily modeled via the geometrically-weighted degree (GWD) statistic [@snijders_new_2006; @hunter_curved_2007]. However, the published literature is ambiguous about how to interpret GWD estimates, and there is little guidance on how to interpret or fix values of the GWD shape-parameter, $\theta_S$. This Shiny application seeks to relieve this ambiguity by demonstrating:

1. how the GWD statistic responds to adding edges to nodes of various degrees, contingent on the value of the shape parameter, $\theta_S$;

1. how the degree distribution of networks of various size and density are shaped by GWD parameter and $\theta_S$ values;

1. how GWD and GWESP -- an ERGM term used to model triadic closure -- interact to affect network centralization and clustering.

The application is bundled as an R package and can be launched by installing and attaching the `gwdegree` package and running the `gwdegree()` function.

# References
