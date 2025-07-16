# Bayesian Source Apportionment of Spatio-Temporal Air Pollution Data
### _Michela Frigeri, Veronica J. Berrocal, Alessandra Guglielmi_

## Table of Contents
- [Abstract](#abstract)
- [Prerequisites](#prerequisites)
- [Simulation study](#simulation-study)
- [PM2.5 speciation data](#applied-study)

## Abstract
Understanding the sources that contribute to fine particulate matter
(PM2.5) is of crucial importance for designing and implementing targeted
air pollution mitigation strategies. Determining what factors contribute to a
pollutant’s concentration goes under the name of source apportionment and
it is a problem long studied by atmospheric scientists and statisticians alike.
In this paper, we propose a Bayesian model for source apportionment, that
advances the literature on source apportionment by allowing estimation of
the number of sources and accounting for spatial and temporal dependence in
the observed pollutants’ concentrations. Taking as example observations of
six species of fine particulate matter observed over the course of a year, we
present a latent functional factor model that expresses the space-time varying 
observations of log concentrations of the six pollutant as a linear combination 
of space-time varying emissions produced by an unknown number of
sources each multiplied by the corresponding source’s relative contribution to
the pollutant. Estimation of the number of sources is achieved by introducing
source-specific shrinkage parameters. Application of the model to simulated
data showcases its ability to retrieve the true number of sources and to reliably
estimate the functional latent factors, whereas application to PM2.5 speciation data 
in California identifies 3 major sources for the six PM2.5 species.

## Prerequisites
Download daily concentration data from the US Environmental Protection Agency (EPA) 
outdoor air quality web portal AirData (available [here](https://aqs.epa.gov/aqsweb/airdata/download_files.html)).
Downloading the pre-generated daily file of speciated PM2.5 concentrations for the year 2021, we retain only the
observations relative to California and the six major components of fine particulate matter:
aluminum, organic carbon, elemental carbon, sulfur, sulfate and nitrate. 

## Simulation study
- <ins>sim_data.R</ins> : R code generating the synthetic data used for the simulation study;
- <ins>bsa_dirichlet.stan</ins> : _Stan_ implementation of the spatio-temporal BSA model;
- <ins>sim_stanMCMC.R</ins> : R code running the _Stan_ code on the simulated data and saving the MCMC posterior samples;
- <ins>sim_posterior.R</ins> : R code providing posterior inference for the simulation study.

## Applied study
### PM2.5 Speciation data
- <ins>PM25_data.R</ins> : R code with the pre-processing of PM2.5 daily speciation data (available [here](https://aqs.epa.gov/aqsweb/airdata/download_files.html); see [**Prerequisites**](#prerequisites)) used in our study;
- <ins>bsa_PM25_dirichlet.stan</ins> : _Stan_ implementation of the spatio-temporal BSA model;
- <ins>PM25_stanMCMC.R</ins> : R code running the _Stan_ code on PM2.5 speciation data and saving the MCMC posterior samples;
- <ins>PM25_posterior.R</ins> : R code providing posterior inference for Bayesian source apportionment of PM2.5 species.
