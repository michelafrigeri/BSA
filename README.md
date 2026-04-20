# Bayesian Source Apportionment of Spatio-Temporal Air Pollution Data
### _Michela Frigeri, Veronica J. Berrocal, Alessandra Guglielmi_

## Table of Contents
- [Abstract](#abstract)
- [Prerequisites](#prerequisites)
- [Simulation study](#simulation-study)
- [PM<sub>2.5</sub> speciation data](#applied-study)

## Abstract
Understanding the sources that contribute to fine particulate matter (PM<sub>2.5</sub>)is of crucial importance for designing and implementing targeted air pollution mitigation strategies. Determining what factors contribute to a pollutant's concentration goes under the name of source apportionment and it is a problem long studied by atmospheric scientists and statisticians alike.
In this paper, we propose a Bayesian model for source apportionment that advances the relevant literature by allowing estimation of the number of sources while simultaneously accounting for spatial and temporal dependence in the observed pollutants' concentrations. Taking as example observations of six species of fine particulate matter observed over the course of a year, we present a latent functional factor model that assumes that the marginal distribution of the space-time varying concentrations of each of the six pollutants is lognormal with median expressed as a linear combination of space-time varying emissions produced by an unknown number of sources, each weighted by the corresponding source's composition profile. Estimation of the number of sources is achieved by introducing source-specific shrinkage parameters.
Application of the model to simulated data showcases its ability to retrieve the true number of sources and to estimate the functional latent factors, whereas application to PM<sub>2.5</sub> speciation data in California identifies three major sources.

## Prerequisites
Download daily concentration data from the US Environmental Protection Agency (EPA) 
outdoor air quality web portal AirData (available [here](https://aqs.epa.gov/aqsweb/airdata/download_files.html)).
Downloading the pre-generated daily file of speciated PM<sub>2.5</sub> concentrations for the year 2021, we retain only the
observations relative to California and the six major components of fine particulate matter:
aluminum, organic carbon, elemental carbon, sulfur, sulfate and nitrate. 

## Simulation study
- <ins>sim_data.R</ins> : R code generating the synthetic data used for the simulation study;
- <ins>bsa_sim.stan</ins> : _Stan_ implementation of the spatio-temporal BSA model;
- <ins>sim_stanMCMC.R</ins> : R code running the _Stan_ code on the simulated data and saving the MCMC posterior samples;
- <ins>sim_posterior.R</ins> : R code providing posterior inference for the simulation study.

## Applied study
### PM<sub>2.5</sub> Speciation data
- <ins>EDA_PM25.R</ins> and <ins>PM25_data.R</ins> : R code with explorative data analysis and pre-processing of PM<sub>2.5</sub> daily speciation data (available [here](https://aqs.epa.gov/aqsweb/airdata/download_files.html); see [**Prerequisites**](#prerequisites)) used in our study;
- <ins>bsa_PM25.stan</ins> : _Stan_ implementation of the spatio-temporal BSA model;
- <ins>PM25_stanMCMC.R</ins> : R code running the _Stan_ code on PM<sub>2.5</sub> speciation data and saving the MCMC posterior samples;
- <ins>PM25_posterior.R</ins> : R code providing posterior inference for Bayesian source apportionment of PM<sub>2.5</sub> species.
