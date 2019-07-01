rm(list=ls())
library(shiny)
library(ggplot2)
library(stringr)

source("waterfall_plot_methods.r", local = TRUE)
source("common.r", local = TRUE)

load("cancerData.RData")

