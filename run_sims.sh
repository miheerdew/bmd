#!/bin/bash
rm -r sims/plots/*
Rscript sims.R
Rscript sims_run_methods.R
Rscript simplots.R
