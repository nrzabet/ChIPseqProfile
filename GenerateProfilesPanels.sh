#!/bin/bash

cd latex/eve
pdflatex mathematicalEstimateAllTFsEveLocus.tex
cd ../../

cd latex/BCDEveStripeEnhancer
pdflatex mathematicalEstimateBCDEveStripeEnhancer.tex
cd ../../

cd latex/BCDAccessibleRegions
pdflatex mathematicalEstimateBCDAccessibleRegionsPWMThreshold070AllLoci.tex
cd ../../

cd latex/CADAccessibleRegions
pdflatex mathematicalEstimateCADAccessibleRegionsPWMThreshold070AllLoci.tex
cd ../../

cd latex/GTAccessibleRegions
pdflatex mathematicalEstimateGTAccessibleRegionsPWMThreshold070AllLoci.tex
cd ../../

cd latex/HBAccessibleRegions
pdflatex mathematicalEstimateHBAccessibleRegionsPWMThreshold070AllLoci.tex
cd ../../

cd latex/KRAccessibleRegions
pdflatex mathematicalEstimateKRAccessibleRegionsPWMThreshold070AllLoci.tex
cd ../../

cd latex/GTAccessibleRegionsJASPAR
pdflatex mathematicalEstimateGTAccessibleRegionsJASPARPWMThreshold070AllLoci.tex
cd ../../

cd latex/HBAccessibleRegionsJASPAR
pdflatex mathematicalEstimateHBAccessibleRegionsJASPARPWMThreshold070AllLoci.tex
cd ../../

cd latex/KRAccessibleRegionsJASPAR
pdflatex mathematicalEstimateKRAccessibleRegionsJASPARPWMThreshold070AllLoci.tex
cd ../../
