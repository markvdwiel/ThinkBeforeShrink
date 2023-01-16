# ThinkBeforeShrink
Annotated R-scripts, Markdowns and Examples to reproduce results in "Think before you shrink..." manuscript.

Repository contains all code to reproduce results in the manuscript: "Think before you shrink: Alternatives to default shrinkage methods can improve prediction accuracy, calibration and coverage" by Mark A. van de Wiel, Gwenaël G.R. Leday, Jeroen Hoogland, Martijn W. Heymans, Erik W. van Zwet, Ailko H. Zwinderman. In order of how the data analyses appear in the manuscript:

1. LinRegrSimIntroExamplePublic: Introductory simulation example for linear regression, that illustrates differential shrinkage of treatment coefficient w.r.t. six other coefficients.
2. HeliusSubsetsPublic: Subset analysis of synthetic Helius data. 
3. LogisticSimPublic: External simulation for binary response setting. 

All folders contain: Annotated R-script for users, R source script with functions and R markdown output (.html) with code and results. In addition, the HeliusSubsetsPublic folder contains the synthetic Helius data.

In addition, the repository contains two simple data analysis examples for linear regression and logistic regression. These focus on usage of the R-packages mgcv, shrinkage, and R-stan for the purpose of fitting low-dimensional problems. We consider variations on ridge regression, including the use (and automatic estimation) of multiple penalties and Bayesian local shrinkage. 
 
