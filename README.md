# glmToJava
R package: takes a coefficient table and converts it to a Java (POJO) file

## Overview
A simple converter to take a GLM coefficient table (gaussian, binomial, poisson, and multinomial) and convert to a pure Java (POJO) file with
no dependencies.  <br />

Note: using R's extensive PMML support in conjunction with Java PMML (https://github.com/jpmml)
provides much more flexibility, but this simple converter was designed to work on any GLM output (R, SAS, H2O)
to ease migration to production -- all you need is a coefficient table in the correct format (see example usage in
R file below).    <br />

This follows the approach in H2O and uses some of the utility code from their exported Java classes, and 
was first created since H2O didn't provide Java output for their GLM model at the time.  While manually typing in GLM coefficients
into a Java class is trivial, this eases the translation for large numbers of predictors (e.g., with an lasso/ridge/elastic net model)
and will automatically set up the correct link functions.  <br />

Compared to wtcooper/glmToJava, this code generates callback methods getClasses() and getFeatureNames() that
define the expected order of the outputs and inputs (respectively).  In particular, this saves time and effort
in that, Rather than building up the features array manually,
one can fetch the values from a dictionary.  This code also handles intercepts differently.
The first column of each element of the coefficients list is expected to be the intercept
with a column name that contains "Intercept" or "intercept" or is blank.  The code will not 
look any further for the intercept, and will give an error if the first column is apparently not the intercept.

## Installation
```R
library(devtools) 
devtools::install_github("twvacek/glmToJava")
```

## Use  <br />

```R
### Fit Gaussian model ###  
x = as.matrix(iris[,2:4])  
y = iris[,1]  
normFit=cv.glmnet(x,y, family="gaussian")  
normCoefs = as.data.frame(t(as.matrix(coef(normFit, s="lambda.min"))))  
  
buildGLMJavaClass(modCoefs=normCoefs, modType="gaussian", filePath="./NormalMod.java", package="glm2java", addMainMeth=TRUE)  
  
  
### Fit Poisson model ###  
poisFit=cv.glmnet(x,y, family="poisson")  
poisCoefs = as.data.frame(t(as.matrix(coef(poisFit, s="lambda.min"))))  
  
# With the default naive threshold (0.5)  
buildGLMJavaClass(modCoefs=poisCoefs, modType="poisson",  filePath="./PoissonMod.java", package="glm2java", addMainMeth=TRUE)  
  
  
### Fit binomial model ###  
x = as.matrix(iris[,1:4])  
y = iris[,5]  
y = as.character(y)  
y[y=="virginica"]="versicolor"  
y=as.factor(y)  

binFit=cv.glmnet(x,y, family="binomial")  
binCoefs = as.data.frame(t(as.matrix(coef(binFit, s="lambda.min"))))  
  
# With the default naive threshold (0.5)  
buildGLMJavaClass(modCoefs=binCoefs, modType="binomial", filePath="./BinaryMod.java", package="glm2java", addMainMeth=TRUE)  
  
  
# With a user provided threshold (0.2)  
buildGLMJavaClass(modCoefs=binCoefs, modType="binomial", filePath="./BinaryMod_2.java", package="glm2java", threshold=0.2, addMainMeth=TRUE)  
 
  
  
### Fit multinomial model ###  
x = as.matrix(iris[,1:4])  
y = iris[,5]  
multiFit = cv.glmnet(x,y, family="multinomial" )  
multiCoefs = lapply(coef(multiFit, s="lambda.min"),   
  FUN=function (x) as.data.frame(t(as.matrix(x))))  
  
buildGLMJavaClass(modCoefs=multiCoefs, modType="multinomial", filePath="./MultinomialMod.java", package="glm2java", addMainMeth=TRUE)  
```
