# glmToJava
R package: takes a coefficient table and converts it to a Java (POJO) file

## Overview
A simple converter to take a GLM coefficient table (normal, binomial, poisson, and multinomial) and convert to a pure Java (POJO) file with
no dependencies.  Note using R's extensive PMML support in conjunction with Java PMML (https://github.com/jpmml)
provides much more flexibility, but this simple converter was designed to work on any GLM output (R, SAS, H2O)
to ease migration to production -- all you need is a coefficient table in the correct format (see example usage in
R file).    <br />

This follows the approach in H2O and uses some of the utility code from their exported Java classes, and 
was first created since H2O didn't provide Java output for their GLM model.  While manually typing in GLM coefficients
into a Java class is trivial, this eases the translation for large numbers of predictors (e.g., with an lasso/ridge/elastic net model)
and will automatically set up the correct link functions.  <br />

## Installation
  > library(devtools)  <br />
  > devtools::install_github("wtcooper/glmToJava")  <br />
  
## Use  <br />

  > -### Fit Gaussian model ###  <br />
  > x = as.matrix(iris[,2:4])  <br />
  > y = iris[,1]  <br />
  > normFit=cv.glmnet(x,y, family="gaussian")  <br />
  > normCoefs = as.data.frame(t(as.matrix(coef(normFit, s="lambda.min"))))  <br />
  >   <br />
  > buildGLMJavaClass(modCoefs=normCoefs, modType="gaussian", filePath="./NormalMod.java", package="glm2java", addMainMeth=TRUE)  <br />
  >   <br />
  >   <br />
  > -### Fit Poisson model ###  <br />
  > poisFit=cv.glmnet(x,y, family="poisson")  <br />
  > poisCoefs = as.data.frame(t(as.matrix(coef(poisFit, s="lambda.min"))))  <br />
  >  <br /> 
  > -## With the default naive threshold (0.5)  <br />
  > buildGLMJavaClass(modCoefs=poisCoefs, modType="poisson",  filePath="./PoissonMod.java", package="glm2java", addMainMeth=TRUE)  <br />
  >  <br /> 
  >  <br /> 
  > -### Fit binomial model ###  <br />
  > x = as.matrix(iris[,1:4])  <br />
  > y = iris[,5]  <br />
  > y = as.character(y)  <br />
  > y[y=="virginica"]="versicolor"  <br />
  > y=as.factor(y)  <br />
  >  <br /> 
  > -# With the default naive threshold (0.5)  <br />
  > buildGLMJavaClass(modCoefs=binCoefs, modType="binomial", filePath="./BinaryMod.java", package="glm2java", addMainMeth=TRUE)  <br />
  >  <br /> 
  >  <br /> 
  > -# With a user provided threshold (0.2)  <br />
  > buildGLMJavaClass(modCoefs=binCoefs, modType="binomial", filePath="./BinaryMod_2.java", package="glm2java", threshold=0.2, addMainMeth=TRUE)  <br />
  >  <br />
  >   <br />
  >   <br />
  > -### Fit multinomial model ###  <br />
  > x = as.matrix(iris[,1:4])  <br />
  > y = iris[,5]  <br />
  > multiFit = cv.glmnet(x,y, family="multinomial" )  <br />
  > multiCoefs = lapply(coef(multiFit, s="lambda.min"),  <br /> 
  >   FUN=function (x) as.data.frame(t(as.matrix(x))))  <br />
  >  <br /> 
  > buildGLMJavaClass(modCoefs=multiCoefs, modType="multinomial", filePath="./MultinomialMod.java", package="glm2java", addMainMeth=TRUE)  <br />
