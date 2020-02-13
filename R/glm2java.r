#' Build a Java class of a GLM model, creating a predict() Java method.
#' the dataset through shiny (more control than 
#' \code{expore_all_data()} provides).
#' 
#' @param modCoefs A data.frame with coefficient names and values. The intercept must be passed as containing case sensitive text "Intercept".  E.g. (Intercept) as from glmnet is fine.
#' @param modType A string of the family distribution (gaussian, poisson, binomial, multinomial).
#' @param filePath Where to write the Java class too.
#' @param package The name of the Java class package.
#' @param threshold For binary classification, the probability threshold to use (default=0.5).
#' @examples
#' data(iris)
#' library(glmnet)
#' 
#' ### Fit Gaussian model ###
#' x = as.matrix(iris[,2:4])
#' y = iris[,1]
#' normFit=cv.glmnet(x,y, family="gaussian")
#' normCoefs = as.data.frame(t(as.matrix(coef(normFit, s="lambda.min"))))
#' 
#' buildGLMJavaClass(modCoefs=normCoefs, modType="gaussian", 
#' 		filePath="./NormalMod.java", package="glm2java", addMainMeth=TRUE)
#' 
#' 
#' ### Fit Poisson model ###
#' poisFit=cv.glmnet(x,y, family="poisson")
#' poisCoefs = as.data.frame(t(as.matrix(coef(poisFit, s="lambda.min"))))
#' 
#' ## With the default naive threshold (0.5)
#' buildGLMJavaClass(modCoefs=poisCoefs, modType="poisson", 
#' 		filePath="./PoissonMod.java", package="glm2java", addMainMeth=TRUE)
#' 
#' 
#' ### Fit binomial model ###
#' x = as.matrix(iris[,1:4])
#' y = iris[,5]
#' y = as.character(y)
#' y[y=="virginica"]="versicolor"
#' y=as.factor(y)
#' 
#' # With the default naive threshold (0.5)
#' buildGLMJavaClass(modCoefs=binCoefs, modType="binomial", 
#' 		filePath="./BinaryMod.java", package="glm2java", addMainMeth=TRUE)
#' 
#' 
#' # With a user provided threshold (0.2)
#' buildGLMJavaClass(modCoefs=binCoefs, modType="binomial", 
#' 		filePath="./BinaryMod_2.java", package="glm2java", threshold=0.2, addMainMeth=TRUE)
#' 
#' 
#' 
#' ### Fit multinomial model ###
#' x = as.matrix(iris[,1:4])
#' y = iris[,5]
#' multiFit = cv.glmnet(x,y, family="multinomial" )
#' multiCoefs = lapply(coef(multiFit, s="lambda.min"), 
#'   FUN=function (x) as.data.frame(t(as.matrix(x))))
#' 
#' buildGLMJavaClass(modCoefs=multiCoefs, modType="multinomial", 
#' filePath="./MultinomialMod.java", package="glm2java", addMainMeth=TRUE)
#' 
#' @export
buildGLMJavaClass <- function (modCoefs, modType="gaussian", filePath, package, threshold=0.5, addMainMeth=FALSE) {
	
	require(stringr)
	require(dplyr)
	
	
	clssNm = str_split(filePath,"/")
	clssNm = str_split(clssNm[[1]][length(clssNm[[1]])],"\\.")[[1]][1]
	
	if (!(modType %in% c("gaussian", "binomial", "poisson", "multinomial"))) {
		cat('Not an acceptable type.  Use "gaussian", "binomial", "poisson", or "multinomial".')
		return()
	}
	
	##Set up the model coefficients as a list  
	if (is.data.frame(modCoefs)) modCoefs=list(target=modCoefs)
	
	labels=names(modCoefs)
	
	
	## Set up the file
	if (file.exists(filePath)) file.remove(filePath)
	
	
	#######################################################
	## Write the package and imports and description
	#######################################################
	
	## write package
	cat(paste("package ",package,";\n\n", sep="") ,file=filePath,append=F)
	
	## write imports  - only used for main function example
	if (addMainMeth) {
		cat("import java.sql.*;\n",file=filePath,append=T)
		cat("import java.util.Arrays;\n",file=filePath,append=T)
	}
	
	cat("\n\n/**Functions to do prediction from GLM*/\n\n",file=filePath,append=T)
	
	cat(paste("public class ",clssNm," {\n\n", sep=""),file=filePath,append=T) 
	
	## Print out the methods description
	string=paste("\t/**\n",
			"\t*  Scores the model for an individual observation (i.e., a 1D array of data for a single claim).\n",
			"\t*  Note: the data must be in the proper order as what the model expects:\n",sep="")
	
	nms=names(modCoefs[[1]] %>% select(-contains("Intercept"), -contains("intercept")))
	nmsString=paste("\t*  ",strwrap(paste(nms, collapse=", ",sep=""),width=80),sep="",collapse="\n")
	string=c(string, paste(nmsString,"\n"))
	string=c(string, 
			paste( 
					"\t*  @param data a vector of data inputs (in correct order listed above)\n",
					"\t*  @param preds an empty vector of length 2 for continuous target and logistic,\n",
					"\t*    or length 1+(number of target classes) for multinomial\n",
					"\t*  @return the predicted model score\n",
					"\t*/\n", sep=""))
	cat(string,file=filePath,append=T)

	######################################################
	## Predicted class names callback
	######################################################
	cat('\tprivate static final String[] classNames  = {""', file=filePath, append=T)
	for (i in 1:length(labels)) {
		cat(paste(' ,"',labels[i],'"',sep=""), file=filePath, append=T)
	}
	cat("};\n", file=filePath, append=T)
	cat("\tpublic static String[] getClassNames(){\n", file=filePath, append=T)
	cat("\t\treturn classNames;\n\t}\n", file=filePath, append=T)
	
	#####################################################
	## Feature names callback
	#####################################################
	cat("\tprivate static final String[] featureNames = {\n", file=filePath, append=T)
	nms = names(modCoefs[[1]])
	cat(paste('\t\t"', nms[2],'"', sep=""), file=filePath, append=T)
	for (i in 3:length(nms)) {
		cat(paste(',\n\t\t"',nms[i],'"', sep=""), file=filePath, append=T)
	}
	cat("\n\t};\n", file=filePath, append=T)
	cat("\tpublic static String[] getFeatureNames(){\n", file=filePath, append=T)
	cat("\t\treturn featureNames;\n\t}\n", file=filePath, append=T)
	
	
	
	
	#######################################################
	## Entry prediction method
	#######################################################
	cat("\tpublic static float[] predict( double[] data, float[] preds) {\n",file=filePath,append=T) 
	cat("\n\t\t//fill in empty preds[] vector\n",file=filePath,append=T) 
	cat("\t\tjava.util.Arrays.fill(preds,0f);\n\n",file=filePath,append=T) 
	
	## Set up class-specific methods call (multiple if multinomial, single if logistic) 
	for (i in 1:length(modCoefs)) {
		cat(paste("\t\t//get linear predictor for ",labels[i],"\n",sep=""),file=filePath,append=T) 
		cat(paste("\t\tpreds[",i,"] += ",labels[i],"Predict(data);\n",sep=""),file=filePath,sep="",append=T)
	}
	cat("\n",file=filePath,sep="",append=T)	
	
	switch(modType,
			
			gaussian = {
				cat("\t\t//Compute response from linear predictor (link fnx=unity)\n",file=filePath,append=T) 
				cat("\t\tpreds[0] = preds[1];\n",file=filePath,append=T) 
			},
			poisson = {
				cat("\t\t//Compute response from linear predictor (link fnx=log)\n",file=filePath,append=T) 
				cat("\t\tpreds[0] = (float) (Math.exp(preds[1]));\n",file=filePath,append=T) 
			},
			binomial = {
				cat("\t\t//Compute response from linear predictor (link fnx=logit)\n",file=filePath,append=T) 
				cat("\t\tpreds[1] = (float) (Math.exp(preds[1])/(1+Math.exp(preds[1])));\n",file=filePath,append=T) 
				cat("\n\t\t//Compute the predicted class based on 0.5 threshold\n",file=filePath,append=T) 
				cat(paste("\t\tpreds[0] = preds[1] < ",threshold," ? 0 : 1;\n",sep=""),file=filePath,append=T) 
			},
			multinomial = {
				# Here, link is:
				# exp(linpred_k)/(exp(linpred_1)+exp(linpred_2)+...+exp(linpred_n))
				cat("\t\t//Compute response from linear predictor for each category\n",file=filePath,append=T) 
				cat("\t\t//Inverse link=exp(linpred_k)/(exp(linpred_1)+exp(linpred_2)+...+exp(linpred_n))\n",file=filePath,append=T) 
				cat("\t\tfloat dsum = 0, maxval = Float.NEGATIVE_INFINITY;\n",file=filePath,append=T)
				cat("\t\tfor(int i=1; i<preds.length; i++) maxval = Math.max(maxval, preds[i]);\n",file=filePath,append=T)
				cat("\t\tfor(int i=1; i<preds.length; i++) dsum += (preds[i]=(float) Math.exp(preds[i] - maxval));\n",file=filePath,append=T)
				cat("\t\tfor(int i=1; i<preds.length; i++) preds[i] = preds[i] / dsum;\n",file=filePath,append=T)
				
				cat("\n\t\t//Compute the predicted class based on highest probability\n",file=filePath,append=T)
				cat("\t\tpreds[0] = getPrediction(preds,data);\n",file=filePath,append=T) 
			}
	)
	
	cat("\t\treturn preds;\n",file=filePath,append=T) 
	cat("\t}\n\n",file=filePath,append=T) 
	
	
	#######################################################
	## Set up the individual prediction for linear predictors
	#######################################################
	
	cat("\n\t//Get linear predictors\n",file=filePath,append=T) 
	for (i in 1:length(modCoefs)) {
		cat(paste("\tprivate static float ", labels[i], "Predict(double[] data) {\n",sep=""),file=filePath,append=T) 
		
		intercept = 0
		
		if(!any(grep("[Ii]ntercept|^$", names(modCoefs[[i]][1])))) {
			stop('first element of coefficients expected to be intercept or blank')
		}
		if(any(grep("[Ii]ntercept", names(modCoefs[[i]][1])))) {
			intercept = modCoefs[[i]][1]
		}
			
		linPredString = paste("float pred = (float) (", intercept,"\t//Intercept\n", sep="")
		
		for (j in 2:length(modCoefs[[i]])) {
			linPredString= paste(linPredString,"\t\t\t + data[",j-2,"]*",unlist(modCoefs[[i]][j]),sep="")
			linPredString= paste(linPredString,"\t//",names(modCoefs[[i]])[j] ,"\n",sep="")
		}
		
		cat(paste("\t\t",linPredString,"\t\t\t);\n",sep=""),file=filePath,append=T) 
		
		cat("\t\treturn pred;\n",file=filePath,append=T) 
		cat("\t}\n\n",file=filePath,append=T) 
	}
	
	
	#######################################################
	## Spit out function to get class prediction based on ranking
	#######################################################
	if (modType=="multinomial") addUtilFnx(filePath)
	
	
	## add a main method if requested
	if (addMainMeth) addMainMethod(modCoefs, clssNm, nms, filePath)
	
	cat("}", file=filePath,append=T)
}



addMainMethod <- function(modCoefs, clssNm, nms, filePath) {
	
	
	
	#######################################################
	# Set up main function
	#######################################################
	
	simpleLower <- function(x) {
		s <- strsplit(x, " ")[[1]]
		paste(tolower(substring(s, 1,1)), substring(s, 2),
				sep="", collapse=" ")
	}
	
	cat("\n\t//Example main method reading from DB and setting up variables in proper order\n",file=filePath, append=T)
	cat("\tpublic static void main(String[] args) {\n",file=filePath, append=T)
	cat("\t\tConnection c = null;\n",file=filePath, append=T)
	cat("\t\tStatement stmt = null;\n",file=filePath, append=T)
	
	cat("\n\t\t//Get the models\n",file=filePath, append=T) 
	cat(paste("\t\t",clssNm," model = new ",clssNm,"(); //get model",sep=""),file=filePath, append=T)
	cat('\n\t\tString dbPath = "data/data.db"; //data location\n',file=filePath, append=T)	
	
	cat("\n\t\ttry {",file=filePath, append=T)
	
	cat('\n\t\t\tClass.forName("org.sqlite.JDBC");\n', file=filePath, append=T)
	cat('\t\t\tc = DriverManager.getConnection("jdbc:sqlite:"+dbPath);\n', file=filePath, append=T)
	cat('\t\t\tc.setAutoCommit(false);\n', file=filePath, append=T)
	cat('\t\t\tstmt = c.createStatement();\n', file=filePath, append=T)
	
	cat('\n\t\t\tResultSet rs = stmt.executeQuery( "SELECT * FROM dataTable;" );\n', file=filePath, append=T)
	cat('\t\t\twhile ( rs.next() ) {\n', file=filePath, append=T)
	
	# Get the database iterator results for each row, set up variable names all as double
	for (nm in nms)
		cat(paste('\t\t\t\tfloat ',simpleLower(nm),' = rs.getFloat("',nm,'");\n', sep=""),file=filePath, append=T)
	
	# get a string of the predictions
	cat(paste('\n\t\t\t\tfloat[] preds = model.predict(new double[]{\n',
					paste('\t\t\t\t\t\t',unlist(lapply(nms, simpleLower)),sep="",collapse=",\n"),
					' },',
					'\n\t\t\t\t\t\tnew float[',
					length(modCoefs)+1,
					']);\n', sep=""),file=filePath, append=T)
	
	cat("\n\n\t\t\t\t//Do something with the predictions\n",file=filePath, append=T)
	cat("\n\t\t\t\t//Spit out to console as example\n",file=filePath, append=T)
	cat('\t\t\t\tString str = Arrays.toString(preds);\n', file=filePath, append=T)
	cat('\t\t\t\tstr = str.replaceAll("\\\\[", "").replaceAll("\\\\]","");\n', file=filePath, append=T)
	cat('\t\t\t\tSystem.out.println(str);\n', file=filePath, append=T)
	
	cat('\t\t\t}\n', file=filePath, append=T)
	cat('\t\t\trs.close();\n', file=filePath, append=T)
	cat('\t\t\tstmt.close();\n', file=filePath, append=T)
	cat('\t\t\tc.close();\n', file=filePath, append=T)
	
	cat('\t\t} catch ( Exception e ) {\n', file=filePath, append=T)
	cat('\t\t\tSystem.err.println( e.getClass().getName() + ": " + e.getMessage() );\n', file=filePath, append=T)
	cat('\t\t\tSystem.exit(0);\n', file=filePath, append=T)
	cat('\t\t}\n', file=filePath, append=T)
	cat('\t\tSystem.out.println("Operation done successfully");\n', file=filePath, append=T)
	cat('\t}\n', file=filePath, append=T)
	
	
	# class ended previously
	
	
	
}



## From 0xdata h2o water.util.ModelUtils
addUtilFnx <- function(filePath){
	
	string=paste("\t/**\n",
			"\t*  From 0xdata h2o water.util.ModelUtils\n",
			"\t*  Utility function to get a best prediction from an array of class\n",
			"\t*  prediction distribution.  It returns index of max value if predicted\n",
			"\t*  values are unique.  In the case of tie, the implementation solve it in\n",
			"\t*  pseudo-random way.\n",
			"\t*  @param preds an array of prediction distribution.  Length of arrays is equal to a number of classes+1.\n",
			"\t*  @return the best prediction (index of class, zero-based)\n",
			"\t*/\n",
			"\tprivate static int getPrediction( float[] preds, double data[] ) {\n",
			"\t\tint best=1, tieCnt=0;   // Best class; count of ties\n",
			"\t\tfor( int c=2; c<preds.length; c++) {\n",
			"\t\t\tif( preds[best] < preds[c] ) {\n",
			"\t\t\t\tbest = c;               // take the max index\n",
			"\t\t\t\ttieCnt=0;               // No ties\n",
			"\t\t\t} else if (preds[best] == preds[c]) {\n",
			"\t\t\t\ttieCnt++;               // Ties\n",
			"\t\t\t}\n",
			"\t\t}\n\n",
			"\t\tif( tieCnt==0 ) return best-1; // Return zero-based best class\n",
			"\t\t// Tie-breaking logic\n",
			"\t\tfloat res = preds[best];    // One of the tied best results\n",
			"\t\tlong hash = 0;              // hash for tie-breaking\n",
			"\t\tif( data != null )\n",
			"\t\t\tfor( double d : data ) hash ^= Double.doubleToRawLongBits(d) >> 6; // drop 6 least significants bits of mantisa (layout of long is: 1b sign, 11b exp, 52b mantisa)\n",
			"\t\tint idx = (int)hash%(tieCnt+1);  // Which of the ties we'd like to keep\n",
			"\t\tfor( best=1; best<preds.length; best++)\n",
			"\t\t\tif( res == preds[best] && --idx < 0 )\n",
			"\t\t\t\treturn best-1;          // Return best\n",
			'\t\tthrow new RuntimeException("do not call");\n',
			"\t}\n\n",sep="")
	
	cat(string,file=filePath,append=T) 
	
	
}
