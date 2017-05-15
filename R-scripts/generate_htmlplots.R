# Script for the generation of additional HTML plots

# MySQL Configuration
MYSQL_HOST = "localhost"
MYSQL_PORT = 3306
MYSQL_DBNAME = "crosspare"
MYSQL_USER = "crosspare"
MYSQL_PASSWORT = "crosspare"

# Path for all generates plots
PLOT_PATH = "html_figures/"

#################################
# Install and load dependencies #
#################################
if (!require("plotly")) install.packages("plotly")
if (!require("RMySQL")) install.packages("RMySQL")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("stringr")) install.packages("stringr")

library(plotly)
library(RMySQL)
library(ggplot2)
library(RColorBrewer)
library(stringr)

###############################
# Definition of plot function #
###############################

createHtmlPlot = function(metricDB, metricPrint, dataset) {
  print(paste(metricDB, metricPrint, dataset))
  
  # fetch data
  mydb = dbConnect(MySQL(), user=MYSQL_USER, password=MYSQL_PASSWORT, host=MYSQL_HOST, port=MYSQL_PORT, dbname=MYSQL_DBNAME)
  
  if( dataset=="ALL" ) {
    sqlStatement <<- paste("SELECT config, resultsView.", metricDB, " FROM resultsView WHERE configurationName NOT LIKE 'SELECTEDJURECZKO%' AND configurationName NOT LIKE 'FILTERJURECZKO%';", sep="")
  } else {
    sqlStatement <<- paste("SELECT config, resultsView.", metricDB, " FROM resultsView WHERE configurationName LIKE '", dataset, "%';", sep="")
  }
  rs = dbSendQuery(mydb, sqlStatement)
  
  results = fetch(rs, n=-1)
  dbClearResult(rs)
  dbDisconnect(mydb)
  #results[,2] = as.numeric(results[,2])
  if( metricDB=="auc" ) {
    results$config = with(results, reorder(config, auc, FUN=mean))
  } else if( metricDB=="fscore") {
    results$config = with(results, reorder(config, fscore, FUN=mean))
  } else if( metricDB=="gscore") {
    results$config = with(results, reorder(config, gscore, FUN=mean))
  } else if( metricDB=="mcc") {
    results$config = with(results, reorder(config, mcc, FUN=mean))
  } else if( metricDB=="recall") {
    results$config = with(results, reorder(config, recall, FUN=mean))
  } else if( metricDB=="precision") {
    results$config = with(results, reorder(config, precision, FUN=mean))
  } else if( metricDB=="error") {
    results$error = 1-results$error
    results$config = with(results, reorder(config, error, FUN=mean))
  }
  # create plot
  getPalette = colorRampPalette(brewer.pal(11, "RdYlGn"))
  #plotObj = ggplot(data=results, aes(x=reorder(results["config"], results[metricDB], FUN=mean), y=auc, fill=reorder(config, auc, FUN=mean))) +
  plotObj = ggplot(data=results, aes_string(x="config", y=metricDB, fill="config")) +
    ggtitle(paste(metricPrint, "ordered by the mean value for", dataset, "data")) + ylab(metricPrint) + xlab("Approach") + 
    geom_boxplot() + stat_summary(fun.y=mean, geom="point", size=1, shape=9, aes(text="mean")) +
    coord_flip() + scale_fill_manual(values = getPalette(142)) +
    guides(fill=guide_legend(title=NULL)) +
    theme(axis.text=element_text(size=8), axis.title=element_text(size=8), title=element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1))
  if( metricDB=="mcc" ) {
    plotObj = plotObj + scale_y_continuous(limits=c(-1.0,1.0))
  } else {
    plotObj = plotObj + scale_y_continuous(limits=c(0,1.0))
  }
  plotlyObj = plotly_build(plotObj)
  for( i in 1:length(plotlyObj$x$data) ) {
    if(!is.null(plotlyObj$x$data[[i]]) && plotlyObj$x$data[[i]]$type=="scatter") {
      plotlyObj$x$data[[i]]$text = paste("mean value:", str_extract_all(plotlyObj$x$data[[i]]$text, "\\d+\\.*\\d*")[[1]][length(str_extract_all(plotlyObj$x$data[[i]]$text, "\\d+\\.*\\d*")[[1]])])
    }
  }
  htmlwidgets::saveWidget(as.widget(plotlyObj), paste(PLOT_PATH, metricPrint, "_", dataset, ".html", sep=""))
}

###################
# Plot generation #
###################

metricsDB = c("auc","fscore","gscore","mcc", "recall", "precision", "error")
metricsPrint = c("AUC", "F-measure", "G-measure", "MCC", "recall", "precision", "accuracy")
datasets = c("JURECZKO", "MDP", "AEEEM_LDHHWCHU", "RELINK", "NETGENE", "SELECTEDJURECZKO", "FILTERJURECZKO", "ALL")

for( i in 1:length(metricsDB) ) {
  for( dataset in datasets ) {
    createHtmlPlot(metricsDB[i], metricsPrint[i], dataset)
  }
}
