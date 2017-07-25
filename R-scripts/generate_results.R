# Script for the result evaluation and generation of plots. 

#########################
# Definition of Globals #
#########################

# You may change these global variables for different purposes.
# You should at least check the MySQL configuration below

# MySQL Configuration -
MYSQL_HOST = "localhost"
MYSQL_PORT = 3306
MYSQL_DBNAME = "crosspare"
MYSQL_USER = "crosspare"
MYSQL_PASSWORT = "crosspare"

# Path for all generates plots
PLOT_PATH = "figures/"

# Path for all generated tables
TABLES_PATH = "tables/"

# If true the plots included in the article are generated
CREATEPLOTS = FALSE

# Defines wether Friedman-Nemenyi or ANOVA with Scott-Knott clustering is used for ranking
# TRUE = Friedman-Nemenyi, FALSE = Scott-Knott
NONPARAMETRIC = TRUE

# 1 for fine-grained Nemenyi ranking from the correction
# Any other number for the three-ranks approach from the original paper
NONPARAMETRICRANKINGMODE = 1

# If true CD charts for Friedman-Nemenyi tests are created
PRINTCDCHARTS = FALSE
CD_EXPORT_WIDTH = 10

#################################
# Install and load dependencies #
#################################
if (!require("RMySQL")) install.packages("RMySQL")
if (!require("ScottKnott")) install.packages("ScottKnott")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("xtable")) install.packages("xtable")
if (!require("reshape")) install.packages("reshape")
if (!require("gdata")) install.packages("gdata")
if (!require("stringr")) install.packages("stringr")

library(RMySQL)
library(ScottKnott)
library(ggplot2)
library(gridExtra)
library(xtable)
library(reshape)
library(gdata)
library(stringr)

# We have to increase the number of nested expressions, because the CD diagrams cannot be plotted otherwise
options(expressions=10000)

########################################################
# Definition of functions with result evaluation logic #
########################################################
evaluateCPDPBenchmark = function(metricNames, datasets) {
  mydb = dbConnect(MySQL(), user=MYSQL_USER, password=MYSQL_PASSWORT, host=MYSQL_HOST, port=MYSQL_PORT, dbname=MYSQL_DBNAME)
  
  dbListTables(mydb)
  dbListFields(mydb, 'crosspare.results')
  
  sqlStatement = "SELECT DISTINCT concat(substring(configurationName, 8), '-', classifier) as config FROM resultsView WHERE configurationName LIKE 'RELINK-%';"
  rs = dbSendQuery(mydb, sqlStatement)
  configurations <- fetch(rs, n=-1)
  
  for( i in 1:length(datasets) ) {
    dataset = datasets[i]
    sqlStatement = paste("SELECT concat(substring(configurationName, ",
                         nchar(dataset)+2,
                         "), '-', classifier) as config, ",
                         paste(metricNames, collapse=", "),
                         " FROM resultsView WHERE configurationName LIKE '",
                         dataset, "-%';", sep = "")
    rs = dbSendQuery(mydb, sqlStatement)
    results = fetch(rs, n=-1)
    results[is.na(results)] = 0
    
    cnts = data.frame(config=unique(results$config), count=0)
    for( j in 1:nrow(results) ) {
      cnts[cnts$config==results$config[j],]$count = cnts[cnts$config==results$config[j],]$count+1
      results$index[j] = cnts[cnts$config==results$config[j],]$count
    }
    
    if( "auc" %in% metricNames ) {
      if( NONPARAMETRIC) {
        # Friedman-Nemenyi Test
        meanRanks = createMeanRankMat(results, "auc")
        if( PRINTCDCHARTS ) {
          printMetricName = "AUC"
          plotObj = plotCDNemenyiFriedman(meanRanks, max(results$index), title=paste("Critical Distance Diagram for ", printMetricName, " and the ", dataset, " data", sep=""))
          ggsave(filename = paste(PLOT_PATH,printMetricName,"_",dataset,".png",sep=""), plot=plotObj, device="png", width = CD_EXPORT_WIDTH, height = CD_EXPORT_WIDTH)
        }
        
        colnumber = ncol(configurations)+1
        if( NONPARAMETRICRANKINGMODE==1 ) {
          for( j in 1:nrow(configurations) ) {
            if( configurations$config[j] %in% meanRanks$config) {
              configurations[j,colnumber] = meanRanks[meanRanks$config==configurations$config[j],]$normRank
            }
          }
        } else {
          cdNemenyi = getNemenyiCD(0.05, length(unique(results$config)), max(results$index))
          maxMeanRank = max(meanRanks$meanRank)
          minMeanRank = min(meanRanks$meanRank)
          
          # with CD within first results
          if( maxMeanRank-minMeanRank<=cdNemenyi ) {
            # all on first rank
            configurations[, colnumber] = 1
          }
          else if(maxMeanRank-minMeanRank<=2*cdNemenyi)  {
            # only two ranks
            firstRank = meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank = meanRanks$config[meanRanks$meanRank<maxMeanRank-cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
          } else {
            # three ranks
            firstRank <<- meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank <<- meanRanks$config[meanRanks$meanRank>minMeanRank+cdNemenyi & meanRanks$meanRank<maxMeanRank-cdNemenyi]
            thirdRank <<- meanRanks$config[meanRanks$meanRank<=minMeanRank+cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
            configurations[configurations$config %in% thirdRank, colnumber] = 1-(length(firstRank)+length(secondRank))/(nrow(meanRanks)-1)
          }
        }
        colnames(configurations)[colnumber] = paste("AUC",dataset,sep="_")
      } else {
        # AOV and Scott-Knott test
        skresult = SK(aov( auc ~ config, results), sig.level=0.05)
        colnumber = ncol(configurations)+1
        higherRanked = 0
        for( j in 1:length(unique(skresult$groups)) ) {
          orderedNms = skresult$nms[skresult$ord]
          configurations[configurations$config %in% orderedNms[skresult$groups==j], colnumber] = 1-(higherRanked/(nrow(configurations)-1))
          higherRanked <- higherRanked+length(orderedNms[skresult$groups==j])
        }
        colnames(configurations)[colnumber] = paste("AUC",dataset,sep="_")
      }
      colnumber = ncol(configurations)+1
      for( j in 1:nrow(configurations) ) {
        configurations[j,colnumber] = mean(results[results$config==configurations$config[j],]$auc, na.rm = TRUE)
      }
      colnames(configurations)[colnumber] = paste("meanAUC",dataset,sep="_")
    }
    
    if( "fscore" %in% metricNames ) {
      if (NONPARAMETRIC ) {
        # Friedman-Nemenyi Test
        meanRanks = createMeanRankMat(results, "fscore")
        if( PRINTCDCHARTS ) {
          printMetricName = "F-Measure"
          plotObj = plotCDNemenyiFriedman(meanRanks, max(results$index), title=paste("Critical Distance Diagram for ", printMetricName, " and the ", dataset, " data", sep=""))
          ggsave(filename = paste(PLOT_PATH,printMetricName,"_",dataset,".png",sep=""), plot=plotObj, device="png", width = CD_EXPORT_WIDTH, height = CD_EXPORT_WIDTH)
        }
        
        colnumber = ncol(configurations)+1
        if( NONPARAMETRICRANKINGMODE==1 ) {
          for( j in 1:nrow(configurations) ) {
            if( configurations$config[j] %in% meanRanks$config) {
              configurations[j,colnumber] = meanRanks[meanRanks$config==configurations$config[j],]$normRank
            }
          }
        } else {
          cdNemenyi = getNemenyiCD(0.05, length(unique(results$config)), max(results$index))
          maxMeanRank = max(meanRanks$meanRank)
          minMeanRank = min(meanRanks$meanRank)
          
          # with CD within first results
          if( maxMeanRank-minMeanRank<=cdNemenyi ) {
            # all on first rank
            configurations[, colnumber] = 1
          }
          else if(maxMeanRank-minMeanRank<=2*cdNemenyi)  {
            # only two ranks
            firstRank = meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank = meanRanks$config[meanRanks$meanRank<maxMeanRank-cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
          } else {
            # three ranks
            firstRank <<- meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank <<- meanRanks$config[meanRanks$meanRank>minMeanRank+cdNemenyi & meanRanks$meanRank<maxMeanRank-cdNemenyi]
            thirdRank <<- meanRanks$config[meanRanks$meanRank<=minMeanRank+cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
            configurations[configurations$config %in% thirdRank, colnumber] = 1-(length(firstRank)+length(secondRank))/(nrow(meanRanks)-1)
          }
        }
        colnames(configurations)[colnumber] = paste("FMEAS",dataset,sep="_")
      }
      else {
        skresult = SK(aov( fscore ~ config, results), sig.level=0.05)
        colnumber = ncol(configurations)+1
        higherRanked = 0
        for( j in 1:length(unique(skresult$groups)) ) {
          orderedNms = skresult$nms[skresult$ord]
          configurations[configurations$config %in% orderedNms[skresult$groups==j], colnumber] = 1-(higherRanked/(nrow(configurations)-1))
          higherRanked <- higherRanked+length(orderedNms[skresult$groups==j])
        }
        colnames(configurations)[colnumber] = paste("FMEAS",dataset,sep="_")
      }
      colnumber = ncol(configurations)+1
      for( j in 1:nrow(configurations) ) {
        configurations[j,colnumber] = mean(results[results$config==configurations$config[j],]$fscore, na.rm = TRUE)
      }
      colnames(configurations)[colnumber] = paste("meanFMEAS",dataset,sep="_")
    }
    
    if( "gscore" %in% metricNames ) {
      if (NONPARAMETRIC ) {
        # Friedman-Nemenyi Test
        meanRanks = createMeanRankMat(results, "gscore")
        if( PRINTCDCHARTS ) {
          printMetricName = "G-measure"
          plotObj = plotCDNemenyiFriedman(meanRanks, max(results$index), title=paste("Critical Distance Diagram for ", printMetricName, " and the ", dataset, " data", sep=""))
          ggsave(filename = paste(PLOT_PATH,printMetricName,"_",dataset,".png",sep=""), plot=plotObj, device="png", width = CD_EXPORT_WIDTH, height = CD_EXPORT_WIDTH)
        }
        
        colnumber = ncol(configurations)+1
        if( NONPARAMETRICRANKINGMODE==1 ) {
          for( j in 1:nrow(configurations) ) {
            if( configurations$config[j] %in% meanRanks$config) {
              configurations[j,colnumber] = meanRanks[meanRanks$config==configurations$config[j],]$normRank
            }
          }
        } else {
          cdNemenyi = getNemenyiCD(0.05, length(unique(results$config)), max(results$index))
          maxMeanRank = max(meanRanks$meanRank)
          minMeanRank = min(meanRanks$meanRank)
          
          # with CD within first results
          if( maxMeanRank-minMeanRank<=cdNemenyi ) {
            # all on first rank
            configurations[, colnumber] = 1
          }
          else if(maxMeanRank-minMeanRank<=2*cdNemenyi)  {
            # only two ranks
            firstRank = meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank = meanRanks$config[meanRanks$meanRank<maxMeanRank-cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
          } else {
            # three ranks
            firstRank <<- meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank <<- meanRanks$config[meanRanks$meanRank>minMeanRank+cdNemenyi & meanRanks$meanRank<maxMeanRank-cdNemenyi]
            thirdRank <<- meanRanks$config[meanRanks$meanRank<=minMeanRank+cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
            configurations[configurations$config %in% thirdRank, colnumber] = 1-(length(firstRank)+length(secondRank))/(nrow(meanRanks)-1)
          }
        }
        colnames(configurations)[colnumber] = paste("GMEAS",dataset,sep="_")
      }
      else {
        skresult = SK(aov( gscore ~ config, results), sig.level=0.05)
        colnumber = ncol(configurations)+1
        higherRanked = 0
        for( j in 1:length(unique(skresult$groups)) ) {
          orderedNms = skresult$nms[skresult$ord]
          configurations[configurations$config %in% orderedNms[skresult$groups==j], colnumber] = 1-(higherRanked/(nrow(configurations)-1))
          higherRanked <- higherRanked+length(orderedNms[skresult$groups==j])
        }
        colnames(configurations)[colnumber] = paste("GMEAS",dataset,sep="_")
      }
      colnumber = ncol(configurations)+1
      for( j in 1:nrow(configurations) ) {
        configurations[j,colnumber] = mean(results[results$config==configurations$config[j],]$gscore, na.rm = TRUE)
      }
      colnames(configurations)[colnumber] = paste("meanGMEAS",dataset,sep="_")
    }
    
    if( "mcc" %in% metricNames ) {
      if (NONPARAMETRIC ) {
        # Friedman-Nemenyi Test
        meanRanks = createMeanRankMat(results, "mcc")
        if( PRINTCDCHARTS ) {
          printMetricName = "MCC"
          plotObj = plotCDNemenyiFriedman(meanRanks, max(results$index), title=paste("Critical Distance Diagram for ", printMetricName, " and the ", dataset, " data", sep=""))
          ggsave(filename = paste(PLOT_PATH,printMetricName,"_",dataset,".png",sep=""), plot=plotObj, device="png", width = CD_EXPORT_WIDTH, height = CD_EXPORT_WIDTH)
        }
        
        colnumber = ncol(configurations)+1
        if( NONPARAMETRICRANKINGMODE==1 ) {
          for( j in 1:nrow(configurations) ) {
            if( configurations$config[j] %in% meanRanks$config) {
              configurations[j,colnumber] = meanRanks[meanRanks$config==configurations$config[j],]$normRank
            }
          }
        } else {
          cdNemenyi = getNemenyiCD(0.05, length(unique(results$config)), max(results$index))
          maxMeanRank = max(meanRanks$meanRank)
          minMeanRank = min(meanRanks$meanRank)
          
          # with CD within first results
          if( maxMeanRank-minMeanRank<=cdNemenyi ) {
            # all on first rank
            configurations[, colnumber] = 1
          }
          else if(maxMeanRank-minMeanRank<=2*cdNemenyi)  {
            # only two ranks
            firstRank = meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank = meanRanks$config[meanRanks$meanRank<maxMeanRank-cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
          } else {
            # three ranks
            firstRank <<- meanRanks$config[meanRanks$meanRank>=maxMeanRank-cdNemenyi]
            secondRank <<- meanRanks$config[meanRanks$meanRank>minMeanRank+cdNemenyi & meanRanks$meanRank<maxMeanRank-cdNemenyi]
            thirdRank <<- meanRanks$config[meanRanks$meanRank<=minMeanRank+cdNemenyi]
            configurations[configurations$config %in% firstRank, colnumber] = 1
            configurations[configurations$config %in% secondRank, colnumber] = 1-length(firstRank)/(nrow(meanRanks)-1)
            configurations[configurations$config %in% thirdRank, colnumber] = 1-(length(firstRank)+length(secondRank))/(nrow(meanRanks)-1)
          }
        }
        colnames(configurations)[colnumber] = paste("MCC",dataset,sep="_")
      }
      else {
        skresult = SK(aov( mcc ~ config, results), sig.level=0.05)
        colnumber = ncol(configurations)+1
        higherRanked = 0
        for( j in 1:length(unique(skresult$groups)) ) {
          orderedNms = skresult$nms[skresult$ord]
          configurations[configurations$config %in% orderedNms[skresult$groups==j], colnumber] = 1-(higherRanked/(nrow(configurations)-1))
          higherRanked <- higherRanked+length(orderedNms[skresult$groups==j])
        }
        colnames(configurations)[colnumber] = paste("MCC",dataset,sep="_")
      }
      colnumber = ncol(configurations)+1
      for( j in 1:nrow(configurations) ) {
        configurations[j,colnumber] = mean(results[results$config==configurations$config[j],]$mcc, na.rm = TRUE)
      }
      colnames(configurations)[colnumber] = paste("meanMCC",dataset,sep="_")
    }
    
    if( CREATEPLOTS ) {
      if( "auc" %in% metricNames ) {
        plotauc = ggplot(data=results, aes(reorder(config, auc, FUN=mean), auc)) +
          ggtitle(paste("AUC ordered by the mean value for", dataset, "data")) +
          geom_boxplot() + stat_summary(fun.y=mean, geom="point", size=3, shape=9) +
          scale_y_continuous(limits=c(0,1.0)) + coord_flip() +
          theme(axis.text=element_text(size=8), axis.title=element_text(size=8), title=element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(filename=paste(PLOT_PATH,"AUC_",dataset,".pdf",sep=""), plot=plotauc)
      }
      
      if( "fscore" %in% metricNames ) {
        plotfscore = ggplot(data=results, aes(reorder(config, fscore, FUN=mean), fscore)) +
          ggtitle(paste("F-measure ordered by the mean value for", dataset, "data")) +
          geom_boxplot() + stat_summary(fun.y=mean, geom="point", size=3, shape=9) +
          scale_y_continuous(limits=c(0,1.0)) + coord_flip() +
          theme(axis.text=element_text(size=8), axis.title=element_text(size=8), title=element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(filename=paste(PLOT_PATH,"FMEASURE_",dataset,".pdf",sep=""), plot=plotfscore) 
      }
      
      if( "gscore" %in% metricNames ) {
        plotgscore = ggplot(data=results, aes(reorder(config, gscore, FUN=mean), gscore)) +
          ggtitle(paste("G-measure ordered by the mean value for", dataset, "data")) +
          geom_boxplot() + stat_summary(fun.y=mean, geom="point", size=3, shape=9) +
          scale_y_continuous(limits=c(0,1.0)) + coord_flip() +
          theme(axis.text=element_text(size=8), axis.title=element_text(size=8), title=element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(filename=paste(PLOT_PATH,"GMEASURE_",dataset,".pdf",sep=""), plot=plotgscore) 
      }
      
      if( "mcc" %in% metricNames ) {
        plotmcc = ggplot(data=results, aes(reorder(config, mcc, FUN=mean), mcc)) +
          ggtitle(paste("MCC ordered by the mean value for", dataset, "data")) +
          geom_boxplot() + stat_summary(fun.y=mean, geom="point", size=3, shape=9) +
          scale_y_continuous(limits=c(0,1.0)) + coord_flip() +
          theme(axis.text=element_text(size=8), axis.title=element_text(size=8), title=element_text(size=8), axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(filename=paste(PLOT_PATH,"MCC_",dataset,".pdf",sep=""), plot=plotmcc)
      }
    }
    
  }
  dbDisconnect(mydb)
  
  skresults = configurations
  rownames(skresults) = skresults$config
  skresults$config = NULL
  
  skresults$MEANRANK = rowMeans(skresults[!startsWith(colnames(skresults), "mean")], na.rm=TRUE)
  skresults$config = as.factor(rownames(skresults))
  for( i in 1:nrow(skresults) ) {
    skresults$approach[i] = strsplit(rownames(skresults)[[i]],"-")[[1]][1]
    splitLength = nchar(strsplit(rownames(skresults)[[i]],"-")[[1]][1])
    skresults$classifier[i] = substring(rownames(skresults)[[i]],splitLength+2)
  }
  
  
  return(skresults)
}

plotBestResults = function(results, rq, metricsString) {
  best = results
  approaches = unique(best$approach)
  for( j in 1:length(approaches)) {
    maxval = max(best[best$approach==approaches[j],]$MEANRANK)
    for( i in nrow(best):1 ) {
      if( best$approach[i]==approaches[j] && best$MEANRANK[i]<maxval ) {
        best = best[-i,]
      }
    }
  }
  
  best$MEANRANK = round(best$MEANRANK, digits=3)
  textsize = 12
  plotBest = ggplot(data=best, aes(reorder(config, MEANRANK, FUN=mean), MEANRANK)) +
    ggtitle(paste("Ranking of approaches using",metricsString)) + ylab("Mean rankscore") + xlab("Approach") +
    geom_bar(stat="identity" ) + geom_label(aes(label=MEANRANK), hjust=0, nudge_y = 0.002) +
    scale_y_continuous(limits=c(0,1.02), breaks=c(0,0.25,0.5,0.75,1)) + coord_flip() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), title=element_text(size=textsize))
  print(plotBest)
  ggsave(filename=paste(PLOT_PATH,"BEST_", rq, ".pdf", sep=""), plot=plotBest)
  return(best)
}

writeResultsTableRQ1 = function(rq1best, rq4results, rq5results) {
  resultsTableFrame = round(rq1best[,1:41], digits=2)
  
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanAUC_JURECZKO, " (", resultsTableFrame$AUC_JURECZKO, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "auc"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanFMEAS_JURECZKO, " (", resultsTableFrame$FMEAS_JURECZKO, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "F-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanGMEAS_JURECZKO, " (", resultsTableFrame$GMEAS_JURECZKO, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "G-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanMCC_JURECZKO, " (", resultsTableFrame$MCC_JURECZKO, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "MCC"
  colnumber = ncol(resultsTableFrame)+1
  
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(round(rq4results[rq4results$config %in% rownames(resultsTableFrame),]$meanAUC_FILTERJURECZKO-rq1best$meanAUC_JURECZKO, digits=2), " / " ,round(rq5results[rq5results$config %in% rownames(resultsTableFrame),]$meanAUC_SELECTEDJURECZKO-rq1best$meanAUC_JURECZKO, digits=2), sep="")
  colnames(resultsTableFrame)[colnumber] = "auc"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(round(rq4results[rq4results$config %in% rownames(resultsTableFrame),]$meanFMEAS_FILTERJURECZKO-rq1best$meanFMEAS_JURECZKO, digits=2), " / " ,round(rq5results[rq5results$config %in% rownames(resultsTableFrame),]$meanFMEAS_SELECTEDJURECZKO-rq1best$meanFMEAS_JURECZKO, digits=2), sep="")
  colnames(resultsTableFrame)[colnumber] = "F-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(round(rq4results[rq4results$config %in% rownames(resultsTableFrame),]$meanGMEAS_FILTERJURECZKO-rq1best$meanGMEAS_JURECZKO, digits=2), " / " ,round(rq5results[rq5results$config %in% rownames(resultsTableFrame),]$meanGMEAS_SELECTEDJURECZKO-rq1best$meanGMEAS_JURECZKO, digits=2), sep="")
  colnames(resultsTableFrame)[colnumber] = "G-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(round(rq4results[rq4results$config %in% rownames(resultsTableFrame),]$meanMCC_FILTERJURECZKO-rq1best$meanMCC_JURECZKO, digits=2), " / " ,round(rq5results[rq5results$config %in% rownames(resultsTableFrame),]$meanMCC_SELECTEDJURECZKO-rq1best$meanMCC_JURECZKO, digits=2), sep="")
  colnames(resultsTableFrame)[colnumber] = "MCC"
  colnumber = ncol(resultsTableFrame)+1
  
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanAUC_MDP, " (", resultsTableFrame$AUC_MDP, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "auc"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanFMEAS_MDP, " (", resultsTableFrame$FMEAS_MDP, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "F-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanGMEAS_MDP, " (", resultsTableFrame$GMEAS_MDP, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "G-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanMCC_MDP, " (", resultsTableFrame$MCC_MDP, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "MCC"
  colnumber = ncol(resultsTableFrame)+1
  
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanAUC_AEEEM_LDHHWCHU, " (", resultsTableFrame$AUC_AEEEM_LDHHWCHU, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "auc"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanFMEAS_AEEEM_LDHHWCHU, " (", resultsTableFrame$FMEAS_AEEEM_LDHHWCHU, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "F-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanGMEAS_AEEEM_LDHHWCHU, " (", resultsTableFrame$GMEAS_AEEEM_LDHHWCHU, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "G-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanMCC_AEEEM_LDHHWCHU, " (", resultsTableFrame$MCC_AEEEM_LDHHWCHU, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "MCC"
  colnumber = ncol(resultsTableFrame)+1
  
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanAUC_NETGENE, " (", resultsTableFrame$AUC_NETGENE, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "auc"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanFMEAS_NETGENE, " (", resultsTableFrame$FMEAS_NETGENE, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "F-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanGMEAS_NETGENE, " (", resultsTableFrame$GMEAS_NETGENE, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "G-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanMCC_NETGENE, " (", resultsTableFrame$MCC_NETGENE, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "MCC"
  colnumber = ncol(resultsTableFrame)+1
  
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanAUC_RELINK, " (", resultsTableFrame$AUC_RELINK, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "auc"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanFMEAS_RELINK, " (", resultsTableFrame$FMEAS_RELINK, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "F-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanGMEAS_RELINK, " (", resultsTableFrame$GMEAS_RELINK, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "G-measure"
  colnumber = ncol(resultsTableFrame)+1
  resultsTableFrame[,colnumber] = paste(resultsTableFrame$meanMCC_RELINK, " (", resultsTableFrame$MCC_RELINK, ")", sep="")
  colnames(resultsTableFrame)[colnumber] = "MCC"
  colnumber = ncol(resultsTableFrame)+1
  
  tableCaption = "Mean results over all products with rankscores in brackets. Bold-faced values are top-ranking for the metric on the data set. For FILTERJURECZKO and SELECTEDJURECKO, we show the difference in the mean values to JURECZKO."
  tablePart1 = print.xtable(xtable(resultsTableFrame[,42:53]))
  tablePart2 = print.xtable(xtable(resultsTableFrame[,54:65]))
  
  tablePart1 = sub("\\centering", "", tablePart1, fixed=TRUE)
  tablePart1 = sub("\\begin{table}[ht]", "\\begin{table*}\n\\scriptsize\\centering\\begin{sideways}", tablePart1, fixed=TRUE)
  tablePart1 = sub("rllllrrrrllll", "|r|llll|llll|llll|", tablePart1, fixed=TRUE)
  tablePart1 = sub("\\hline", "\\hline\n & \\multicolumn{4}{c|}{JURECZKO} & \\multicolumn{4}{c|}{FILTERJURECZKO / SELECTEDJURECZKO} & \\multicolumn{4}{c|}{MDP} \\\\\n\\hline", tablePart1, fixed=TRUE)
  tablePart1 = sub("\\end{tabular}\n\\end{table}\n", "", tablePart1, fixed=TRUE)
  tablePart2 = sub("\\begin{table}[ht]\n\\centering\n\\begin{tabular}{rllllllllllll}\n  \\hline\n", "\\hline\n & \\multicolumn{4}{c|}{AEEEM} & \\multicolumn{4}{c|}{NETGENE} & \\multicolumn{4}{c|}{RELINK} \\\\\n\\hline", tablePart2, fixed=TRUE)
  tablePart2 = sub(".*\\\\hline\\\n & \\\\\\multicolumn","\\\\hline\n & \\\\multicolumn",tablePart2, fixed=FALSE)
  tablePart2 = sub("\\end{table}", paste("\\end{sideways}\n\\caption{", tableCaption, "}\n\\label{tbl:results}\n\\end{table*}", sep=""), tablePart2, fixed=TRUE)
  tableStr = paste(tablePart1, tablePart2, sep="")
  tableStr = gsub("NaN (NA)", "-", tableStr, fixed=TRUE)
  tableStr = gsub("NaN / NaN", "-", tableStr, fixed=TRUE)
  tableStr = gsub("&  &", "& - &", tableStr, fixed=TRUE)
  tableStr = gsub("0 ", "0.00 ", tableStr, fixed=TRUE)
  tableStr = gsub("auc.1", "\\emph{AUC}", tableStr, fixed=TRUE)
  tableStr = gsub("auc.2", "\\emph{AUC}", tableStr, fixed=TRUE)
  tableStr = gsub("auc.3", "\\emph{AUC}", tableStr, fixed=TRUE)
  tableStr = gsub("auc.4", "\\emph{AUC}", tableStr, fixed=TRUE)
  tableStr = gsub("auc.5", "\\emph{AUC}", tableStr, fixed=TRUE)
  tableStr = gsub("auc ", "\\emph{AUC} ", tableStr, fixed=TRUE)
  tableStr = gsub("F-measure.1", "\\emph{F-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("F-measure.2", "\\emph{F-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("F-measure.3", "\\emph{F-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("F-measure.4", "\\emph{F-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("F-measure.5", "\\emph{F-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("F-measure ", "\\emph{F-measure} ", tableStr, fixed=TRUE)
  tableStr = gsub("G-measure.1", "\\emph{G-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("G-measure.2", "\\emph{G-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("G-measure.3", "\\emph{G-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("G-measure.4", "\\emph{G-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("G-measure.5", "\\emph{G-measure}", tableStr, fixed=TRUE)
  tableStr = gsub("G-measure ", "\\emph{G-measure} ", tableStr, fixed=TRUE)
  tableStr = gsub("MCC.1", "\\emph{MCC}", tableStr, fixed=TRUE)
  tableStr = gsub("MCC.2", "\\emph{MCC}", tableStr, fixed=TRUE)
  tableStr = gsub("MCC.3", "\\emph{MCC}", tableStr, fixed=TRUE)
  tableStr = gsub("MCC.4", "\\emph{MCC}", tableStr, fixed=TRUE)
  tableStr = gsub("MCC.5", "\\emph{MCC}", tableStr, fixed=TRUE)
  tableStr = gsub("MCC ", "\\emph{MCC} ", tableStr, fixed=TRUE)
  tableStr = gsub("([0-9]\\.[0-9][0-9]? \\(1\\))", "\\\\textbf\\{\\1\\}", tableStr, perl=TRUE )
  write(tableStr, file=paste(TABLES_PATH, "resultsTable.tex", sep=""))
}

compareDatasets = function(set1, set2, metric) {
  mydb = dbConnect(MySQL(), user=MYSQL_USER, password=MYSQL_PASSWORT, host=MYSQL_HOST, port=MYSQL_PORT, dbname=MYSQL_DBNAME)
  
  sqlStatement = paste("SELECT q1.config, q1.mean1, q2.mean2 FROM ",
                       "(SELECT concat(substring(configurationName, ",nchar(set1)+2,"), '-', classifier) as config, avg(", metric,
                       ") as mean1 FROM resultsView WHERE configurationName LIKE '",set1, "-%' GROUP BY configurationName, classifier) q1, ",
                       "(SELECT concat(substring(configurationName, ",nchar(set2)+2,"), '-', classifier) as config, avg(", metric,
                       ") as mean2 FROM resultsView WHERE configurationName LIKE '",set2, "-%' GROUP BY configurationName, classifier) q2 ",
                       "WHERE q1.config=q2.config", sep = "")
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  
  cat(paste("metric:", metric, "\n"))
  cat(paste("mean", set1, ":", mean(results$mean1), "\n"))
  cat(paste("mean", set2, ":", mean(results$mean2)), "\n")
  cat(paste("mean(abs(mean1-mean2)):", mean(abs(results$mean1-results$mean2)), "\n"))
  cat(paste("sd(abs(mean1-mean2)):", sd(abs(results$mean1-results$mean2)), "\n"))
  print(wilcox.test(results$mean1, results$mean2))
  dbDisconnect(mydb)
}

evalRQ2 = function() {
  mydb = dbConnect(MySQL(), user=MYSQL_USER, password=MYSQL_PASSWORT, host=MYSQL_HOST, port=MYSQL_PORT, dbname=MYSQL_DBNAME)
  
  sqlStatement = "SELECT count1, count2, count1/count2 FROM (SELECT count(*) as count1 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' AND recall>=0.75 AND results.precision>=0.75 AND error<=0.25) q1, (SELECT count(*) as count2 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%') q2;"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat(paste("results with baselines for Zimmermann\n"))
  cat(paste("number of single results that fulfill the criterion:       ", results[1,1], "\n"))
  cat(paste("number of single results that do not fulfill the criterion:", results[1,2], "\n"))
  cat(paste("rate of fulfilling the criterion:                          ", results[1,3], "\n"))
  cat("\n")
  
  sqlStatement = "SELECT count1, count2, count1/count2 FROM (SELECT count(*) as count1 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' and configurationName NOT LIKE '%CV' AND classifier!='FIX' AND classifier!='RANDOM' AND recall>=0.75 AND results.precision>=0.75 AND error<=0.25) q1, (SELECT count(*) as count2 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' and configurationName NOT LIKE '%CV' AND classifier!='FIX' AND classifier!='RANDOM') q2;"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat(paste("results without baselines for Zimmermann\n"))
  cat(paste("number of single results that fulfill the criterion:       ", results[1,1], "\n"))
  cat(paste("number of single results that do not fulfill the criterion:", results[1,2], "\n"))
  cat(paste("rate of fulfilling the criterion:                          ", results[1,3], "\n"))
  cat("\n")
  
  sqlStatement = "SELECT count1, count2, count1/count2 FROM (SELECT count(*) as count1 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' AND recall>=0.70 AND results.precision>=0.5) q1, (SELECT count(*) as count2 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%') q2;"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat(paste("results with baselines for He\n"))
  cat(paste("number of single results that fulfill the criterion:       ", results[1,1], "\n"))
  cat(paste("number of single results that do not fulfill the criterion:", results[1,2], "\n"))
  cat(paste("rate of fulfilling the criterion:                          ", results[1,3], "\n"))
  cat("\n")
  
  sqlStatement = "SELECT count1, count2, count1/count2 FROM (SELECT count(*) as count1 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' and configurationName NOT LIKE '%CV' AND classifier!='FIX' AND classifier!='RANDOM' AND recall>=0.7 AND results.precision>=0.5) q1, (SELECT count(*) as count2 FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' and configurationName NOT LIKE '%CV' AND classifier!='FIX' AND classifier!='RANDOM') q2;"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat(paste("results without baselines for He\n"))
  cat(paste("number of single results that fulfill the criterion:       ", results[1,1], "\n"))
  cat(paste("number of single results that do not fulfill the criterion:", results[1,2], "\n"))
  cat(paste("rate of fulfilling the criterion:                          ", results[1,3], "\n"))
  cat("\n")
  
  sqlStatement = "SELECT res1.configurationName, res1.classifier, cnt1 as cnt, cnt1/cnt2 as rate FROM (SELECT configurationName, classifier, count(*) as cnt1 FROM results WHERE configurationName  NOT LIKE 'F%' AND recall>=0.75 AND results.precision>=0.75 AND error<=0.25 GROUP BY configurationName, classifier) as res1, (SELECT configurationName, classifier, count(*) as cnt2 FROM results WHERE configurationName  NOT LIKE 'F%' GROUP BY configurationName, classifier) as res2 WHERE res1.configurationName=res2.configurationName AND res1.classifier=res2.classifier ORDER BY rate DESC;"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat("distinct products with at least 100 classes that fulfill the criterion:\n")
  cat("configuration\t\t\tclassifier\t\tcount\trate\n")
  for( i in 1:nrow(results) ) {
    cat(paste(str_pad(results[i,1], width=25, side="right"), str_pad(results[i,2], width=20, side="right"), results[i,3], results[i,4], "\n", sep="\t"))
  }
  cat("\n")
  
  sqlStatement = "SELECT distinct productName FROM results WHERE configurationName  NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' AND recall>=0.75 AND results.precision>=0.75 AND error<=0.25 AND testsize>=100 GROUP BY configurationName, classifier"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat("distinct products with at least 100 classes that fulfill the criterion:\n")
  cat(paste(results[,1], collapse = "\n"))
  cat("\n\n")
  
  sqlStatement = "SELECT distinct productName FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' AND recall>=0.75 AND results.precision>=0.75 AND error<=0.25 AND testsize<100 GROUP BY configurationName, classifier"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat("distinct products with less than 100 classes that fulfill the criterion:\n")
  cat(paste(results[,1], collapse = "\n"))
  cat("\n\n")
  
  sqlStatement = "SELECT productName as Product, testsize, sum(cnt1) as cnt, prec FROM (SELECT productName, testsize, if(count(*)=1, 1, count(*)/10) as cnt1, 1.0 as recall, (tp+fn)/testsize as 'prec' FROM results WHERE configurationName NOT LIKE 'F%' AND configurationName NOT LIKE 'S%' AND configurationName NOT LIKE '%CV' AND classifier!='FIX' AND classifier!='RANDOM' AND recall>=0.75 AND results.precision>=0.75 AND error<=0.25 GROUP BY configurationName, classifier, productName) as res  GROUP BY productName, testsize UNION SELECT DISTINCT productName as Product, testsize, 0 as cnt, (tp+fn)/testsize as prec FROM results WHERE classifier='FIX' and (tp+fn)/testsize>=0.75 and productName!='pbeans1.csv' ORDER BY Product"
  rs = dbSendQuery(mydb, sqlStatement)
  results = fetch(rs, n=-1)
  cat("distinct products with less than 100 classes that fulfill the criterion:\n")
  cat(paste(results[,1], collapse = "\n"))
  cat("\n\n")
  
  rownames(results) = results$Product
  results$Product = NULL
  results$Data = 'JURECZKO'
  results = results[,c(4,1,2,3)]
  
  captionStr = "Products where any \\ac{CPDP} approach fulfills the criterion by Zimmermann\\etal~\\cite{Zimmermann2009} and the \\emph{precision} of the trivial prediction FIX."
  tableStr = print.xtable(xtable(results, caption=captionStr))
  #tableStr = sub("\\begin{table}", "\\begin{table*}", tableStr, fixed=TRUE)
  tableStr = sub("rlrrr", "|r|l|l|l|l|", tableStr, fixed=TRUE)
  tableStr = sub("\\hline\n", "\\hline\n \\textbf{Product}", tableStr, fixed=TRUE)
  tableStr = sub("Data", "\\textbf{Data Set}", tableStr, fixed=TRUE)
  tableStr = sub("testsize", "\\textbf{\\#Inst.}", tableStr, fixed=TRUE)
  tableStr = sub("cnt", "\\textbf{\\#Appr.}", tableStr, fixed=TRUE)
  tableStr = sub("prec \\\\ \n", "\\emph{precision} \\\\ &&&& FIX \\\\ \\hline\n", tableStr, fixed=TRUE)
  tableStr = gsub(".csv", "", tableStr, fixed=TRUE)
  tableStr = sub("openintents & JURECZKO", "openintents & RELINK", tableStr, fixed=TRUE)
  tableStr = sub("\\end{table}", "\\label{tbl:productsSuccess}\n\\end{table}", tableStr, fixed=TRUE)
  write(tableStr, paste(TABLES_PATH,"productsSuccess.tex", sep=""))
  
  dbDisconnect(mydb)
}

getNemenyiCD <- function(pval, nConfigs, nData) {
  return(qtukey(1 - pval, nConfigs, Inf) / sqrt((nConfigs*(nConfigs+1))/(12*nData)))
}

plotCDNemenyiFriedman <- function(meanRanks, nData, pval=0.05, title="" ) {
  nConfigs = nrow(meanRanks)
  cdNemenyi = getNemenyiCD(pval, nConfigs, nData)
  
  roundedMaxMeanRank = ceiling(max(meanRanks$meanRank))
  roundedMinMeanRank = floor(min(meanRanks$meanRank))
  maxMeanRank = max(meanRanks$meanRank)
  minMeanRank = min(meanRanks$meanRank)
  halfEntries = ceiling(nrow(meanRanks)/2)
  
  red = 0.8
  green = 0
  stepsize = 1.6/max(meanRanks$rank)
  colors = list()
  for( j in 1:max(meanRanks$rank) ) {
    colors[j] = rgb(green, red, 0)
    if( green+stepsize>0.8 ) {
      green = 0.8
      red = red-stepsize
      if( red<0 ) {
        red = 0
      }
    } else {
      green = green+stepsize
      if( green>0.8 ) {
        green = 0.8
      }
    }
  }
  meanRanks <- meanRanks[order(meanRanks$meanRank, decreasing = FALSE),]
  g <- ggplot(data=meanRanks)
  # create lines for configurations
  for( j in 1:nrow(meanRanks)) {
    xval = meanRanks$meanRank[j]
    if( j<=halfEntries ) {
      yend = j
      xlinestart = roundedMinMeanRank
      xtext = roundedMinMeanRank-0.2/roundedMaxMeanRank
      hjust = "right"
    } else {
      yend = j-2*(j-halfEntries)
      xlinestart= roundedMaxMeanRank
      xtext = roundedMaxMeanRank+0.2/roundedMaxMeanRank
      hjust="left"
    }
    if( meanRanks$meanRank[j] <= minMeanRank+cdNemenyi ) {
      color = "#00cc00"
    } 
    else if( meanRanks$meanRank[j] >= maxMeanRank-cdNemenyi ) {
      color = "#cc0000"
    }
    else {
      color = "#0000cc"
    }
    color = colors[[meanRanks$rank[j]]]
    g <- g + geom_segment(x=xval, y=0, xend=xval, yend=-yend, color = color)
    g <- g + geom_segment(x=xlinestart, y=-yend, xend=xval, yend=-yend, color = color)
    g <- g + geom_text(label=meanRanks$config[j], x=xtext, y=-yend, hjust=hjust, color=color)
  }
  # create axis
  if( roundedMaxMeanRank>100 ) {
    breaks = 1:(roundedMaxMeanRank/5)*5
  } else if( roundedMaxMeanRank>80 ) {
    breaks = 1:(roundedMaxMeanRank/4)*4
  } else if( roundedMaxMeanRank>60 ) {
    breaks = 1:(roundedMaxMeanRank/3)*3
  } else if( roundedMaxMeanRank>40 ) {
    breaks = 1:(roundedMaxMeanRank/2)*2
  } else {
    breaks = 1:roundedMaxMeanRank
  }
  breaks <<- breaks
  g <- g + geom_segment(x=0,y=0,xend=roundedMaxMeanRank,yend=0)
  g <- g + geom_text(label="Mean Rank", x=roundedMinMeanRank-0.2/roundedMaxMeanRank, y=0.01*halfEntries, hjust="right", vjust="bottom")
  for( j in roundedMinMeanRank:roundedMaxMeanRank ) {
    if( j %in% breaks ) {
      g <- g + geom_segment(x=j,y=0,xend=j,yend=0.01*halfEntries)
      g <- g + geom_text(label=j, x=j, y=0.015*halfEntries, vjust="bottom")
    }
  }
  # add critical distance
  critDistY = 3
  g <- g + geom_segment(x=minMeanRank,y=critDistY,xend=minMeanRank+cdNemenyi, yend=critDistY)
  g <- g + geom_segment(x=minMeanRank+cdNemenyi,y=critDistY-0.005*halfEntries, xend=minMeanRank+cdNemenyi, yend=critDistY+0.005*halfEntries)
  g <- g + geom_segment(x=minMeanRank,y=critDistY-0.005*halfEntries, xend=minMeanRank, yend=critDistY+0.005*halfEntries)
  g <- g + geom_text(label=paste("Critical Distance =", round(cdNemenyi, digits=3)), x=minMeanRank-0.2/roundedMaxMeanRank, y=critDistY, hjust="right")
  # set scales
  g <- g + scale_y_continuous(limits=c(-halfEntries,1), breaks=NULL, labels=NULL)
  g <- g + scale_x_continuous(limits=c(roundedMinMeanRank-1,roundedMaxMeanRank+1), breaks=breaks, labels=NULL)
  g <- g + theme(axis.ticks.x = element_blank())
  # set title
  g <- g + ggtitle(title)
  g <- g + coord_cartesian(xlim=c(roundedMinMeanRank-0.2*roundedMaxMeanRank,roundedMaxMeanRank+0.2*roundedMaxMeanRank))
  return(g)
}

createMeanRankMat <- function(results, column) {
  nemenyi.groups = factor(as.character(results$config))
  nemenyi.blocks = factor(results$index)
  nemenyi.y = results[[column]]
  nemenyi.n = length(levels(nemenyi.blocks))
  nemenyi.k = length(levels(nemenyi.groups))
  nemenyi.y = nemenyi.y[order(nemenyi.groups, nemenyi.blocks)]
  nemenyi.mat = matrix(nemenyi.y, nrow = nemenyi.n, ncol = nemenyi.k, byrow = FALSE)
  for (nemenyi.i in 1:length(nemenyi.mat[, 1])) nemenyi.mat[nemenyi.i, ] <- rank(nemenyi.mat[nemenyi.i, ])
  nemenyi.mnsum = data.frame(meanRank=colMeans(nemenyi.mat))
  nemenyi.mnsum$config = levels(nemenyi.groups)
  # create ranks
  # break if difference > cd
  nemenyi.cd = getNemenyiCD(0.05, length(unique(results$config)), max(results$index))
  nemenyi.mnsum = nemenyi.mnsum[order(-nemenyi.mnsum$meanRank),]
  currentRank = 1
  nemenyi.mnsum$rank[1] = 1
  for( i in 2:nrow(nemenyi.mnsum ))  {
    if(nemenyi.mnsum$meanRank[i-1]-nemenyi.mnsum$meanRank[i]>nemenyi.cd) {
      currentRank = currentRank+1
    }
    nemenyi.mnsum$rank[i] = currentRank
  }
  nemenyi.mnsum$normRank = 1-(nemenyi.mnsum$rank-1)/(max(nemenyi.mnsum$rank)-1)
  return(nemenyi.mnsum)
}

#####################
# Result evaluation #
#####################

# Performance metrics used
# Script can only handle auc, fscore, gscore, and mcc
# To use other metrics available in the MySQL database, 
# the functions above must be modified, because we do
# some string matching with the metric names.
metricNamesRQ1 = c("auc","fscore","gscore","mcc")

# Datasets that are used for different research questions. 
datasetsRQ1 = c("JURECZKO", "MDP", "AEEEM_LDHHWCHU", "RELINK", "NETGENE")
datasetsRQ3 = c("FILTERJURECZKO")
datasetsRQ4 = c("SELECTEDJURECZKO")

rq1results = evaluateCPDPBenchmark(metricNamesRQ1, datasetsRQ1)
rq1best = plotBestResults(rq1results, "RQ1", "AUC, F-measure, G-Measure, and MCC")
rq3results = evaluateCPDPBenchmark(metricNamesRQ1, datasetsRQ3)
rq4results = evaluateCPDPBenchmark(metricNamesRQ1, datasetsRQ4)
writeResultsTableRQ1(rq1best, rq3results, rq4results)
evalRQ2()

cat("comparing JURECZKO and FILTERJURECZKO\n")
compareDatasets("JURECZKO", "FILTERJURECZKO", "auc")
compareDatasets("JURECZKO", "FILTERJURECZKO", "fscore")
compareDatasets("JURECZKO", "FILTERJURECZKO", "gscore")
compareDatasets("JURECZKO", "FILTERJURECZKO", "mcc")

cat("comparing JURECZKO and SELECTEDJURECZKO\n")
compareDatasets("JURECZKO", "SELECTEDJURECZKO", "auc")
compareDatasets("JURECZKO", "SELECTEDJURECZKO", "fscore")
compareDatasets("JURECZKO", "SELECTEDJURECZKO", "gscore")
compareDatasets("JURECZKO", "SELECTEDJURECZKO", "mcc")

cat("comparing AEEEM and AEEEM_LDHH")
compareDatasets("AEEEM", "AEEEM_LDHH", "auc")
compareDatasets("AEEEM", "AEEEM_LDHH", "fscore")
compareDatasets("AEEEM", "AEEEM_LDHH", "gscore")
compareDatasets("AEEEM", "AEEEM_LDHH", "mcc")

cat("comparing AEEEM and AEEEM_WCHU")
compareDatasets("AEEEM", "AEEEM_WCHU", "auc")
compareDatasets("AEEEM", "AEEEM_WCHU", "fscore")
compareDatasets("AEEEM", "AEEEM_WCHU", "gscore")
compareDatasets("AEEEM", "AEEEM_WCHU", "mcc")

cat("comparing AEEEM and AEEEM_LDHHWCHU")
compareDatasets("AEEEM", "AEEEM_LDHHWCHU", "auc")
compareDatasets("AEEEM", "AEEEM_LDHHWCHU", "fscore")
compareDatasets("AEEEM", "AEEEM_LDHHWCHU", "gscore")
compareDatasets("AEEEM", "AEEEM_LDHHWCHU", "mcc")
