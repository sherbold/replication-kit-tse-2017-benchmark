Introduction
============
Within this archive you find the replication package for the article "A Comparative Study to Benchmark Cross-project Defect Prediction Approaches" by Steffen Herbold, Alexander Trautsch, and Jens Grawbowski for currently under review at IEEE Transactions on Software Engineering. The aim of this replication package is to allow other researchers to replicate our results with minimal effort, as well as to provide additional results that could not be included in the article directly. 

Requirements
============
- MySQL
- Java 8
- R 3.3.2
- Optional: JDK and Ant for compilation of the source code.
- Tested on Windows 10, but does not contain OS specific code. 

Contents
========
- additional-visualizations folder with Critical Distance (CD) diagrams and boxplots, that could not be included in the article.
- benchmark-execution folder with all files required to execute the benchmark and replicate the results.
- R-scripts folder with scripts for the evaluation of results and generation of plots.
- sql folder with SQL code for the setup of a local MySQL database for the benchmark results.
- raw-results folder with the raw results of the benchmark stored as a Comma Separated Value (CSV) file.

How to use
==========
There are three ways to use this replication kit.
1. Look at the additional visualizations. For this, simply download the folder and open the boxplots.html/cd-diagrams.html file.
2. Get access to the raw results of the benchmark. For this, you may either use the [raw CSV data](raw-results/crosspare-results-replication-kit.csv) directly, or import the data into a MySQL data base.
3. Replicate the results of the benchmark, including running the evaluation scripts that execute the statistical tests and create plots/tables for the results.

In the following, we will explain:
- how to setup your local MySQL database for using this replication kit;
- how to import the raw benchmark results into the MySQL database;
- how the replicate the raw results; and
- how to run the evaluation scripts. 

Setup local a MySQL database
----------------------------

1. Install MySQL - any version should work, we do not use special features.
2. Create a database called "crosspare".
3. Create a user called "crosspare" with password "crosspare".
4. Create the table results and the view resultsView by executing the SQL code in the [setup-db.sql](sql/setup-db.sql) script.

It is also possible to a different database name, user name, and password. In that case, the [mysql.cred](benchmark-execution/mysql.cred) muts be adopted for the execution of the benchmark and the global variables at the beginning of the R scripts for the evaluation and plot generation. 

The database contains more columns for metrics than are reported on within the paper, among others:
- error rate
- true negative rate
- the confusion matrix itself (tp, fp, tn, fn columns).

**There are several other columns within the as well, however, they should be ignored as they are experimental and where not tested properly**. Most of them result rather from testing what is possible with CrossPare and where hacked in rather quickly without any double checking, e.g., AUCEC to see if the evaluation framework would be able to take effort into account. However, this works only for some data sets, and even there buggy. *Please also note that this comment holds only for the version of CrossPare contained in the replication package. Further development since April 2017 may have fixed, removed or added the calculation of metrics.*

Raw benchmark results
------------------------------------------
The raw data is available in two formats.
- The [db-inserts.sql](raw-results/db-inserts.sql) SQL script with SQL insert statements to populate the results table.
- The [crosspare-results-replication-kit.csv](raw-results/crosspare-results-replication-kit.csv) Comma Separated Values (CSV) file for exploration of the results without the requirement of a MySQL database. 

Replicate raw results
---------------------
It is also possible to execute the benchmark and replicate all raw results from scratch. To execute the replication, you can simply use the batch scripts for Linux and Windows. Please note that they try to create a Java Virtual Machine (JVM) with access to roughly 30 GB of heap space. Sometimes, this is still not enough, in which case the execution will simply crash. However, you can just restart the execution, and the calculation of results will continue where it left of. The only solution to completely remove the crashs would be lots more memory. 

Please note, that executing the full benchmark if the database is already populated with the raw results we provided will not work. Ideally, nothing will happen, because the benchmark execution checks which results are already available, and only executes missing parts. More likely, the execution will crash because the MySQL database will at some point reject the connection. This happens because of the rapid checking if results are already available, which will lead to thousands of queries within seconds, if the database is already populated. 

The replication package uses the tool [CrossPare](https://github.com/sherbold/CrossPare/). CrossPare implements many approaches proposed for cross-project defect prediction and supports loading data from various published data sets. The experiments themselves are defined using XML files, that describe which data is loaded and how the defect prediction models are trained. All experiment configurations are contained in the [config](benchmark-execution/benchmark/config) folder. The [config-infeasible](benchmark-execution/benchmark/config-infeasible) folder contains the configurations that could not be executed due to memory and runtime constraints. 

The crosspare.jar gets the paths to the config folders as parameters and loads all the experiment configurations. Then, these experiment configurations are executed. The experiments themselves nearly all only use a single CPU. However, CrossPare is implemented in such a way that all CPU cores a utilized by executing multiple experiments in parallel. Thus, CPU usage is usually at 100% for all cores at all times, while the benchmark is running. Due to the large number of experiment configurations, this may take some weeks. You can use the exec_all scripts to executed all experiments or modify the exec_single_config scripts to execute only a single experiment configuration. 

The results of the execution are stored in the MySQL database. Additionally, the results-csv folder contains on-the-fly CSV outputs for each experiment configuration. However, these are overriden and may be invalid in case crosspare crashes. The MySQL database, on the other hand, always stays consistent. 

Run evaluation scripts
----------------------
The replication kit also contains our scripts for the evaluation of the results. The scripts contain R code. We recommend using [RStudio](https://www.rstudio.com/) for the execution of the code. The scripts try to install some libraries from [CRAN](http://cran.us.r-project.org/). Make sure that you have a internet connection available and that the CRAN mirror is properly configured in your R environment. 

The [generate_results.R](R-scripts/generate_results.R) performs most of the this. The script contains the logic for the statistical comparison of the results, as well as the generation of all plots and results tables contained in the article. At the beginning of the file, global variables define some options for the execution, e.g., if the plots should be generated, where they should be stored, and the connection details for the MySQL database. 

Additionally, the script [generate_htmlplots.R](R-scripts/generate_htmlplots.R) contains code that generates the [boxplots](additional-visualizations/boxplots.html) which are contained in this replication kit. Same as above, the script also contains some global variable at the beginning of the script for the location and the connection details for the MySQL database. 

Building CrossPare from Source
==============================
The source code for CrossPare provided within the replication kit in the [crosspare-src](crosspare-src) folder. The crosspare.jar file can be built using the provided Ant script by calling "ant dist". The build results will then appear in the folder "dist". The source code is equal to [revision cf4239f689d315c9cf84e4c05287006c673dc911](https://github.com/sherbold/CrossPare/tree/cf4239f689d315c9cf84e4c05287006c673dc911).

License
=======
This replication package and CrossPare are licensed under the Apache License, Version 2.0. 