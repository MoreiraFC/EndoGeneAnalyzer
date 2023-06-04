# EndoGeneAnalyzer
EndoGeneAnalyzer is a dynamic R shiny tool that simplifies and assists in selecting reference genes in scientific studies and performing differential gene expression analysis for RT-qPCR data. 



The tool has been developed to be intuitive, interactive, and efficient, with several steps to guide the user. In the first step, the user enters the data with the option to choose between the supported file formats: .xls/.xlsx or .txt/.csv. This flexibility makes it easy to import data from different sources. After loading the data, the user selects the targets of interest for analysis (non-reference genes). The next step focuses on evaluating the best reference set of genes based on mean variation using descriptive statistical data such as gene standard deviation, the sum of square differences between the mean of each group and the gene mean, and the sum of squares differences between the standard deviation of each group and the gene standard deviation. It also calculates stability metrics using NormFinder, which helps to identify the genes that best fit the study conditions. One of the critical innovations of EndoGeneAnalyzer is the ability to analyze the stability of reference genes. In this step, the tool allows the user to identify and remove outliers, which are samples with ∆Ct mean values above or below a user-defined threshold (default = 2 standard deviations), providing flexibility in the analysis.
Finally, the EndoGeneAnalyzer can perform differential expression analysis using the target ∆Ct and the mean ∆Ct of the set of reference genes. This step allows accurate and efficient comparisons between different groups or conditions, further, delivering a fold change result

Data Upload

This step is critical to ensure correct information is used during the analysis. The input file must contain the following columns: i) the first column with the sample names; ii) the following columns with the mean Ct values of specific targets and reference genes for the sample; and iii) the last column with information about the groups or conditions to which the samples belong.
The data can be imported in two ways: i) Excel tables (.xls / .xlsx), in this option, the tool does not require modification of the decimal separator; and ii) text tables (.txt / .csv), in this option the default decimal separator is dot(.) and it is necessary to configure the text delimiter.
Finally, after verifying the correct formatting of the table, the user needs to click on the "Confirm Data Table" button to proceed with the analysis process. This final step ensures that the tool correctly recognizes and validates the data provided.

Data Summary

The selection of target genes (non-reference) is essential as it is crucial to guide the analysis,  these genes are related to the research objectives. Once the target-gene(s) have been identified, the user must click the "Update Target Gene" button to confirm the selection. This step ensures that the selected genes are processed and included in the analysis.

Reference Genes and Sample Stability: Removing Outliers.

Outliers are atypical data values that can be identified in RT-qPCR data; experimental errors are the leading cause of these occurrences. These errors are related to environmental conditions, instrument calibration problems, or other sources of uncontrolled variation that may occur during the experiment. It is essential to be aware of outliers and to understand the potential impact of their removal on results (18).
EndoGeneAnalyzer identifies outliers per group for each gene and their removal can be easily performed using an available function. By default, the tool considers a sample as an outlier if the mean ∆Ct is greater or less than 2 standard deviations from the mean of the group/condition to which the sample belongs for the reference gene. This value can be configured according to the user's preferences. The tool offers two methods of outlier removal: i) removal of all outliers; ii) removal of those that directly interfere with the mean Ct values of the reference genes.

Reference Genes Analysis

This is a crucial step in the tool's operation, as it provides information about the reference genes and their variation between the different groups or conditions studied. At this stage, significant changes in the reference genes are observed, especially in the mean values between the analyzed groups.
The first table generated is the "Reference of genes per group", which presents information about the variation observed between the groups or conditions studied for each reference gene or the averages of the group of reference genes The statistical tests used are Wilcoxon-Mann-Whitney (2 groups) or Kruskall-Wallis/Dunn (3 or more groups). At this stage, it is expected that there will be no significant changes (p-value < in reference genes between the studied groups or conditions.
The tool also provides the "Gene Reference Descriptive Statistics" table, which presents three fundamental values for assessing the reference genes: gene standard deviation, sum of squared differences between the mean of each group and the gene mean, and sum of squared differences between the standard deviation of each group and the gene standard deviation. The formulas used to calculate the sum of squared differences are as follows: n is the number of groups; µi is the mean Ct of the group; µg is the mean Ct of the gene; σi is the standard deviation of the group and σg  is the standard deviation of the gene.
sum.mean.square.diff=i=0ngroup(µi-µg)2    sum.SD.square.diff=i=0ngroup(σi-σg)2 
In addition, NormFinder provides information on the stability and suitability of reference, since this software is integrated into our tool's interface.

Differential Analysis

EndoGeneAnalyzer allows for comparing gene expression differences among the investigated groups using ∆Ct. ∆Ct is calculated as the difference between the target gene and the mean of the reference genes. For 2 groups, two statistical tests integrated into the tool are available: Pearson t-test and Wilcoxon-Mann-Whitney Rank Sum test. For the comparison between 3 or more groups ANOVA/Tukey and Kruskall-Wallis/Dunn are applied.
The system also calculates the Shapiro test for the normality of each group and the Fold-Change using the formula 2-∆∆CT. This metric quantifies the difference in expression between two groups, considering the relative variation of ∆Ct values.
