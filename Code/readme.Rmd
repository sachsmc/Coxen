---
title: "readme"
author: "Michael Sachs"
date: "Wednesday, April 09, 2014"
output: html_document
---

### Background

This is an attempt to statistically evaluate the Co-expression Extrapolation (COXEN) algorithm. The method is described [in this paper](http://www.pnas.org/content/104/32/13086). The COXEN algorithm has been applied in several studies to develop a gene expression summary score that can be used to predict treatment success. The key feature of the approach is that the gene expression score is developed *in vitro* by extrapolating gene expression profiles of drug sensitivity on the NCI-60 cell lines to a new type of cancer cell. The success of the algorithm is enthusiastically portrayed by the authors. The goal of this study is to thoroughly evaluate the statistical operating characteristics of the approach.

The COXEN approach relies on the NCI-60 cell line data. The NCI-60 is a set of human cancer cells on which thousands of drug compounds have been tested. For each compound, the drug and dose dependent sensitivity is measured experimentally on each type of cancer cell. Further, each cell type is molecularly profiled using a broad array of assays. Protein expression, gene expression, DNA methylation, microRNA, metabolomics, and more are all measured on each cell type. We assume that the new cell panel has the same expression measurements, but does not have the drug sensitivity measurements. 

### Description of COXEN

Considering that background, the steps in the COXEN approach for evaluating a compound $D$ can be described as follows:

  1. Select two equal subsets of cell types from the NCI-60 that exhibit high and low activity to drug $D$. The authors select the 12 most sensitive and 12 most resistent cell types, as measured by GI50, which is the drug concentration that causes 50% cell growth inhibition as compared to control. 
  2. Use the sigificance analysis of microarrays ([Tusher et al. 2001](http://www.pnas.org/content/98/9/5116.full)) with FDR 0.1 to identify a subset of $k$ probe sets measured on the subset of the NCI-60 that are differentially expressed between the sensitive and resistent cell types. 
  3. (This part is hairy in the methods description). Using the subset of probes identified in the previous step, the goal is to identify a set of probes in the new cell panel (call this panel BLA-40). The $k \times k$ correlation matrix between the subset of probes is computed for each of the 2 cell panels, denote them $U$ and $V$. Then the coexpression coefficient for probe $j$ is defined as 
  \[
  rc(j) = \frac{\sum_{k=1}^n(U_{kj} - \overline{U}_j)(V_{kj} - \overline{V}_j)}{\sqrt{\sum_{k=1}^n(U_{kj} - \overline{U}_j)^2 \sum_{k=1}^n(V_{kj} - \overline{V}_j)^2}}. 
  \]
  Note that this is simply the correlation between the columns of $U$ and $V$.
  4. A subset of probes is selected on the basis of the estimated $rc(j)$ values. The authors suggest taking the top 2 percent. 
  5. Using the subset of probes selected at step 4, a model is estimated to predict the activity using the NCI-60 data as the training set. The risk score or linear predictor resulting from this model is the COXEN score. 
  6. The COXEN score is evaluated after experimentally measuring drug activity in the BLA-40 cell line. 
  
### Aims

The stated advantage of the COXEN is that a risk score can be developed in the pre-clinical setting, and then immediately validated on a test set using specimens from a clinical trial. It obviates the need to split human subject data into training and test sets. 

Our goals are to statistically evaluate the COXEN algorithm as a procedure for developing a risk score. Specifically,

  - In terms of feature selection, how frequently does the COXEN select the correct probes? 
  - Similarly, how frequently does the COXEN fail to select correct probes? 
  - How frequently does COXEN select the incorrect probes?
  - Is there any sort of overfitting going on?
  - Can we do better with a more straightforward approach to feature selection?
  
### A general approach

We can generalize and more clearly describe the Coxen algorithm in fewer steps. Therefore we can examine alternatives and identify similar methods for comparison. 

Given a set of $K$ features measured on $n_1$ supervised observations $X$ and $n_2$ unsupervised observations $V$, 

  1. Estimate the univariate associations $\beta_k$ for $k = 1, \ldots, K$ between the features and the outcome among the $n_1$ subjects. 
  2. Select the top $J$ features, according to $\beta_k$. 
  3. Using the $J$ features, identifiy a further subset $J_2$ according to some metric $f_j(X_{n_1 \times J}, V_{n_2 \times J})$. 
  4. Using the $J_2$ features identified in the previous step, estimate a risk score in the supervised data set using whatever model you choose. 
  
Clearly, there are several decisions to make in these steps. First, one must decide on a measure of the association $\beta_k$. The Coxen people use the significance analysis of microarrays approach on a curated subset of the supervised data. Then, the metric $f_j(\cdot)$ need to be defined. Coxen uses the correlation of the correlation matrix of the features. Then the number of features selected needs to be specified. To optimize, one would typically select features by cross-validation. The Coxen folks set arbitrary thresholds. In general, this is a reasonable approach, but clearly there is room for improvement by optimizing these arbitrary decisions and retaining all available data. 
  
  
  
  
  
