# Pattern-Recognition-on-Biomedical-Data

Master Thesis - Abstract

Motivation: The large amount of genomic data makes necessary the use
of machine learning and data analysis techniques in order to extract useful
information and discover hidden patterns. In this project, we apply biclustering for biomarker detection on real asthma datasets.
Data: The datasets that we work with are 4. The first consists of clinical
information for subjects that can be organized in 5 categories: general, blood,
lung function, sputum, biopsy. The second and third contain gene expression
and DNA methylation data, respectively. The last one includes micro RNA
expression data.
Methods: We use a biclustering method, called ’FABIA: Factor Analysis for
Bicluster Acquisition’. FABIA is a generative multiplicative model which
is based on linear dependencies between gene expression and conditions
to form biclusters. t-SNE algorithm is also used for the visualization of
the data in 2 dimensions. We run FABIA algorithm multiple times and we
rank the resulting biclusters. For the evaluation of the biclusters, 4 quality
measures are being used. The first, information content, shows the amount
of information each bicluster contains about the data. The next ones are the
variance and mean squared residue (MSR). The fourth quality measure is
the virtual error, which shows the tendency that genes follow under a set of
conditions.
Results: After applying FABIA multiple times on the datasets and getting
the biclusters, we choose to examine only the robust ones. We consider a
bicluster as robust if it has average overlap percentage more than 80% over
the runs. These biclusters gave us combinations of particular features from
the first dataset that may indicate the existence of asthma.
