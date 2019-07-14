"""A module containing functions for detecting edges representing conditional
dependence, between nodes in a multivariate time series.

References
    [1] https://ieeexplore.ieee.org/document/7084657/
    [2] http://www.math.tau.ac.il/~ybenja/MyPapers/benjamini_yekutieli_ANNSTAT2001.pdf
    [3] https://academic.oup.com/biomet/article-abstract/93/2/399/220995?redirectedFrom=fulltext
    [4] Eichler test statistic paper.
    """

def MultHypTest(data, alpha, control, stat):
    """Uses the method from [1].

    Inputs:
        data: A dictionary object containing the following fields
                  raw: A (m, n) numpy array containing the raw time series data
                       where m is the spatial dimension size and n is the
                       temporal dimension size.
                  pdgm: A (m, m, n) numpy array containing an estimator of the
                        spectral matrix for the data X. This is different
                        depending on the `stat` used.
              Note that the `raw` entry can be empty if the `pdgm` entry is
              already filled. If `pdgm` is empty than we calculate it depending
              on the `stat` option, from the `raw` data.
        alpha: The probabilitiy used in controlling the type I errors. It has
               different meanings for different control types (mentioned
               below).
        control: This is the control method for how to accept/reject tests
                 to control type I errors. It can be either
                     - FDR (False discovery rate). This is defined as
                           E[V/R]
                       where V is total false positives (i.e. reject but should
                       have accepted) and R is total rejections and E[]
                       represents expectation.

                       In this case we use the Benjamini-Hochberg procedure to
                       control
                           E[V/R] <= (m_0/m) * alpha <= alpha
                       where m_0 are the true null hypotheses (unknown) and
                       m is the total number of tests.

                       Note the above inequality holds if the m tests are 
                       independent and some dependent cases [2].

                     - FWER (Familywise error rate). This is defined as
                           Pr(V >= 1)
                       where V is as above i.e. the probabilitiy of making at
                       least one incorrect rejection.

                       In this case we use the Holm-Bonferroni method as in [1]
                       and ensures that
                           Pr(V >= 1) <= alpha.

                        Note that this doesn't require independence of the
                        tests.

                     - FDR_dep (False discovery rate with dependence). This
                       method controls the FDR 
                           E[V/R] <= alpha
                       when the tests are dependent. However it will be much
                       more conservative than simply the 'FDR' method i.e. make
                       more type II errors.
        stat: The multiple hypothesis testing framework allows different
              statistics to be used as long as they are (asymptotically) normal
              under the null hypothesis.
              We allow the use of 2 such statistics:
                  - Matsuda: The test statistic derived in [3].
                  - Eichler: The test statistic derived in [4].
                  
