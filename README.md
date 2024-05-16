# CFDscoRe

This is an implementation of the CFD score from Doench et al 2016. There is no reference implementation of the CFD score calculation, and the available implementations are either incomplete (focusing only on the mismatch penalties) or are known to contain a bug in the calculation of DNA-bulge containing alignments.

# Installation

To install this package clone it, change to the cloned directory, and start an R session. From R run:

```
devtools::install()
```

To do this the R package `devtools` must be installed.

# Usage

For a given guide/target-site alignment, e.g.:

<pre>
 guide: CA-TGCCGTGTGTACCATGAC
          *                 *
target: CAGTGCCATGTGTACCATCAGNGG
</pre>


The CFD score of the alignment can be calculated by calling:

```
CFDscoRe::cfd_score('CA-TGCCGTGTGTACCATGAC', 'CAGTGCCATGTGTACCATCAG', 'GG')
```
