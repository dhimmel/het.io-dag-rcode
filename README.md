# R Code for Hetio Disease-Gene Prediction Study

This repository contains R code for the study:

> **Heterogeneous Network Edge Prediction: A Data Integration Approach to Prioritize Disease-Associated Genes**  
Daniel S. Himmelstein, Sergio E. Baranzini  
*PLOS Computational Biology* (2015-07-09) <https://doi.org/98q>  
DOI: [10.1371/journal.pcbi.1004259](https://doi.org/10.1371/journal.pcbi.1004259) · PMID: [26158728](https://www.ncbi.nlm.nih.gov/pubmed/26158728) · PMCID: [PMC4497619](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4497619)

## Other resources

Datasets related to this study are available at <https://github.com/dhimmel/het.io-dag-data>.
Predictions from this study can be viewed at <https://het.io/disease-genes/>.
Python code related to this study is available at <https://github.com/dhimmel/het.io-dag-pycode>.

## Documentation

This repository is provided for historical and archiving purposing.
There is not much documentation on how to use the code.
We are happy address specific questions via GitHub Issues.

## History

This repository originally used mercurial and was hosted on BitBucket at `https://bitbucket.org/dhimmel/rhetio`.
Unfortunately, BitBucket [deleted](https://community.atlassian.com/t5/Bitbucket-articles/What-to-do-with-your-Mercurial-repos-when-Bitbucket-sunsets/ba-p/1155380) mercurial repos in 2020 without an automated migration path.
On 2020-09-10, dhimmel realized this source code was no longer online.

The repo could not be retrieved from BitBucket, but dhimmel had a local copy.
The local copy matched the version archived by [Software Heritage](https://archive.softwareheritage.org/swh:1:dir:ba7e1b33e1a615b2d8e42b179e82568166e463c4;origin=https://bitbucket.org/dhimmel/rhetio;visit=swh:1:snp:1a4eadb899ff36036b41d3112cff08d7e2d3fabb;anchor=swh:1:rev:314c83bdf3012973d526fd02faa2e9bd3d1b261c;path=//).
Using the [hg-git](https://hg-git.github.io/) mercurial plugin, dhimmel pushed the repo to <https://github.com/dhimmel/het.io-dag-rcode>.

The local repo was dirty, so dhimmel [committed](https://github.com/dhimmel/het.io-dag-rcode/commit/952bd5ccc0aebff3f4ca887babfe4404f15af8bd) and [reverted](https://github.com/dhimmel/het.io-dag-rcode/commit/e1b05928a4e7fbc5ac63abbc551b28d33075752f) those changes.
Going forward the plan is to switch to using git to manage this repository.
