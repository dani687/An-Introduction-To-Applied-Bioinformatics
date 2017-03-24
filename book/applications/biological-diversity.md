```python
>>> # DELETE ME!
... import qiime2.plugins
>>> qiime2.plugins.available_plugins()
{'qiime2.plugins.alignment',
 'qiime2.plugins.composition',
 'qiime2.plugins.dada2',
 'qiime2.plugins.demux',
 'qiime2.plugins.diversity',
 'qiime2.plugins.emperor',
 'qiime2.plugins.feature_classifier',
 'qiime2.plugins.feature_table',
 'qiime2.plugins.phylogeny',
 'qiime2.plugins.taxa',
 'qiime2.plugins.types'}
```

---
deletable: true
editable: true
...


# Studying Microbial Diversity <link src='2bb2cf'/>

---
deletable: true
editable: true
...

Unicellular organisms, often referred to as “microbes”, represent the vast majority of the diversity of life on Earth. Microbes perform an amazing array of biological functions and rarely live or act alone, but rather exist in complex communities composed of many interacting species. We now know that the traditional approach for studying microbial communities (or microbiomes), which relied on microbial culture (in other words, being able to grow the microbes in the lab), is insufficient because we don’t know the conditions required for the growth of most microbes. Recent advances that have linked microbiomes to processes ranging from global (for example, the cycling of biologically essential nutrients such as carbon and nitrogen) to personal (for example, human disease, including obesity and cancer) have thus relied on “culture independent” techniques. Identification now relies on sequencing fragments of microbial genomes, and using those fragments as “molecular fingerprints” that allow researchers to profile which microbes are present in an environment. Currently, the bottleneck in microbiome analysis is not DNA sequencing, but rather interpreting the large quantities of DNA sequence data that are generated: often on the order of tens to hundreds of gigabytes. This chapter will integrate many of the topics we've covered in previous chapters to introduce how we study communities of microorganisms using their DNA sequences.

---
deletable: true
editable: true
...

## Getting started: the feature table <link src='23ef6e'/>

---
deletable: true
editable: true
...

From a bioinformatics perspective, studying biological diversity is centered around a few key pieces of information:

---
deletable: true
editable: true
...

* A table of the frequencies of certain biological features (e.g., species or OTUs) on a per sample basis.
* *Sample metadata* describing exactly what each of the samples is, as well as any relevant technical information.
* *Feature metadata* describing each of the features. This can be taxonomic information, for example, but we'll come back to this when we discuss features in more detail (this will be completed as part of [#105](https://github.com/caporaso-lab/An-Introduction-To-Applied-Bioinformatics/issues/105)).
* Optionally, information on the relationships between the biological features, typically in the form of a phylogenetic tree where tips in the tree correspond to OTUs in the table.

---
deletable: true
editable: true
...

None of these are trivial to generate (defining OTUs was described in the [OTU clustering chapter](../algorithms/5-sequence-mapping-and-clustering.ipynb), building trees in the [Phylogenetic reconstruction chapter](../algorithms/3-phylogeny-reconstruction.ipynb), and there is a lot of active work on standardized ways to describe samples in the form of metadata, for example [Yilmaz et al (2011)](http://www.nature.com/nbt/journal/v29/n5/full/nbt.1823.html) and the [isa-tab](http://isa-tools.org/) project. For this discussion we're going to largely ignore the complexities of generating each of these, so we can focus on how we study diversity once we know what features were observed in each of our samples.

---
deletable: true
editable: true
...

The sample by feature frequency table is central to investigations of biological diversity, and we'll refer to it here as the *feature table*. The Genomics Standards Consortium has recognized the [Biological Observation Matrix](http://www.biom-format.org) ([McDonald et al. (2011) *Gigascience*](http://www.gigasciencejournal.com/content/1/1/7)) file format definition as a community standard for representing those tables, so you may sometimes hear these referred to as *BIOM tables*. Since the features in these tables are very commonly OTUs, you may also hear these referred to as *OTU tables*.

---
deletable: true
editable: true
...

The basic data that goes into a feature table is the list of sample ids, the list of feature ids, and the frequency matrix, which describes how many times each feature was observed in each sample. We can build and display a feature table as follows:

---
deletable: true
editable: true
...

```python
>>> %pylab inline
...
>>> from __future__ import division
>>> import numpy as np
>>> import pandas as pd
...
>>> sample_ids = ['sampleA', 'sampleB', 'sampleC']
>>> feature_ids = ['feature-1', 'feature-2', 'feature-3', 'feature-4', 'feature-5']
>>> data = array([[1, 3, 0, 1, 0],
...               [0, 2, 0, 4, 4],
...               [0, 0, 6, 2, 1]])
...
>>> table1 = pd.DataFrame(data, columns=feature_ids, index=sample_ids)
>>> table1
Populating the interactive namespace from numpy and matplotlib
         feature-1  feature-2  feature-3  feature-4  feature-5
sampleA          1          3          0          1          0
sampleB          0          2          0          4          4
sampleC          0          0          6          2          1
```

---
deletable: true
editable: true
...

If we want the feature frequency vector for sample `A` from the above table, we use the pandas API to get that as follows:

---
deletable: true
editable: true
...

```python
>>> table1.loc['sampleA']
feature-1    1
feature-2    3
feature-3    0
feature-4    1
feature-5    0
Name: sampleA, dtype: int64
```

---
deletable: true
editable: true
...

**TODO**: Trees in Newick format; sample metadata in TSV format, and loaded into a pandas DataFrame.

---
deletable: true
editable: true
...

Before we start looking at what we can do with this data once we have it, let's discuss some terminology.

```python
>>> import qiime2
>>> table1 = qiime2.Artifact.import_data('FeatureTable[Frequency]', table1)
```

```python
>>> table1.view(pd.DataFrame)
         feature-1  feature-2  feature-3  feature-4  feature-5
sampleA        1.0        3.0        0.0        1.0        0.0
sampleB        0.0        2.0        0.0        4.0        4.0
sampleC        0.0        0.0        6.0        2.0        1.0
```

```python
>>> import biom
>>> table1.view(biom.Table)
5 x 3 <class 'biom.table.Table'> with 9 nonzero entries (60% dense)
```

---
deletable: true
editable: true
...

## Terminology <link src='e67fbb'/>

---
deletable: true
editable: true
...

There are literally hundreds of metrics of biological diversity. Here is some terminology that is useful for classifying these metrics.

---
deletable: true
editable: true
...

**Alpha versus beta diversity**

---
deletable: true
editable: true
...

 * $\alpha$ (i.e., within sample) diversity: Who is there? How many are there?
 * $\beta$ (i.e., between sample) diversity: How similar are pairs of samples?

---
deletable: true
editable: true
...

**Quantitative versus qualitative metrics**

---
deletable: true
editable: true
...

 * qualitative metrics only account for whether an organism is present or absent
 * quantitative metrics account for abundance

---
deletable: true
editable: true
...

**Phylogenetic versus non-phylogenetic metrics**

---
deletable: true
editable: true
...

 * non-phylogenetic metrics treat all OTUs as being equally related
 * phylogenetic metrics incorporate evolutionary relationships between the OTUs

---
deletable: true
editable: true
...

In the next sections we'll look at some metrics that cross these different categories. As new metrics are introduced, try to classify each of them into one class for each of the above three categories.

---
deletable: true
editable: true
...

## Measuring alpha diversity <link src='200e93'/>

---
deletable: true
editable: true
...

The first type of metric that we'll look at will be alpha diversity, and we'll specifically focus on *richness* here. Richness refers to how many different *types* of organisms are present in a sample: for example, if we're interested in species richness of plants in the Sonoran Desert and the Costa Rican rainforest, we could go to each, count the number of different species of plants that we observe, and have a basic measure of species richness in each environment.

---
deletable: true
editable: true
...

An alternative type of alpha diversity measure would be *evenness*, and would tell us how even or uneven the distribution of species abundances are in a given environment. If, for example, the most abundant plant in the Sonoran desert was roughly as common as the least abundant plant (not the case!), we would say that the evenness of plant species was high. On the other hand, if the most abundant plant was thousands of times more common than the least common plant (probably closer to the truth), then we'd say that the evenness of plant species was low. We won't discuss evenness more here, but you can find coverage of this topic (as well as many of the others presented here) in [Measuring Biological Diversity](http://www.amazon.com/Measuring-Biological-Diversity-Anne-Magurran/dp/0632056339).

---
deletable: true
editable: true
...

Let's look at two metrics of alpha diversity: observed species, and phylogenetic diversity.

---
deletable: true
editable: true
...

### Observed species (or Observed OTUs) <link src='e3f5c8'/>

---
deletable: true
editable: true
...

Observed species, or Observed OTUs as it's more accurately described, is about as simple of a metric as can be used to quantify alpha diversity. With this metric, we simply count the OTUs that are observed in a given sample. Note that this is a qualitative metric: we treat each OTU as being observed or not observed - we don't care how many times it was observed.

---
deletable: true
editable: true
...

Let's define a new table for this analysis:

```python
>>> data = array([[1, 1, 5],
...               [1, 2, 0],
...               [3, 1, 0],
...               [0, 2, 0],
...               [0, 0, 0],
...               [0, 0, 3],
...               [0, 0, 1]]).T
>>> print(data)
[[1 1 3 0 0 0 0]
 [1 2 1 2 0 0 0]
 [5 0 0 0 0 3 1]]
```

---
deletable: true
editable: true
...

```python
>>> sample_ids = ['A', 'B', 'C']
>>> feature_ids = ['B1','B2','B3','B4','B5','A1','E2']
>>> data = array([[100, 100, 300, 0, 0, 0, 0],
...               [100, 200, 100, 200, 0, 0, 0],
...               [500, 0, 0, 0, 0, 300, 100]])
...
>>> table2 = pd.DataFrame(data, columns=feature_ids, index=sample_ids)
>>> table2 = qiime2.Artifact.import_data('FeatureTable[Frequency]', table2)
>>> table2.view(pd.DataFrame)
      B1     B2     B3     B4   B5     A1     E2
A  100.0  100.0  300.0    0.0  0.0    0.0    0.0
B  100.0  200.0  100.0  200.0  0.0    0.0    0.0
C  500.0    0.0    0.0    0.0  0.0  300.0  100.0
```

---
deletable: true
editable: true
...

Our sample $A$ has an observed OTU frequency value of 3, sample $B$ has an observed OTU frequency of 4, and sample $C$ has an observed OTU frequency of 3. Note that this is different than the total counts for each column (which would be 5, 6, and 9 respectively). Based on the observed OTUs metric, we could consider samples $A$ and $C$ to have even OTU richness, and sample $B$ to have 33% higher OTU richness.

---
deletable: true
editable: true
...

We could compute this in python as follows:

---
deletable: true
editable: true
...

```python
>>> def observed_otus(table, sample_id):
...     return sum([e > 0 for e in table.loc[sample_id]])
```

---
deletable: true
editable: true
...

```python
>>> print(observed_otus(table2.view(pd.DataFrame), 'A'))
3
```

---
deletable: true
editable: true
...

```python
>>> print(observed_otus(table2.view(pd.DataFrame), 'B'))
4
```

---
deletable: true
editable: true
...

```python
>>> print(observed_otus(table2.view(pd.DataFrame), 'C'))
3
```

---
deletable: true
editable: true
...

### Even sampling <link src='2466b2'/>

---
deletable: true
editable: true
...

Imagine again that we're going out to count plants in the Sonoran Desert and the Costa Rican rainforest. We're interested in getting an idea of the plant richness in each environment. In the Sonoran Desert, we survey a square kilometer area, and count 150 species of plants. In the rainforest, we survey a square meter, and count 15 species of plants. So, clearly the plant species richness in the Sonoran Desert is higher, right? What's wrong with this comparison?

---
deletable: true
editable: true
...

The problem is that we've expended a lot more sampling effort in the desert than we did in the rainforest, so it shouldn't be surprising that we observed more species there. If we expended the same effort in the rainforest, we'd probably observe a lot more than 15 or 150 plant species, and we'd have a more sound comparison.

---
deletable: true
editable: true
...

In sequencing-based studies of microorganism richness, the analog of sampling area is sequencing depth. If we collect 100 sequences from one sample, and 10,000 sequences from another sample, we can't directly compare the number of observed OTUs or the phylogenetic diversity of these because we expended a lot more sampling effort on the sample with 10,000 sequences than on the sample with 100 sequences. The way this is typically handled is by randomly subsampling sequences from the sample with more sequences until the sequencing depth is equal to that in the sample with fewer sequences. If we randomly select 100 sequences at random from the sample with 10,000 sequences, and compute the alpha diversity based on that random subsample, we'll have a better idea of the relative alpha diversities of the two samples.

---
deletable: true
editable: true
...

```python
>>> sample_ids = ['A', 'B', 'C']
>>> feature_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
>>> data = array([[50, 35, 100, 15, 0],
...               [4, 200, 2, 400, 40],
...               [0, 0, 1, 1, 1]])
...
>>> bad_table = pd.DataFrame(data, columns=feature_ids, index=sample_ids)
>>> bad_table
   OTU1  OTU2  OTU3  OTU4  OTU5
A    50    35   100    15     0
B     4   200     2   400    40
C     0     0     1     1     1
```

---
deletable: true
editable: true
...

```python
>>> print(observed_otus(bad_table, 'A'))
4
```

---
deletable: true
editable: true
...

```python
>>> print(observed_otus(bad_table, 'B'))
5
```

---
deletable: true
editable: true
...

```python
>>> print(observed_otus(bad_table, 'C'))
3
```

---
deletable: true
editable: true
...

```python
>>> print(bad_table.sum(axis=1))
A    200
B    646
C      3
dtype: int64
```

```python
>>> bad_table = qiime2.Artifact.import_data("FeatureTable[Frequency]", bad_table)
```

```python
>>> import qiime2.plugins.feature_table
>>> rarefied_result = qiime2.plugins.feature_table.actions.rarefy(bad_table, 200)
>>> rarefied_table = rarefied_result.rarefied_table
```

```python
>>> rarefied_table.view(pd.DataFrame)
   OTU1  OTU2   OTU3   OTU4  OTU5
A  50.0  35.0  100.0   15.0   0.0
B   1.0  60.0    1.0  127.0  11.0
```

```python
>>> rarefied_table.view(pd.DataFrame).sum(axis=1)
A    200.0
B    200.0
dtype: float64
```

---
deletable: true
editable: true
...

**TODO**: Add alpha rarefaction discussion.

```python
>>> table2.view(pd.DataFrame).sum(axis=1)
A    500.0
B    600.0
C    900.0
dtype: float64
```

```python
>>> rarefy_result = qiime2.plugins.feature_table.actions.rarefy(table2, 500)
>>> table2_rarefied = rarefy_result.rarefied_table
>>> table2_rarefied.view(pd.DataFrame)
      B1     B2     B3     B4     A1    E2
A  100.0  100.0  300.0    0.0    0.0   0.0
B   78.0  160.0   90.0  172.0    0.0   0.0
C  285.0    0.0    0.0    0.0  169.0  46.0
```

```python
>>> table2_rarefied.view(pd.DataFrame).sum(axis=1)
A    500.0
B    500.0
C    500.0
dtype: float64
```

```python
>>> import qiime2.plugins.diversity
...
>>> alpha_result = qiime2.plugins.diversity.actions.alpha(table=table2_rarefied, metric='observed_otus')
>>> alpha_result.alpha_diversity.view(pd.Series)
A    3
B    4
C    3
Name: observed_otus, dtype: int64
```

QIIME2 can save you from performing analyses incorrectly... 
```

$ qiime2.plugins.diversity.actions.alpha(table2, 'observed_otus')


---------------------------------------------------------------------------
TypeError                                 Traceback (most recent call last)
<ipython-input-74-d63413ad9748> in <module>()
      1 import qiime2.plugins.diversity
      2 
----> 3 qiime2.plugins.diversity.actions.alpha(table2, 'observed_otus')

<decorator-gen-130> in alpha(table, metric)

/Users/gregcaporaso/miniconda3/envs/iab-latest/lib/python3.5/site-packages/qiime2-2017.2.0-py3.5.egg/qiime2/sdk/action.py in callable_wrapper(*args, **kwargs)
    149             user_input.update(kwargs)
    150 
--> 151             self.signature.check_types(**user_input)
    152             output_types = self.signature.solve_output(**user_input)
    153 

/Users/gregcaporaso/miniconda3/envs/iab-latest/lib/python3.5/site-packages/qiime2-2017.2.0-py3.5.egg/qiime2/core/type/signature.py in check_types(self, **kwargs)
    280             if kwargs[name] not in spec.qiime_type:
    281                 raise TypeError("Argument to input %r is not a subtype of"
--> 282                                 " %r." % (name, spec.qiime_type))
    283 
    284         for name, spec in self.parameters.items():

TypeError: Argument to input 'table' is not a subtype of FeatureTable[Frequency] % Properties(['uniform-sampling']).
```

---
deletable: true
editable: true
...

#### A limitation of OTU counting <link src='a3bd5f'/>

---
deletable: true
editable: true
...

Imagine that we have the same table, but some additional information about the OTUs in the table. Specifically, we've computed the following phylogenetic tree. And, for the sake of illustration, imagine that we've also assigned taxonomy to each of the OTUs and found that our samples contain representatives from the archaea, bacteria, and eukaryotes (their labels begin with `A`, `B`, and `E`, respectively).

---
deletable: true
editable: true
...

<img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/pd_calc_tree.png">

---
deletable: true
editable: true
...

Pairing this with the table we defined above (displayed again in the cell below), given what you now know about these OTUs, which would you consider the most diverse? Are you happy with the $\alpha$ diversity conclusion that you obtained when computing the number of observed OTUs in each sample?

---
deletable: true
editable: true
...

```python
>>> table2_rarefied.view(pd.DataFrame)
      B1     B2     B3     B4     A1    E2
A  100.0  100.0  300.0    0.0    0.0   0.0
B   78.0  160.0   90.0  172.0    0.0   0.0
C  285.0    0.0    0.0    0.0  169.0  46.0
```

---
deletable: true
editable: true
...

### Phylogenetic Diversity (PD) <link src='c407f8'/>

---
deletable: true
editable: true
...

Phylogenetic Diversity (PD) is a metric that was developed by Dan Faith in the early 1990s (find the original paper [here](http://www.sciencedirect.com/science/article/pii/0006320792912013)). Like many of the measures that are used in microbial community ecology, it wasn't initially designed for studying microbial communities, but rather communities of "macro-organisms" (macrobes?). Some of these metrics, including PD, do translate well to microbial community analysis, while some don't translate as well. (For an illustration of the effect of sequencing error on PD, where it is handled well, versus its effect on the Chao1 metric, where it is handled less well, see Figure 1 of [Reeder and Knight (2010)](http://www.nature.com/nmeth/journal/v7/n9/full/nmeth0910-668b.html)).

---
deletable: true
editable: true
...

PD is relatively simple to calculate. It is computed simply as the sum of the branch length in a phylogenetic tree that is "covered" or represented in a given sample. Let's look at an example to see how this works.

---
deletable: true
editable: true
...

First, let's define a phylogenetic tree using the Newick format (which is described [here](http://evolution.genetics.washington.edu/phylip/newicktree.html), and more formally defined [here](http://evolution.genetics.washington.edu/phylip/newick_doc.html)). We'll then load that up using [scikit-bio](http://scikit-bio.org)'s [TreeNode](http://scikit-bio.org/generated/skbio.core.tree.TreeNode.html#skbio.core.tree.TreeNode) object.

---
deletable: true
editable: true
...

```python
>>> import os.path
>>> import tempfile
...
>>> newick_tree = ('(((B1:0.2,B2:0.3):0.3,((B3:0.5,B4:0.3):0.2,B5:0.9):0.3):0.35,'
...                '(((A1:0.2,A2:0.3):0.3,(E1:0.3,E2:0.4):0.7):0.2):0.05)root;')
...
>>> with tempfile.TemporaryDirectory() as newick_dir:
...     newick_fp = os.path.join(newick_dir, 'tree.nwk')
...     open(newick_fp, 'w').write(newick_tree)
...     unrooted_tree = qiime2.Artifact.import_data('Phylogeny[Unrooted]', newick_dir)
```

```python
>>> import qiime2.plugins.phylogeny
>>> midpoint_root_result = qiime2.plugins.phylogeny.actions.midpoint_root(unrooted_tree)
>>> rooted_tree = midpoint_root_result.rooted_tree
```

```python
>>> import skbio
>>> rooted_tree.view(skbio.TreeNode)
<TreeNode, name: root, internal node count: 9, tips count: 9>
```

---
deletable: true
editable: true
...

I'll now define a couple of functions that we'll use to compute PD.

---
deletable: true
editable: true
...

```python
>>> def get_observed_nodes(tree, table, sample_id, verbose=False):
...     observed_otus = [obs_id for obs_id in table.columns
...                 if table[obs_id][sample_id] > 0]
...     observed_nodes = set()
...     # iterate over the observed OTUs
...     for otu in observed_otus:
...         t = tree.find(otu)
...         observed_nodes.add(t)
...         if verbose:
...             print(t.name, t.length, end=' ')
...         for internal_node in t.ancestors():
...             if internal_node.length is None:
...                 # we've hit the root
...                 if verbose:
...                     print('')
...             else:
...                 if verbose and internal_node not in observed_nodes:
...                     print(internal_node.length, end=' ')
...                 observed_nodes.add(internal_node)
...     return observed_nodes
...
>>> def phylogenetic_diversity(tree, table, sample_id, verbose=False):
...     observed_nodes = get_observed_nodes(tree, table, sample_id, verbose=verbose)
...     result = sum(o.length for o in observed_nodes)
...     return result
```

---
deletable: true
editable: true
...

And then apply those to compute the PD of our three samples. For each computation, we're also printing out the branch lengths of the branches that are observed *for the first time* when looking at a given OTU. When computing PD, we include the length of each branch only one time.

---
deletable: true
editable: true
...

```python
>>> pd_A = phylogenetic_diversity(rooted_tree.view(skbio.TreeNode),
...                               table2_rarefied.view(pd.DataFrame),
...                               'A',
...                               verbose=True)
>>> print(pd_A)
B1 0.2 0.3 0.25 
B2 0.3 
B3 0.5 0.2 0.3 
2.05
```

---
deletable: true
editable: true
...

```python
>>> pd_B = phylogenetic_diversity(rooted_tree.view(skbio.TreeNode),
...                               table2_rarefied.view(pd.DataFrame),
...                               'B',
...                               verbose=True)
>>> print(pd_B)
B1 0.2 0.3 0.25 
B2 0.3 
B3 0.5 0.2 0.3 
B4 0.3 
2.35
```

---
deletable: true
editable: true
...

```python
>>> pd_C = phylogenetic_diversity(rooted_tree.view(skbio.TreeNode),
...                               table2_rarefied.view(pd.DataFrame),
...                               'C',
...                               verbose=True)
>>> print(pd_C)
B1 0.2 0.3 0.25 
A1 0.2 0.3 0.2 0.05 0.1 
E2 0.4 0.7 
2.7
```

```python
>>> alpha_phylogenetic_result = qiime2.plugins.diversity.actions.alpha_phylogenetic(
...                                 table=table2_rarefied,
...                                 phylogeny=rooted_tree,
...                                 metric='faith_pd')
>>> faith_pd = alpha_phylogenetic_result.alpha_diversity
>>> faith_pd.view(pd.Series)
A    2.05
B    2.35
C    2.70
Name: faith_pd, dtype: float64
```

---
deletable: true
editable: true
...

How does this result compare to what we observed above with the Observed OTUs metric? Based on your knowledge of biology, which do you think is a better representation of the relative diversities of these samples?

---
deletable: true
editable: true
...

## Measuring beta diversity <link src='cb4608'/>

---
deletable: true
editable: true
...

$\beta$-diversity (canonically pronounced *beta diversity*) refers to **between sample diversity**, and is typically used to answer questions of the form: is sample $A$ more similar in composition to sample $B$ or sample $C$? In this section we'll explore two (of tens or hundreds) of metrics for computing pairwise dissimilarity of samples to estimate $\beta$ diversity.

---
deletable: true
editable: true
...

### Distance metrics <link src='eac92f'/>

---
deletable: true
editable: true
...

#### Bray-Curtis <link src='dac934'/>

---
deletable: true
editable: true
...

The first metric that we'll look at is a quantitative non-phylogenetic $\beta$ diversity metric called Bray-Curtis. The Bray-Curtis dissimilarity between a pair of samples, $j$ and $k$, is defined as follows:

---
deletable: true
editable: true
...

$BC_{jk} = \frac{ \sum_{i} | X_{ij} - X_{ik}|} {\sum_{i} (X_{ij} + X_{ik})}$

---
deletable: true
editable: true
...

$i$ : feature (e.g., OTUs)

---
deletable: true
editable: true
...

$X_{ij}$ : frequency of feature $i$ in sample $j$

---
deletable: true
editable: true
...

$X_{ik}$ : frequency of feature $i$ in sample $k$

---
deletable: true
editable: true
...

This could be implemented in python as follows:

---
deletable: true
editable: true
...

```python
>>> def bray_curtis_distance(table, sample1_id, sample2_id):
...     numerator = 0
...     denominator = 0
...     sample1_counts = table.loc[sample1_id]
...     sample2_counts = table.loc[sample2_id]
...     for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
...         numerator += abs(sample1_count - sample2_count)
...         denominator += sample1_count + sample2_count
...     return numerator / denominator
```

---
deletable: true
editable: true
...

```python
>>> table2_rarefied.view(pd.DataFrame)
      B1     B2     B3     B4     A1    E2
A  100.0  100.0  300.0    0.0    0.0   0.0
B   78.0  160.0   90.0  172.0    0.0   0.0
C  285.0    0.0    0.0    0.0  169.0  46.0
```

---
deletable: true
editable: true
...

Let's now apply this to some pairs of samples:

---
deletable: true
editable: true
...

```python
>>> print(bray_curtis_distance(table2_rarefied.view(pd.DataFrame), 'A', 'B'))
0.464
```

---
deletable: true
editable: true
...

```python
>>> print(bray_curtis_distance(table2_rarefied.view(pd.DataFrame), 'A', 'C'))
0.8
```

---
deletable: true
editable: true
...

```python
>>> print(bray_curtis_distance(table2_rarefied.view(pd.DataFrame), 'B', 'C'))
0.844
```

---
deletable: true
editable: true
...

```python
>>> print(bray_curtis_distance(table2_rarefied.view(pd.DataFrame), 'A', 'A'))
0.0
```

---
deletable: true
editable: true
...

```python
>>> print(bray_curtis_distance(table2_rarefied.view(pd.DataFrame), 'C', 'B'))
0.844
```

---
deletable: true
editable: true
...

```python
>>> beta_result = qiime2.plugins.diversity.actions.beta(table=table2_rarefied,
...                                                     metric='braycurtis')
>>> bray_curtis_distance_matrix = beta_result.distance_matrix
>>> _ = bray_curtis_distance_matrix.view(skbio.DistanceMatrix).plot(cmap='viridis')
```

```python
>>> print(bray_curtis_distance_matrix.view(skbio.DistanceMatrix))
3x3 distance matrix
IDs:
'A', 'B', 'C'
Data:
[[ 0.     0.464  0.8  ]
 [ 0.464  0.     0.844]
 [ 0.8    0.844  0.   ]]
```

---
deletable: true
editable: true
...

#### Unweighted UniFrac <link src='2682aa'/>

**THIS TEXT NEEDS A LOT OF UPDATES TO USE THE NEW TREE!!**

---
deletable: true
editable: true
...

Just as phylogenetic alpha diversity metrics can be more informative than non-phylogenetic alpha diversity metrics, phylogenetic beta diversity metrics offer advantages over non-phylogenetic metrics such as Bray-Curtis. The most widely applied phylogenetic beta diversity metric as of this writing is unweighted UniFrac. UniFrac was initially presented in [Lozupone and Knight, 2005, Applied and Environmental Microbiology](http://aem.asm.org/content/71/12/8228.abstract), and has been widely applied in microbial ecology since (and the illustration of UniFrac computation presented below is derived from a similar example originally developed by Lozupone and Knight).

---
deletable: true
editable: true
...

The unweighted UniFrac distance between a pair of samples `A` and `B` is defined as follows:

---
deletable: true
editable: true
...

$U_{AB} = \frac{unique}{observed}$

---
deletable: true
editable: true
...

where:

---
deletable: true
editable: true
...

$unique$ : the unique branch length, or branch length that only leads to OTU(s) observed in sample $A$ or sample $B$

---
deletable: true
editable: true
...

$observed$ : the total branch length observed in either sample $A$ or sample $B$

---
deletable: true
editable: true
...

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_d0.png" align=right/></div>

---
deletable: true
editable: true
...

To illustrate how UniFrac distances are computed, before we get into actually computing them, let's look at a few examples. In these examples, imagine that we're determining the pairwise UniFrac distance between two samples: a red sample, and a blue sample. If a red box appears next to an OTU, that indicates that it's observed in the red sample; if a blue box appears next to the OTU, that indicates that it's observed in the blue sample; if a red and blue box appears next to the OTU, that indicates that the OTU is present in both samples; and if no box is presented next to the OTU, that indicates that it's present in neither sample.

---
deletable: true
editable: true
...

To compute the UniFrac distance between a pair of samples, we need to know the sum of the branch length that was observed in either sample (the *observed* branch length), and the sum of the branch length that was observed only in a single sample (the *unique* branch length). In these examples, we color all of the *observed* branch length. Branch length that is unique to the red sample is red, branch length that is unique to the blue sample is blue, and branch length that is observed in both samples is purple. Unobserved branch length is black (as is the vertical branches, as those don't contribute to branch length - they are purely for visual presentation).

---
deletable: true
editable: true
...

In the tree on the right, all of the OTUs that are observed in either sample are observed in both samples. As a result, all of the observed branch length is purple. The unique branch length in this case is zero, so **we have a UniFrac distance of 0 between the red and blue samples**.

---
deletable: true
editable: true
...

<hr>

---
deletable: true
editable: true
...

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_d1.png" align=right/></div>

---
deletable: true
editable: true
...

On the other end of the spectrum, in the second tree, all of the OTUs in the tree are observed either in the red sample, or in the blue sample. All of the observed branch length in the tree is either red or blue, meaning that if you follow a branch out to the tips, you will observe only red or blue samples. In this case the unique branch length is equal to the observed branch length, so **we have a UniFrac distance of 1 between the red and blue samples**.

---
deletable: true
editable: true
...

<hr>

---
deletable: true
editable: true
...

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_d0.5.png" align=right/></div>

---
deletable: true
editable: true
...

Finally, most of the time we're somewhere in the middle. In this tree, some of our branch length is unique, and some is not. For example, OTU 1 is only observed in our red sample, so the terminal branch leading to OTU 1 is red (i.e., unique to the red sample). OTU 2 is only observed in our blue sample, so the terminal branch leading to OTU 2 is blue (i.e., unique to the blue sample). However, the internal branch leading to the node connecting OTU 1 and OTU 2 leads to OTUs observed in both the red and blue samples (i.e., OTU 1 and OTU 2), so is purple (i.e, observed branch length, but not unique branch length). In this case, **we have an intermediate UniFrac distance between the red and blue samples, maybe somewhere around 0.5**.

---
deletable: true
editable: true
...

<hr>
<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_with_distances.png" align=right/></div>

---
deletable: true
editable: true
...

Let's now compute the Unweighted UniFrac distances between some samples. Imagine we have the following tree, paired with our table below (printed below, for quick reference).

---
deletable: true
editable: true
...

```python
>>> table2_rarefied.view(pd.DataFrame)
      B1     B2     B3     B4     A1    E2
A  100.0  100.0  300.0    0.0    0.0   0.0
B   78.0  160.0   90.0  172.0    0.0   0.0
C  285.0    0.0    0.0    0.0  169.0  46.0
```

---
deletable: true
editable: true
...

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_with_distances_ab.png" align=right/></div>

---
deletable: true
editable: true
...

First, let's compute the unweighted UniFrac distance between samples $A$ and $B$. The *unweighted* in *unweighted UniFrac* means that this is a qualitative diversity metric, meaning that we don't care about the abundances of the OTUs, only whether they are present in a given sample ($frequency > 0$) or not present ($frequency = 0$).

---
deletable: true
editable: true
...

Start at the top right branch in the tree, and for each branch, determine if the branch is observed, and if so, if it is also unique. If it is observed then you add its length to your observed branch length. If it is observed and unique, then you also add its length to your unique branch length.

---
deletable: true
editable: true
...

For samples $A$ and $B$, I get the following (in the tree on the right, red branches are those observed in $A$, blue branches are those observed in $B$, and purple are observed in both):

---
deletable: true
editable: true
...

$unique_{ab} = 0.5 + 0.75 = 1.25$

---
deletable: true
editable: true
...

$observed_{ab} = 0.5 + 0.5 + 0.5 + 1.0 + 1.25 + 0.75 + 0.75 = 5.25$

---
deletable: true
editable: true
...

$uu_{ab} = \frac{unique_{ab}}{observed_{ab}} = \frac{1.25}{5.25} = 0.238$

---
deletable: true
editable: true
...

As an exercise, now compute the UniFrac distances between samples $B$ and $C$, and samples $A$ and $C$, using the above table and tree. When I do this, I get the following distance matrix.

---
deletable: true
editable: true
...

```python
>>> beta_result = qiime2.plugins.diversity.actions.beta_phylogenetic(table=table2_rarefied,
...                                                                  phylogeny=rooted_tree,
...                                                                  metric='unweighted_unifrac')
>>> unweighted_unifrac_distance_matrix = beta_result.distance_matrix
>>> _ = unweighted_unifrac_distance_matrix.view(skbio.DistanceMatrix).plot(cmap='viridis')
```

```python
>>> print(unweighted_unifrac_distance_matrix.view(skbio.DistanceMatrix))
3x3 distance matrix
IDs:
'A', 'B', 'C'
Data:
[[ 0.          0.12765957  0.8125    ]
 [ 0.12765957  0.          0.8255814 ]
 [ 0.8125      0.8255814   0.        ]]
```

---
deletable: true
editable: true
...

```python
>>> def unweighted_unifrac(tree, table, sample_id1, sample_id2, verbose=False):
...     observed_nodes1 = get_observed_nodes(tree, table, sample_id1, verbose=verbose)
...     observed_nodes2 = get_observed_nodes(tree, table, sample_id2, verbose=verbose)
...     observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2)
...     shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2)
...     unique_branch_length = observed_branch_length - shared_branch_length
...     unweighted_unifrac = unique_branch_length / observed_branch_length
...     return unweighted_unifrac
...
>>> print(unweighted_unifrac(rooted_tree.view(skbio.TreeNode), table2_rarefied.view(pd.DataFrame), 'A', 'B'))
>>> print(unweighted_unifrac(rooted_tree.view(skbio.TreeNode), table2_rarefied.view(pd.DataFrame), 'A', 'C'))
>>> print(unweighted_unifrac(rooted_tree.view(skbio.TreeNode), table2_rarefied.view(pd.DataFrame), 'B', 'C'))
0.1276595744680852
0.8125
0.8255813953488372
```

---
deletable: true
editable: true
...

#### Even sampling <link src='200e13'/>

---
deletable: true
editable: true
...

**TODO**: Add discussion on necessity of even sampling

---
deletable: true
editable: true
...

### Interpreting distance matrices <link src='2be688'/>

---
deletable: true
editable: true
...

In the previous section we computed distance matrices that contained the pairwise distances between a few samples. You can look at those distance matrices and get a pretty good feeling for what the patterns are. For example, what are the most similar samples? What are the most dissimilar samples?

---
deletable: true
editable: true
...

What if instead of three samples though, we had more. Here's a screenshot from a distance matrix containing data on 105 samples (this is just the first few rows and columns):

---
deletable: true
editable: true
...

<img src='https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/example_big_dm.png', width=800>

---
deletable: true
editable: true
...

Do you have a good feeling for the patterns here? What are the most similar samples? What are the most dissimilar samples?

---
deletable: true
editable: true
...

Chances are, you can't just squint at that table and understand what's going on (but if you can, I'm hiring!). The problem is exacerbated by the fact that in modern microbial ecology studies we may have thousands or tens of thousands of samples, not "just" hundreds as in the table above. We need tools to help us take these raw distances and convert them into something that we can interpret. In this section we'll look at some techniques, one of which we've covered previously, that will help us interpret large distance matrices.

---
deletable: true
editable: true
...

<hr>

---
deletable: true
editable: true
...

One excellent paper that includes a comparison of several different strategies for interpreting beta diversity results is [Costello *et al.* Science (2009) Bacterial Community Variation in Human Body Habitats Across Space and Time](https://www.sciencemag.org/content/326/5960/1694.full). In this study, the authors collected microbiome samples from 7 human subjects at about 25 sites on their bodies, at four different points in time.

---
deletable: true
editable: true
...

Figure 1 shows several different approaches for comparing the resulting UniFrac distance matrix (this image is linked from the *Science* journal website - copyright belongs to *Science*):

---
deletable: true
editable: true
...

<img src="https://www.sciencemag.org/content/326/5960/1694/F1.large.jpg" width=800>

---
deletable: true
editable: true
...

Let's generate a small distance matrix representing just a few of these body sites, and figure out how we'd generate and interpret each of these visualizations. The values in the distance matrix below are a subset of the unweighted UniFrac distance matrix representing two samples each from three body sites from the Costello *et al.* (2009) study.

---
deletable: true
editable: true
...

```python
>>> sample_ids = ['A', 'B', 'C', 'D', 'E', 'F']
>>> _columns = ['body site', 'individual']
>>> _md = [['gut', 'subject 1'],
...        ['gut', 'subject 2'],
...        ['tongue', 'subject 1'],
...        ['tongue', 'subject 2'],
...        ['skin', 'subject 1'],
...        ['skin', 'subject 2']]
...
>>> human_microbiome_sample_md = pd.DataFrame(_md, index=sample_ids, columns=_columns)
>>> human_microbiome_sample_md
  body site individual
A       gut  subject 1
B       gut  subject 2
C    tongue  subject 1
D    tongue  subject 2
E      skin  subject 1
F      skin  subject 2
```

---
deletable: true
editable: true
...

```python
>>> dm_data = array([[0.00, 0.35, 0.83, 0.83, 0.90, 0.90],
...                  [0.35, 0.00, 0.86, 0.85, 0.92, 0.91],
...                  [0.83, 0.86, 0.00, 0.25, 0.88, 0.87],
...                  [0.83, 0.85, 0.25, 0.00, 0.88, 0.88],
...                  [0.90, 0.92, 0.88, 0.88, 0.00, 0.50],
...                  [0.90, 0.91, 0.87, 0.88, 0.50, 0.00]])
...
>>> human_microbiome_dm = DistanceMatrix(dm_data, sample_ids)
>>> print(human_microbiome_dm)
6x6 distance matrix
IDs:
'A', 'B', 'C', 'D', 'E', 'F'
Data:
[[ 0.    0.35  0.83  0.83  0.9   0.9 ]
 [ 0.35  0.    0.86  0.85  0.92  0.91]
 [ 0.83  0.86  0.    0.25  0.88  0.87]
 [ 0.83  0.85  0.25  0.    0.88  0.88]
 [ 0.9   0.92  0.88  0.88  0.    0.5 ]
 [ 0.9   0.91  0.87  0.88  0.5   0.  ]]
```

---
deletable: true
editable: true
...

#### Distribution plots and comparisons <link src='8fcf92'/>

---
deletable: true
editable: true
...

First, let's look at the analysis presented in panels E and F. Instead of generating bar plots here, we'll generate box plots as these are more informative (i.e., they provide a more detailed summary of the distribution being investigated). One important thing to notice here is the central role that the sample metadata plays in the visualization. If we just had our sample ids (i.e., letters ``A`` through ``F``) we wouldn't be able to group distances into *within* and *between* sample type categories, and we therefore couldn't perform the comparisons we're interested in.

---
deletable: true
editable: true
...

```python
>>> def within_between_category_distributions(dm, md, md_category):
...     within_category_distances = []
...     between_category_distances = []
...     for i, sample_id1 in enumerate(dm.ids):
...         sample_md1 = md[md_category][sample_id1]
...         for sample_id2 in dm.ids[:i]:
...             sample_md2 = md[md_category][sample_id2]
...             if sample_md1 == sample_md2:
...                 within_category_distances.append(dm[sample_id1, sample_id2])
...             else:
...                 between_category_distances.append(dm[sample_id1, sample_id2])
...     return within_category_distances, between_category_distances
```

---
deletable: true
editable: true
...

```python
>>> within_category_distances, between_category_distances = within_between_category_distributions(human_microbiome_dm, human_microbiome_sample_md, "body site")
>>> print(within_category_distances)
>>> print(between_category_distances)
[0.34999999999999998, 0.25, 0.5]
[0.82999999999999996, 0.85999999999999999, 0.82999999999999996, 0.84999999999999998, 0.90000000000000002, 0.92000000000000004, 0.88, 0.88, 0.90000000000000002, 0.91000000000000003, 0.87, 0.88]
```

---
deletable: true
editable: true
...

```python
>>> import seaborn as sns
>>> ax = sns.boxplot(data=[within_category_distances, between_category_distances])
>>> ax.set_xticklabels(['same body habitat', 'different body habitat'])
>>> ax.set_ylabel('Unweighted UniFrac Distance')
>>> _ = ax.set_ylim(0.0, 1.0)
/Users/caporaso/miniconda3/envs/iab/lib/python3.4/site-packages/matplotlib/__init__.py:872: UserWarning: axes.color_cycle is deprecated and replaced with axes.prop_cycle; please use the latter.
  warnings.warn(self.msg_depr % (key, alt_key))
/Users/caporaso/miniconda3/envs/iab/lib/python3.4/site-packages/matplotlib/__init__.py:892: UserWarning: axes.color_cycle is deprecated and replaced with axes.prop_cycle; please use the latter.
  warnings.warn(self.msg_depr % (key, alt_key))
```

---
deletable: true
editable: true
...

```python
>>> from skbio.stats.distance import anosim
>>> anosim(human_microbiome_dm, human_microbiome_sample_md, 'body site')
method name               ANOSIM
test statistic name            R
sample size                    6
number of groups               3
test statistic                 1
p-value                    0.054
number of permutations       999
Name: ANOSIM results, dtype: object
```

---
deletable: true
editable: true
...

If we run through these same steps, but base our analysis on a different metadata category where we don't expect to see any significant clustering, you can see that we no longer get a significant result.

---
deletable: true
editable: true
...

```python
>>> within_category_distances, between_category_distances = within_between_category_distributions(human_microbiome_dm, human_microbiome_sample_md, "individual")
>>> print(within_category_distances)
>>> print(between_category_distances)
[0.82999999999999996, 0.84999999999999998, 0.90000000000000002, 0.88, 0.91000000000000003, 0.88]
[0.34999999999999998, 0.85999999999999999, 0.82999999999999996, 0.25, 0.92000000000000004, 0.88, 0.90000000000000002, 0.87, 0.5]
```

---
deletable: true
editable: true
...

```python
>>> ax = sns.boxplot(data=[within_category_distances, between_category_distances])
>>> ax.set_xticklabels(['same person', 'different person'])
>>> ax.set_ylabel('Unweighted UniFrac Distance')
>>> _ = ax.set_ylim(0.0, 1.0)
/Users/caporaso/miniconda3/envs/iab/lib/python3.4/site-packages/matplotlib/__init__.py:892: UserWarning: axes.color_cycle is deprecated and replaced with axes.prop_cycle; please use the latter.
  warnings.warn(self.msg_depr % (key, alt_key))
```

---
deletable: true
editable: true
...

```python
>>> anosim(human_microbiome_dm, human_microbiome_sample_md, 'individual')
method name                 ANOSIM
test statistic name              R
sample size                      6
number of groups                 2
test statistic           -0.333333
p-value                      0.889
number of permutations         999
Name: ANOSIM results, dtype: object
```

---
deletable: true
editable: true
...

Why do you think the distribution of distances between people has a greater range than the distribution of distances within people in this particular example?

---
deletable: true
editable: true
...

Here we used ANOSIM testing whether our with and between category groups differ. This test is specifically designed for distance matrices, and it accounts for the fact that the values are not independent of one another. For example, if one of our samples was very different from all of the others, all of the distances associated with that sample would be large. It's very important to choose the appropriate statistical test to use. One free resource for helping you do that is [*The Guide to Statistical Analysis in Microbial Ecology (GUSTAME)*](http://mb3is.megx.net/gustame). If you're getting started in microbial ecology, I recommend spending some time studying GUSTAME.

---
deletable: true
editable: true
...

#### Hierarchical clustering <link src='09f456'/>

---
deletable: true
editable: true
...

Next, let's look at a hierarchical clustering analysis, similar to that presented in panel G above. Here I'm applying the UPGMA functionality implemented in [SciPy](http://www.scipy.org/scipylib/index.html) to generate a tree which we visualize with a dendrogram. However the tips in this tree don't represent sequences or OTUs, like they did when we [covered UPGMA for phylogenetic reconstruction](alias://73d028), but instead they represent samples, and samples with a smaller branch length between them are more similar in composition than samples with a longer branch length between them. (Remember that only horizontal branch length is counted - vertical branch length is just to aid in the organization of the dendrogram.)

---
deletable: true
editable: true
...

```python
>>> from scipy.cluster.hierarchy import average, dendrogram
>>> lm = average(human_microbiome_dm.condensed_form())
>>> d = dendrogram(lm, labels=human_microbiome_dm.ids, orientation='right',
...                link_color_func=lambda x: 'black')
```

---
deletable: true
editable: true
...

Again, we can see how the data really only becomes interpretable in the context of metadata:

---
deletable: true
editable: true
...

```python
>>> labels = [human_microbiome_sample_md['body site'][sid] for sid in sample_ids]
>>> d = dendrogram(lm, labels=labels, orientation='right',
...                link_color_func=lambda x: 'black')
```

---
deletable: true
editable: true
...

```python
>>> labels = [human_microbiome_sample_md['individual'][sid] for sid in sample_ids]
>>> d = dendrogram(lm, labels=labels, orientation='right',
...                link_color_func=lambda x: 'black')
```

---
deletable: true
editable: true
...

### Ordination <link src='b1cdbe'/>

---
deletable: true
editable: true
...

Finally, let's look at ordination, similar to that presented in panels A-D. The basic idea behind ordination is dimensionality reduction: we want to take high-dimensionality data (a distance matrix) and represent that in a few (usually two or three) dimensions. As humans, we're very bad at interpreting high dimensionality data directly: with ordination, we can take an $n$-dimensional data set (e.g., a distance matrix of shape $n \times n$, representing the distances between $n$ biological samples) and reduce that to a 2-dimensional scatter plot similar to that presented in panels A-D above.

---
deletable: true
editable: true
...

Ordination is a technique that is widely applied in ecology and in bioinformatics, but the math behind some of the methods such as *Principal Coordinates Analysis* is fairly complex, and as a result I've found that these methods are a black box for a lot of people. Possibly the most simple ordination technique is one called Polar Ordination. Polar Ordination is not widely applied because it has some inconvenient features, but I find that it is useful for introducing the idea behind ordination. Here we'll work through a simple implementation of ordination to illustrate the process, which will help us to interpret ordination plots. In practice, you will use existing software, such as [scikit-bio](http://scikit-bio.org)'s [ordination module](http://scikit-bio.org/maths.stats.ordination.html).

---
deletable: true
editable: true
...

An excellent site for learning more about ordination is [Michael W. Palmer's Ordination Methods page](http://ordination.okstate.edu/).

---
deletable: true
editable: true
...

#### Polar ordination <link src='538e18'/>

---
deletable: true
editable: true
...

First, let's print our distance matrix again so we have it nearby.

---
deletable: true
editable: true
...

```python
>>> print(human_microbiome_dm)
6x6 distance matrix
IDs:
'A', 'B', 'C', 'D', 'E', 'F'
Data:
[[ 0.    0.35  0.83  0.83  0.9   0.9 ]
 [ 0.35  0.    0.86  0.85  0.92  0.91]
 [ 0.83  0.86  0.    0.25  0.88  0.87]
 [ 0.83  0.85  0.25  0.    0.88  0.88]
 [ 0.9   0.92  0.88  0.88  0.    0.5 ]
 [ 0.9   0.91  0.87  0.88  0.5   0.  ]]
```

---
deletable: true
editable: true
...

Polar ordination works in a few steps:

---
deletable: true
editable: true
...

**Step 1.** Identify the largest distance in the distance matrix.

---
deletable: true
editable: true
...

**Step 2.** Define a line, with the two samples contributing to that distance defining the endpoints.

---
deletable: true
editable: true
...

**Step 3.** Compute the location of each other sample on that axis as follows:

---
deletable: true
editable: true
...

$a = \frac{D^2 + D1^2 - D2^2}{2 \times D}$

---
deletable: true
editable: true
...

where:

---
deletable: true
editable: true
...

$D$ is distance between the endpoints

---
deletable: true
editable: true
...

$D1$ is distance between the current sample and endpoint 1

---
deletable: true
editable: true
...

$D2$ is distance between sample and endpoint 2.

---
deletable: true
editable: true
...

**Step 4.** Find the next largest distance that could be used to define an *uncorrelated axis*. (This step can be labor-intensive to do by hand - usually you would compute all of the axes, along with correlation scores. I'll pick one for the demo, and we'll wrap up by looking at all of the axes.)

---
deletable: true
editable: true
...

Here is what steps 2 and 3 look like in Python:

---
deletable: true
editable: true
...

```python
>>> def compute_axis_values(dm, endpoint1, endpoint2):
...     d = dm[endpoint1, endpoint2]
...     result = {endpoint1: 0, endpoint2: d}
...     non_endpoints = set(dm.ids) - set([endpoint1, endpoint2])
...     for e in non_endpoints:
...         d1 = dm[endpoint1, e]
...         d2 = dm[endpoint2, e]
...         result[e] = (d**2 + d1**2 - d2**2) / (2 * d)
...     return d, [result[e] for e in dm.ids]
```

---
deletable: true
editable: true
...

```python
>>> d, a1_values = compute_axis_values(human_microbiome_dm, 'B', 'E')
>>> for sid, a1_value in zip(human_microbiome_dm.ids, a1_values):
...     print(sid, a1_value)
A 0.0863586956522
B 0
C 0.441086956522
D 0.431793478261
E 0.92
F 0.774184782609
```

---
deletable: true
editable: true
...

```python
>>> d, a2_values = compute_axis_values(human_microbiome_dm, 'D', 'E')
>>> for sid, a2_value in zip(human_microbiome_dm.ids, a2_values):
...     print(sid, a2_value)
A 0.371193181818
B 0.369602272727
C 0.0355113636364
D 0
E 0.88
F 0.737954545455
```

---
deletable: true
editable: true
...

```python
>>> from pylab import scatter
>>> ord_plot = scatter(a1_values, a2_values, s=40)
```

---
deletable: true
editable: true
...

And again, let's look at how including metadata helps us to interpret our results.

---
deletable: true
editable: true
...

First, we'll color the points by the body habitat that they're derived from:

---
deletable: true
editable: true
...

```python
>>> colors = {'tongue': 'red', 'gut':'yellow', 'skin':'blue'}
>>> c = [colors[human_microbiome_sample_md['body site'][e]] for e in human_microbiome_dm.ids]
>>> ord_plot = scatter(a1_values, a2_values, s=40, c=c)
```

---
deletable: true
editable: true
...

And next we'll color the samples by the person that they're derived from. Notice that this plot and the one above are identical except for coloring. Think about how the colors (and therefore the sample metadata) help you to interpret these plots.

---
deletable: true
editable: true
...

```python
>>> person_colors = {'subject 1': 'red', 'subject 2':'yellow'}
>>> person_c = [person_colors[human_microbiome_sample_md['individual'][e]] for e in human_microbiome_dm.ids]
>>> ord_plot = scatter(a1_values, a2_values, s=40, c=person_c)
```

---
deletable: true
editable: true
...

#### Determining the most important axes in polar ordination <link src='fb483b'/>

---
deletable: true
editable: true
...

Generally, you would compute the polar ordination axes for all possible axes. You could then order the axes by which represent the largest differences in sample composition, and the lowest correlation with previous axes. This might look like the following:

---
deletable: true
editable: true
...

```python
>>> from scipy.stats import spearmanr
...
>>> data = []
>>> for i, sample_id1 in enumerate(human_microbiome_dm.ids):
...     for sample_id2 in human_microbiome_dm.ids[:i]:
...         d, axis_values = compute_axis_values(human_microbiome_dm, sample_id1, sample_id2)
...         r, p = spearmanr(a1_values, axis_values)
...         data.append((d, abs(r), sample_id1, sample_id2, axis_values))
...
>>> data.sort()
>>> data.reverse()
>>> for i, e in enumerate(data):
...     print("axis %d:" % i, end=' ')
...     print("\t%1.3f\t%1.3f\t%s\t%s" % e[:4])
axis 0: 	0.920	1.000	E	B
axis 1: 	0.910	0.943	F	B
axis 2: 	0.900	0.928	E	A
axis 3: 	0.900	0.886	F	A
axis 4: 	0.880	0.543	E	D
axis 5: 	0.880	0.429	F	D
axis 6: 	0.880	0.429	E	C
axis 7: 	0.870	0.371	F	C
axis 8: 	0.860	0.543	C	B
axis 9: 	0.850	0.486	D	B
axis 10: 	0.830	0.429	C	A
axis 11: 	0.830	0.406	D	A
axis 12: 	0.500	0.232	F	E
axis 13: 	0.350	0.143	B	A
axis 14: 	0.250	0.493	D	C
```

---
deletable: true
editable: true
...

So why do we care about axes being uncorrelated? And why do we care about explaining a lot of the variation? Let's look at a few of these plots and see how they compare to the plots above, where we compared axes 1 and 4.

---
deletable: true
editable: true
...

```python
>>> ord_plot = scatter(data[0][4], data[1][4], s=40, c=c)
```

---
deletable: true
editable: true
...

```python
>>> ord_plot = scatter(data[0][4], data[13][4], s=40, c=c)
```

---
deletable: true
editable: true
...

```python
>>> ord_plot = scatter(data[0][4], data[14][4], s=40, c=c)
```

---
deletable: true
editable: true
...

#### Interpreting ordination plots <link src='40e0a6'/>

---
deletable: true
editable: true
...

There are a few points that are important to keep in mind when interpreting ordination plots. Review each one of these in the context of polar ordination to figure out the reason for each.

---
deletable: true
editable: true
...

**Directionality of the axes is not important (e.g., up/down/left/right)**

---
deletable: true
editable: true
...

One thing that you may have notices as you computed the polar ordination above is that the method is *not symmetric*: in other words, the axis values for axis $EB$ are different than for axis $BE$. In practice though, we derive the same conclusions regardless of how we compute that axis: in this example, that samples cluster by body site.

---
deletable: true
editable: true
...

```python
>>> d, a1_values = compute_axis_values(human_microbiome_dm, 'E', 'B')
>>> d, a2_values = compute_axis_values(human_microbiome_dm, 'E', 'D')
>>> d, alt_a1_values = compute_axis_values(human_microbiome_dm, 'B', 'E')
```

---
deletable: true
editable: true
...

```python
>>> ord_plot = scatter(a1_values, a2_values, s=40, c=c)
```

---
deletable: true
editable: true
...

```python
>>> ord_plot = scatter(alt_a1_values, a2_values, s=40, c=c)
```

---
deletable: true
editable: true
...

Some other important features:

---
deletable: true
editable: true
...

* Numerical scale of the axis is generally not useful
* The order of axes is generally important (first axis explains the most variation, second axis explains the second most variation, ...)
* Most techniques result in uncorrelated axes.
* Additional axes can be generated (third, fourth, ...)

---
deletable: true
editable: true
...

## Tools for using ordination in practice: scikit-bio, pandas, and matplotlib <link src='098854'/>

---
deletable: true
editable: true
...

As I mentioned above, polar ordination isn't widely used in practice, but the features that it illustrates are common to ordination methods. One of the most widely used ordination methods used to study biological diversity is Principal Coordinates Analysis or PCoA, which is implemented in [scikit-bio](http://scikit-bio.org/)'s [``ordination`` module](http://scikit-bio.org/maths.stats.ordination.html) (among many other packages).

---
deletable: true
editable: true
...

In this section, we're going to make use of three python third-party modules to apply PCoA and visualize the results 3D scatter plots. The data we'll use here is the full unweighted UniFrac distance matrix from a study of soil microbial communities across North and South America (originally published in [Lauber *et al.* (2009)](http://www.ncbi.nlm.nih.gov/pubmed/19502440)). We're going to use [pandas](http://pandas.pydata.org/) to manage the metadata, [scikit-bio](http://scikit-bio.org/) to manage the distance matrix and compute PCoA, and [matplotlib](http://matplotlib.org/) to visualize the results.

---
deletable: true
editable: true
...

First, we'll load sample metadata into a [pandas DataFrame](http://pandas.pydata.org/pandas-docs/dev/generated/pandas.DataFrame.html). These are really useful for loading and working with the type of tabular information that you'd typically store in a spreadsheet or database table. (Note that one thing I'm doing in the following cell is tricking pandas into thinking that it's getting a file as input, even though I have the information represented as tab-separated lines in a multiline string. [python's StringIO](https://docs.python.org/2/library/stringio.html) is very useful for this, and it's especially convenient in your unit tests... which you're writing for all of your code, right?) Here we'll load the tab-separated text, and then print it.

---
deletable: true
editable: true
...

```python
>>> from iab.data import lauber_soil_sample_md
>>> lauber_soil_sample_md
              pH                                         ENVO biome   Latitude
CF3.141691  3.56    ENVO:Temperate broadleaf and mixed forest biome  42.116667
PE5.141692  3.57                        ENVO:Tropical humid forests -12.633333
BF2.141708  3.61    ENVO:Temperate broadleaf and mixed forest biome  41.583333
CF2.141679  3.63    ENVO:Temperate broadleaf and mixed forest biome  41.933333
CF1.141675  3.92    ENVO:Temperate broadleaf and mixed forest biome  42.158333
HF2.141686  3.98    ENVO:Temperate broadleaf and mixed forest biome  42.500000
BF1.141647  4.05    ENVO:Temperate broadleaf and mixed forest biome  41.583333
PE4.141683  4.10                        ENVO:Tropical humid forests -13.083333
PE2.141725  4.11                        ENVO:Tropical humid forests -13.083333
PE1.141715  4.12                        ENVO:Tropical humid forests -13.083333
PE6.141700  4.12                        ENVO:Tropical humid forests -12.650000
TL3.141709  4.23                                     ENVO:shrubland  68.633333
HF1.141663  4.25    ENVO:Temperate broadleaf and mixed forest biome  42.500000
PE3.141731  4.25                        ENVO:Tropical humid forests -13.083333
BB1.141690  4.30                        ENVO:Tropical humid forests  44.870000
MP2.141695  4.38    ENVO:Temperate broadleaf and mixed forest biome  49.466667
MP1.141661  4.56                          ENVO:Temperate grasslands  49.466667
TL1.141653  4.58                                     ENVO:grassland  68.633333
BB2.141659  4.60    ENVO:Temperate broadleaf and mixed forest biome  44.866667
LQ3.141712  4.67                        ENVO:Tropical humid forests  18.300000
CL3.141664  4.89    ENVO:Temperate broadleaf and mixed forest biome  34.616667
LQ1.141701  4.89                        ENVO:Tropical humid forests  18.300000
HI4.141735  4.92  ENVO:Tropical and subtropical grasslands, sava...  20.083333
SN1.141681  4.95    ENVO:Temperate broadleaf and mixed forest biome  36.450000
LQ2.141729  5.03                        ENVO:Tropical humid forests  18.300000
CL4.141667  5.03                          ENVO:Temperate grasslands  34.616667
DF3.141696  5.05    ENVO:Temperate broadleaf and mixed forest biome  35.966667
BZ1.141724  5.12                                        ENVO:forest  64.800000
SP2.141678  5.13                          ENVO:Temperate grasslands  36.616667
IE1.141648  5.27                          ENVO:Temperate grasslands  41.800000
...          ...                                                ...        ...
DF2.141726  6.84    ENVO:Temperate broadleaf and mixed forest biome  35.966667
GB1.141665  6.84                                     ENVO:grassland  39.333333
SR1.141680  6.84                                     ENVO:shrubland  34.700000
SA1.141670  6.90                                        ENVO:forest  35.366667
SR3.141674  6.95                                     ENVO:shrubland  34.683333
KP4.141733  7.10                                     ENVO:shrubland  39.100000
GB3.141652  7.18                                        ENVO:forest  39.316667
CA1.141704  7.27                                        ENVO:forest  36.050000
BP1.141702  7.53                                     ENVO:grassland  43.750000
GB2.141732  7.57                                        ENVO:forest  39.316667
MT1.141719  7.57                                        ENVO:forest  46.800000
JT1.141699  7.60                                     ENVO:shrubland  33.966667
MD2.141689  7.65                                     ENVO:shrubland  34.900000
SF1.141728  7.71                                     ENVO:shrubland  35.383333
CM1.141723  7.85                          ENVO:Temperate grasslands  33.300000
MD3.141707  7.90                                     ENVO:shrubland  34.900000
KP3.141658  7.92                                     ENVO:shrubland  39.100000
SB1.141730  7.92                                     ENVO:shrubland  34.466667
RT1.141654  7.92                                     ENVO:shrubland  31.466667
SR2.141673  8.00                                     ENVO:shrubland  34.683333
CR1.141682  8.00                          ENVO:Temperate grasslands  33.933333
CA2.141685  8.02                                     ENVO:shrubland  36.050000
RT2.141710  8.07                          ENVO:Temperate grasslands  31.466667
MD5.141688  8.07                                     ENVO:shrubland  35.200000
SA2.141687  8.10                                     ENVO:shrubland  35.366667
GB5.141668  8.22                                     ENVO:shrubland  39.350000
SV1.141649  8.31                                     ENVO:shrubland  34.333333
SF2.141677  8.38                                     ENVO:shrubland  35.383333
SV2.141666  8.44                                     ENVO:grassland  34.333333
MD4.141660  8.86                                     ENVO:shrubland  35.200000

[89 rows x 3 columns]
```

---
deletable: true
editable: true
...

Just as one simple example of the many things that pandas can do, to look up a value, such as the pH of sample ``MT2.141698``, we can do the following. If you're interesting in learning more about pandas, [*Python for Data Analysis*](http://shop.oreilly.com/product/0636920023784.do) is a very good resource.

---
deletable: true
editable: true
...

```python
>>> lauber_soil_sample_md['pH']['MT2.141698']
6.6600000000000001
```

---
deletable: true
editable: true
...

Next we'll load our distance matrix. This is similar to ``human_microbiome_dm_data`` one that we loaded above, just a little bigger. After loading, we can visualize the resulting ``DistanceMatrix`` object for a summary.

---
deletable: true
editable: true
...

```python
>>> from iab.data import lauber_soil_unweighted_unifrac_dm
>>> _ = lauber_soil_unweighted_unifrac_dm.plot(cmap='Greens')
```

---
deletable: true
editable: true
...

Does this visualization help you to interpret the results? Probably not. Generally we'll need to apply some approaches that will help us with interpretation. Let's use ordination here. We'll run Principal Coordinates Analysis on our ``DistanceMatrix`` object. This gives us a matrix of coordinate values for each sample, which we can then plot. We can use ``scikit-bio``'s implementation of PCoA as follows:

---
deletable: true
editable: true
...

```python
>>> from skbio.stats.ordination import pcoa
...
>>> lauber_soil_unweighted_unifrac_pc = pcoa(lauber_soil_unweighted_unifrac_dm)
```

---
deletable: true
editable: true
...

What does the following ordination plot tell you about the relationship between the similarity of microbial communities taken from similar and dissimilar latitudes?

---
deletable: true
editable: true
...

```python
>>> _ = lauber_soil_unweighted_unifrac_pc.plot(lauber_soil_sample_md, 'Latitude', cmap='Greens', title="Samples colored by Latitude", axis_labels=('PC1', 'PC2', 'PC3'))
```

---
deletable: true
editable: true
...

If the answer to the above question is that there doesn't seem to be much association, you're on the right track. We can quantify this, for example, by testing for correlation between pH and value on PC 1.

---
deletable: true
editable: true
...

```python
>>> from scipy.stats import spearmanr
>>> spearman_rho, spearman_p = spearmanr(lauber_soil_unweighted_unifrac_pc.samples['PC1'],
...                                      lauber_soil_sample_md['Latitude'][lauber_soil_unweighted_unifrac_pc.samples.index])
>>> print('rho: %1.3f' % spearman_rho)
>>> print('p-value: %1.1e' % spearman_p)
rho: 0.158
p-value: 1.4e-01
```

---
deletable: true
editable: true
...

In the next plot, we'll color the points by the pH of the soil sample they represent. What does this plot suggest about the relationship between the similarity of microbial communities taken from similar and dissimilar pH?

---
deletable: true
editable: true
...

```python
>>> _ = lauber_soil_unweighted_unifrac_pc.plot(lauber_soil_sample_md, 'pH', cmap='Greens', title="Samples colored by pH", axis_labels=('PC1', 'PC2', 'PC3'))
```

---
deletable: true
editable: true
...

```python
>>> from scipy.stats import spearmanr
>>> spearman_rho, spearman_p = spearmanr(lauber_soil_unweighted_unifrac_pc.samples['PC1'],
...                                      lauber_soil_sample_md['pH'][lauber_soil_unweighted_unifrac_pc.samples.index])
>>> print('rho: %1.3f' % spearman_rho)
>>> print('p-value: %1.1e' % spearman_p)
rho: -0.958
p-value: 1.9e-48
```

---
deletable: true
editable: true
...

Taken together, these plots and statistics suggest that soil microbial community composition is much more closely associated with pH than it is with latitude: the key result that was presented in [Lauber *et al.* (2009)](http://www.ncbi.nlm.nih.gov/pubmed/19502440).

---
deletable: true
editable: true
...

## PCoA versus PCA: what's the difference? <link src='163769'/>

---
deletable: true
editable: true
...

You may have also heard of a method related to PCoA, called Principal Components Analysis or PCA. There is, however, an important key difference between these methods. PCoA, which is what we've been working with, performs ordination with a distance matrix as input. PCA on the other hand performs ordination with sample by feature frequency data, such as the OTU tables that we've been working with, as input. It achieves this by computing Euclidean distance (see [here](http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.euclidean.html#scipy.spatial.distance.euclidean)) between the samples and then running PCoA. So, if your distance metric is Euclidean, PCA and PCoA are the same. In practice however, we want to be able to use distance metrics that work better for studying biological diversity, such as Bray-Curtis or UniFrac. Therefore we typically compute distances with whatever metric we want, and then run PCoA.

---
deletable: true
editable: true
...

## Are two different analysis approaches giving me the same result? <link src='371f0d'/>

---
deletable: true
editable: true
...

A question that comes up frequently, often in method comparison, is whether two different approaches for analyzing some data giving the consistent results. This could come up, for example, if you were comparing DNA sequence data from the same samples generated on the 454 Titanium platform with data generated on the Illumina MiSeq platform to see if you would derive the same biological conclusions based on either platform. This was done, for example, in [Additional Figure 1](http://genomebiology.com/2011/12/5/R50/additional) of [*Moving Pictures of the Human Microbiome*](http://genomebiology.com/content/12/5/R50). Similarly, you might wonder if two different OTU clustering methods or beta diversity metrics would lead you to the same biological conclusion. Let's look at one way that you might address this question.

---
deletable: true
editable: true
...

Imagine you ran three different beta diversity metrics on your feature table: unweighted UniFrac, Bray-Curtis, and weighted UniFrac (the quantitative analog of unweighted UniFrac), and then generated the following PCoA plots.

---
deletable: true
editable: true
...

```python
>>> _ = lauber_soil_unweighted_unifrac_pc.plot(lauber_soil_sample_md, 'pH', cmap='Greens',
...                                                title="Unweighted UniFrac, samples colored by pH",
...                                                axis_labels=('PC1', 'PC2', 'PC3'))
```

---
deletable: true
editable: true
...

```python
>>> from iab.data import lauber_soil_bray_curtis_dm
...
>>> lauber_soil_bray_curtis_pcoa = pcoa(lauber_soil_bray_curtis_dm)
...
>>> _ = lauber_soil_bray_curtis_pcoa.plot(lauber_soil_sample_md, 'pH', cmap='Greens',
...                                         title="Bray-Curtis, samples colored by pH",
...                                         axis_labels=('PC1', 'PC2', 'PC3'))
```

---
deletable: true
editable: true
...

```python
>>> from iab.data import lauber_soil_weighted_unifrac_dm
...
>>> lauber_soil_weighted_unifrac_pcoa = pcoa(lauber_soil_weighted_unifrac_dm)
...
>>> _ = lauber_soil_weighted_unifrac_pcoa.plot(lauber_soil_sample_md, 'pH', cmap='Greens',
...                                              title="Weighted UniFrac, samples colored by pH",
...                                              axis_labels=('PC1', 'PC2', 'PC3'))
/Users/caporaso/miniconda3/envs/iab/lib/python3.4/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:102: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.010291669756329344 and the largest is 3.8374200744108204.
  RuntimeWarning
```

---
deletable: true
editable: true
...

Specifically, what we want to ask when comparing these results is **given a pair of ordination plots, is their shape (in two or three dimensions) the same?** The reason we care is that we want to know, **given a pair of ordination plots, would we derive the same biological conclusions regardless of which plot we look at?**

---
deletable: true
editable: true
...

We can use a [Mantel test](http://scikit-bio.org/docs/latest/generated/generated/skbio.stats.distance.mantel.html) for this, which is a way of testing for correlation between distance matrices.

---
deletable: true
editable: true
...

```python
>>> from skbio.stats.distance import mantel
...
>>> r, p, n = mantel(lauber_soil_unweighted_unifrac_dm, lauber_soil_weighted_unifrac_dm, method='spearman', strict=False)
>>> print("Mantel r: %1.3f" % r)
>>> print("p-value: %1.1e" % p)
>>> print("Number of samples compared: %d" % n)
Mantel r: 0.906
p-value: 1.0e-03
Number of samples compared: 88
```

---
deletable: true
editable: true
...

```python
>>> r, p, n = mantel(lauber_soil_unweighted_unifrac_dm, lauber_soil_bray_curtis_dm, method='spearman', strict=False)
>>> print("Mantel r: %1.3f" % r)
>>> print("p-value: %1.1e" % p)
>>> print("Number of samples compared: %d" % n)
Mantel r: 0.930
p-value: 1.0e-03
Number of samples compared: 88
```

---
deletable: true
editable: true
...

```python
>>> r, p, n = mantel(lauber_soil_weighted_unifrac_dm, lauber_soil_bray_curtis_dm, method='spearman', strict=False)
>>> print("Mantel r: %1.3f" % r)
>>> print("p-value: %1.1e" % p)
>>> print("Number of samples compared: %d" % n)
Mantel r: 0.850
p-value: 1.0e-03
Number of samples compared: 88
```

---
deletable: true
editable: true
...

The way that we'd interpret these results is that, although the plots above look somewhat different from one another, the underlying data (the distances between samples) are highly correlated across the different diversity metrics. As a result, we'd conclude that with any of these three diversity metrics we'd come to the conclusion that samples that are more similar in pH are more similar in their microbial community composition.

---
deletable: true
editable: true
...

We could apply this same approach, for example, if we had clustered sequences into OTUs with two different approaches. For example, if we used *de novo* OTU picking and open reference OTU picking, we could compute UniFrac distance matrices based on each resulting feature table, and then compare those distance matrices with a Mantel test. This approach was applied in [Rideout et al 2014](https://peerj.com/articles/545/) to determine which OTU clustering methods would result in different biological conclusions being drawn from a data set.

---
deletable: true
editable: true
...

### Procrustes analysis <link src='baaa8e'/>

---
deletable: true
editable: true
...

A related approach, but which I think is less useful as it compares PCoA plots directly (and therefore a summary of the distance data, rather than the distance data itself) is called Procrustes analysis (you can read about the origin of the name [here](http://en.wikipedia.org/wiki/Procrustes)). Procrustes analysis takes two coordinate matrices as input and effectively tries to find the best superimposition of one on top of the other. The transformations that are applied are as follows:

---
deletable: true
editable: true
...

* Translation (the mean of all points is set to 1 on each dimension)
* Scaling (root mean square distance of all points from the origin is 1 on each dimension)
* Rotation (choosing one set of points as the reference, and rotate the other to minimize the sum of squares distance (SSD) between the corresponding points)

---
deletable: true
editable: true
...

The output is a pair of *transformed coordinate matrices*, and an $M^{2}$ statistic which represents how dissimilar the coordinate matrices are to each other (so a small $M^{2}$ means that the coordinate matrices, and the plots, are more similar). [Procrustes analysis is implemented in scikit-bio](http://scikit-bio.org/generated/skbio.maths.stats.spatial.procrustes.html).

---
deletable: true
editable: true
...

## Where to go from here <link src='fc527a'/>

---
deletable: true
editable: true
...

If you're interested in learning more about the topics presented in this chapter, I recommend [Measuring Biological Diversity](http://www.amazon.com/Measuring-Biological-Diversity-Anne-Magurran/dp/0632056339) by Anne E. Magurran, and the [QIIME tutorials](http://qiime.org/tutorials/index.html). The [QIIME software package](http://www.qiime.org) is designed for performing the types of analyses described in this chapter.

---
deletable: true
editable: true
...

## Acknowledgements <link src='0d9866'/>

---
deletable: true
editable: true
...

Much of content in this section is based on knowledge that I gained through years of working with [Rob Knight](https://knightlab.ucsd.edu/) and the rest of the [QIIME](http://qiime.org/) [development group](https://github.com/biocore/qiime/graphs/contributors). Thanks everyone, I'm looking forward to many more years of productive, fun and exciting work together!
