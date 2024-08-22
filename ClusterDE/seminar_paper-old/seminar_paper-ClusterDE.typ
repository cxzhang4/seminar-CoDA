// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = [
  #line(start: (25%,0%), end: (75%,0%))
]

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): block.with(
    fill: luma(230), 
    width: 100%, 
    inset: 8pt, 
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.amount
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == "string" {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == "content" {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

#show figure: it => {
  if type(it.kind) != "string" {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    new_title_block +
    old_callout.body.children.at(1))
}

#show ref: it => locate(loc => {
  let target = query(it.target, loc).first()
  if it.at("supplement", default: none) == none {
    it
    return
  }

  let sup = it.supplement.text.matches(regex("^45127368-afa1-446a-820f-fc64c546b2c5%(.*)")).at(0, default: none)
  if sup != none {
    let parent_id = sup.captures.first()
    let parent_figure = query(label(parent_id), loc).first()
    let parent_location = parent_figure.location()

    let counters = numbering(
      parent_figure.at("numbering"), 
      ..parent_figure.at("counter").at(parent_location))
      
    let subcounter = numbering(
      target.at("numbering"),
      ..target.at("counter").at(target.location()))
    
    // NOTE there's a nonbreaking space in the block below
    link(target.location(), [#parent_figure.at("supplement") #counters#subcounter])
  } else {
    it
  }
})

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      block(
        inset: 1pt, 
        width: 100%, 
        block(fill: white, width: 100%, inset: 8pt, body)))
}



#let article(
  title: none,
  authors: none,
  date: none,
  abstract: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: (),
  fontsize: 11pt,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: "1",
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)

  if title != none {
    align(center)[#block(inset: 2em)[
      #text(weight: "bold", size: 1.5em)[#title]
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[Abstract] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}
#show: doc => article(
  title: [ClusterDE Seminar Paper],
  authors: (
    ( name: [Carson Zhang],
      affiliation: [],
      email: [] ),
    ),
  abstract: [In typical differential expression analysis, a clustering algorithm is applied to scRNA-seq data, and then a differential expression test is conducted in order to identify genes that are differentially expressed between the clusters. However, this procedure constitutes "double dipping", as it first clusters the data to identify cell types, and then uses those same clusters to identify cell-type marker genes. This leads to an inflated FDR for DE genes. Song et al.~\(2023) propose ClusterDE, a post-clustering DE method that controls the FDR of DE genes. ClusterDE generates a synthetic null dataset that preserves the structure of the real data, computes differences between this null dataset and the real data, then performs FDR control on the results. Simulations and real data analysis demonstrate that ClusterDE controls the FDR and identifies cell-type marker genes as top DE genes, successfully distinguishing them from housekeeping genes.

],
  font: ("Arial",),
  fontsize: 11pt,
  sectionnumbering: "1.1.a",
  toc: true,
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)


= Introduction
<introduction>
== Background
<background>
== Biology
<biology>
TODO: decide which of the following biology background is necessary to include in the presentation.

=== RNA
<rna>
RNA carries the genetic information specific in DNA. There are two main types of RNA.

#strong[Non-coding RNA] performs some biological function.

#strong[Messenger RNA] forms a template for protein production \(it codes for a protein which performs some biological function).

TODO: explain the jump from RNA to the UMI count matrix.

TODO: define UMI.

=== scRNA-seq
<scrna-seq>
TODO: give an explanation of scRNA-seq data collection and analysis.

== Double-dipping
<double-dipping>
Biologists are interested in determining the cell types of the cells in a sample. \(TODO: give a concrete example.)

A natural way to do this might look like this.

- Cluster the sample. We hope that each cluster represents a distinct cell type.
- For each gene, perform a hypothesis test comparing its expression levels between the two clusters.
- Applying a multiple testing correction procedure \(e.g.~Benjamini-Hochberg), determine which null hypotheses to reject: which genes you will say are differentially expressed.

TODO: explain the double-dipping problem in differential expression analysis in more detail.

However, this leads to double-dipping because we used the data twice.

+ First, we clustered the cells. The algorithm chose these two specific clusters because, in some sense, these two subsets of cells look the most different in terms of their gene expression levels. That is, the algorithm found variation in gene expression levels that drove this particular clustering.

+ Next, we test each gene to see how different the expression levels are between the two clusters. However, we already know that this variation is there: that’s how we got the clusters in the first place. So we are testing for variation that we already found through clustering, and therefore, we will find differentially expressed genes when we shouldn’t.

When we cluster a sample of a single cell type into two clusters, we cluster based on differences in gene expression that occurred due to random noise. This means we forced the genes to exhibit different distributions in the two clusters, so we are likely to find this forced variation instead of true variation in the dataset.

TODO: give a concrete example of double dipping.

= Method
<method>
== ClusterDE method: 4 steps
<clusterde-method-4-steps>
In words, the ClusterDE method can be broken up into four steps:

+ Generate a synthetic null dataset that mimics the structure \(in particular, the gene-gene correlation structure) of the original data.

+ Separately partition the synthetic null data and the target data \(real data) into two clusters.

+ Separately for the null and target data, perform hypothesis tests for differentially expressed genes between the two clusters. For each gene, compute some sort of difference between the scores on the two datasets.

+ Output a subset of the significant results from step 3 as potential cell-type marker genes.

== ClusterDE method
<clusterde-method>
#figure([
#box(width: 1304.9096989966556pt, image("ClusterDE_illustration.png"))
], caption: figure.caption(
position: bottom, 
[
A graphical illustration of ClusterDE.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


It is important to note that ClusterDE "does not provide an automatic decision about whether two clusters should be merged". Its outputs are potential DE genes, and therefore it does not directly measure the quality of a given clustering. These potential cell-type marker genes enable researchers to gain biological insights into the clusters, and they empower researchers to further explore the functional and molecular characteristics of the clusters.

== 1. Synthetic null generation
<synthetic-null-generation>
The synthetic null generation consists of three steps, as described in the following figure.

#figure([
#box(width: 1254.3411371237457pt, image("ClusterDE_supp_null_generation.png"))
], caption: figure.caption(
position: bottom, 
[
Null generation steps.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


==== 1. Model the null distribution in terms of the Gaussian copula.
<model-the-null-distribution-in-terms-of-the-gaussian-copula.>
TODO: #link("https://blog.carsonzhang.com/posts/probability_integral_transform/")[explain the probability integral transform];.

===== Sklar’s Theorem
<sklars-theorem>
#strong[Theorem] \(Sklar’s Theorem): Let $bold(X)$ be a $m$-dimensional random vector with joint cumulative distribution function $F$ and marginal distribution functions $F_j , j = 1 , . . . , m$. The joint CDF can be expressed as

$ F (x_1 , . . . , x_m) = C (F_1 (x_1) , . . . , F_m (x_m)) $

with associated probability density \(or mass) function

$ f (x_1 , . . . , x_m) = c (F_1 (x_1) , . . . , F_m (x_m)) f_1 (x_1) . . . f_m (x_m) $

for a $d$-dimensional copula $C$ with copula density $c$.

The inverse also holds: the copula corresponding to a multivariate CDF $F$ with marginal distribution functions $F_j , j = 1 , . . . , m$ can be expressed as

$ C (u_1 , . . . , u_m) = F (F_1^(- 1) (u_1) , . . . , F_m^(- 1) (u_m)) $ , and the copula density \(or mass) function is

$ c (u_1 , . . . , u_m) = frac(f (F_1^(- 1) (u_1) , . . . , F_m^(- 1) (u_m)), f_1 (F_1^(- 1) (u_1)) . . . f_m (F_m^(- 1) (u_m))) $ .

===== Why Sklar’s Theorem
<why-sklars-theorem>
Sklar’s Theorem allows us to choose $C$, the specific copula function that we will use to approximate $F$. For convenience, we choose $C$ to be a multivariate Gaussian distribution, since we know how to sample from it. Therefore, we have found a way to estimate the joint MVNB distribution of genes.

$ C (bold(u) ; bold(R)) = Phi_(bold(R)) (Phi^(- 1) (u_1) , . . . , Phi^(- 1) (u_m)) $

Now, our goal is to estimate the parameters $mu_j , sigma_j_(j = 1)^m$ and $bold(R)$.

==== 2. Fit the null model to the real data.
<fit-the-null-model-to-the-real-data.>
Recall that the power of the copula is that it allows us to consider the marginal distributions separately from the covariance structure. Therefore, we can proceed as follows.

For each gene $j$, estimate ${ mu_j , sigma_j }$ using maximum likelihood. These are the marginal distributions.

For the entire dataset, use the Gaussian copula to model the dependence structure.

- Transform the raw data \(counts) to the CDF values of the counts. Given the parameters for each marginal distribution that we just generated, compute the CDFs specified by those parameters.
- Transform the discrete CDF into a continuous $U (0 , 1)$ variable. Do this by computing $U_(i j) = V_(i j) hat(F)_j (Y_(i j)) + (1 - V_(i j)) hat(F)_j (Y_(i j))$.
- Transform the CDF values to standard Gaussian random variables. That is, compute $Phi^(- 1) (U_(i j))$ for each $j = 1 , . . , m$.
- Fit a $m$-dimensional multivariate Gaussian distribution to this data to compute $bold(hat(R))$. $bold(hat(R))$ is the sample correlation matrix of this data.
- Sample from this multivariate Gaussian distribution: $N_m (bold(0) , bold(hat(R)))$.

==== 3. Sample from the fitted null model.
<sample-from-the-fitted-null-model.>
- Generate a sample of size $n$ from $N_m (bold(0) , bold(hat(R)))$. \(Denote them as $tilde(Z)_(i j)$ because each is standard Gaussian.)

$ mat(delim: "[", tilde(Z)_11, dots.h, tilde(Z)_(1 m); dots.v, dots.down, ; tilde(Z)_(n 1), , tilde(Z)_(n m)) $

- Convert them to negative binomial count vectors. First we apply the marginal CDFs to transform them to $U (0 , 1)$ random variables. Then, we apply the inverse of each negative binomial CDF to transform them to the proper counts \(negative binomial random variables) for the corresponding gene.

$ mat(delim: "[", hat(F)_1^(- 1) (Phi (tilde(Z)_11)), dots.h, hat(F)_m^(- 1) (Phi (tilde(Z)_(1 m))); dots.v, dots.down, ; hat(F)_1^(- 1) (Phi (tilde(Z)_(n 1))), , hat(F)_m^(- 1) (Phi (tilde(Z)_(n m)))) $

== 2. Clustering
<clustering>
ClusterDE allows any clustering algorithm. Note that it only handles the case of two clusters, so if you started out with more clusters, you should identify a particular pair of interest. In the #strong[Practical guidelines for ClusterDE usage] subsection, steps 1 and 2 describe how an analyst should proceed.

+ Given $gt.eq 2$ clusters, identify 2 clusters of interest. Generally, this will be a pair for which you suspect the clustering is spurious \(i.e.~you think the two clusters actually come from the same cell type, so they are strong candidates to be merged into a single cluster).

+ Filter the data so that you only consider the subset of cells that come from those two clusters.

==== UMAP
<umap>
UMAP is common in scRNA-seq analysis.

TODO: summarize UMAP.

==== Louvain
<louvain>
The example analyses in the presentation use the default Seurat clustering procedure, which uses the Louvain algorithm.

TODO: describe the Seurat clustering pipeline.

TODO: summarize the Louvain algorithm.

=== 3. DE analysis \(testing)
<de-analysis-testing>
ClusterDE allows any DE test.

TODO: choose and summarize common DE tests.

Let $P_1 , . . . , P_m$ be the p-values computed by the $m$ DE tests on the target data. Define the target DE score $S_j := - log_10 P_j$. Likewise for the synthetic null data.

TODO: motivate the formula for the DE score. It is essentially the information content of the p-value.

The final outputs of step 3: $m$ target DE scores $S_1 , . . . , S_m$; $m$ null DE scores $tilde(S)_1 , . . . , tilde(S)_m$.

=== 4. FDR control
<fdr-control>
Given the target and null DE scores, compute a contrast score for gene $j$ as $C_j := S_j - tilde(S)_j$.

In reality, we will have some more positive contrast scores than negative ones.

In the case that there truly are two cell types, then we should definitely expect some true DE genes, and should therefore expect more positive contrast scores than negative ones.

Recall the contrast score formula: $C_j := S_j - tilde(S)_j$. When this is positive, then the DE test on the target data was more surprising under the null, which means it outputted a smaller p-value. This #emph[should] be the case if the null is false. This is why the authors say that, ideally, slightly less than 50% of all genes’ contrast scores should be negative: sometimes, the null is false, and when that is the case, we will see additional values in the right tail \(and only in the right tail) of the contrast score distribution.

#figure([
#box(width: 50%,image("ClusterDE_illustration-fdr_control.png"))
], caption: figure.caption(
position: bottom, 
[
FDR control.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


Motivation for the cutoff: under the null hypothesis, we have a symmetric distribution of contrast scores. Assume the negative tail is differences due to random noise. Under the null hypothesis, the right tail would be the same size.

Our real right tail will be bigger \(we hope), since there are true differentially expressed genes which we will discover. And each of these will appear in the right tail. So any instances where the null is violated will appear only in the right tail.

We choose the minimum $t$ to maximize our discoveries, i.e.~to "use" all of the false discoveries that we can.

== Results
<results>
=== Simulation
<simulation>
=== Real data example
<real-data-example>
== Appendix
<appendix>



