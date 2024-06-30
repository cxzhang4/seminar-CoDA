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
=== Biology
<biology>
TODO: give the basic biology background necessary to understand the paper.

==== RNA
<rna>
RNA carries the genetic information specific in DNA. There are two main types: - #strong[non-coding RNA] performs some biological function - #strong[messenger RNA] forms a template for protein production \(it codes for a protein which performs some biological function).

TODO: explain the jump from RNA to the UMI count matrix.

TODO: define UMI.

=== scRNA-seq
<scrna-seq>
TODO: give an explanation of scRNA-seq data collection and analysis.

=== Double-dipping
<double-dipping>
TODO: explain the double-dipping problem in differential expression analysis.

= Method
<method>
In words, the ClusterDE method can be broken up into four steps:

+ Generate a synthetic null dataset that mimics the structure \(in particular, the gene-gene correlation structure) of the original data.

+ Separately partition the synthetic null data and the target data \(real data) into two clusters.

+ Separately for the null and target data, perform hypothesis tests for differentially expressed genes between the two clusters. For each gene, compute some sort of difference between the scores on the two datasets.

+ Output a subset of the significant results from step 3 as potential cell-type marker genes.

#figure([
#box(width: 1304.9096989966556pt, image("ClusterDE_illustration.png"))
], caption: figure.caption(
position: bottom, 
[
An graphical illustration of ClusterDE.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


It is important to note that ClusterDE "does not provide an automatic decision about whether two clusters should be merged". Its outputs are potential DE genes, and therefore it does not directly measure the quality of a given clustering. These potential cell-type marker genes enable researchers to gain biological insights into the clusters, and they empower researchers to further explore the functional and molecular characteristics of the clusters.

== Synthetic null generation
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


+ Model the null distribution in terms of the Gaussian copula.

Now, our goal is to estimate the parameters $mu_j , sigma_j_(j = 1)^m$ and $bold(R)$.

#block[
#set enum(numbering: "1.", start: 2)
+ Fit the null model to the real data.

+ Sample from the fitted null model.
]

== 2. Clustering
<clustering>
ClusterDE allows any clustering algorithm. Note that it only handles the case of two clusters, so if you started out with more clusters, you should identify a particular pair of interest. In the #strong[Practical guidelines for ClusterDE usage] subsection, steps 1 and 2 describe how an analyst should proceed.

+ Given $gt.eq 2$ clusters, identify 2 clusters of interest. Generally, this will be a pair for which you suspect the clustering is spurious \(i.e.~you think the two clusters actually come from the same cell type, so they are strong candidates to be merged into a single cluster).

+ Filter the data so that you only consider the subset of cells that come from those two clusters.

TODO: describe the Seurat clustering pipeline.

=== UMAP
<umap>
UMAP is common.

TODO: summarize UMAP.

=== Louvain
<louvain>
The example analyses in the presentation use the default Seurat clustering procedure, which uses the Louvain algorithm.

TODO: summarize the Louvain algorithm.

== 3. DE analysis \(testing)
<de-analysis-testing>
ClusterDE allows any DE test.

TODO: choose and summarize common DE tests.

Let $P_1 , . . . , P_m$ be the p-values computed by the $m$ DE tests on the target data. Define the target DE score $S_j := - log_10 P_j$. Likewise for the synthetic null data.

The final outputs of step 3: $m$ target DE scores $S_1 , . . . , S_m$; $m$ null DE scores $tilde(S)_1 , . . . , tilde(S)_m$.

== 4. FDR control
<fdr-control>
Given the target and null DE scores, compute a contrast score for gene $j$ as $C_j := S_j - tilde(S)_j$.

= Results
<results>
== Simulation
<simulation>
== Real data example
<real-data-example>
= Appendix
<appendix>



