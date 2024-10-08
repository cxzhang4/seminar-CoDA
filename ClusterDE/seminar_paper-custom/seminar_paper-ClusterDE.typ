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
  title: [ClusterDE: a post-clustering differential expression method],
  authors: (
    ( name: [Carson Zhang],
      affiliation: [LMU Munich],
      email: [carson.zhang\@campus.lmu.de] ),
    ),
  abstract: [In typical differential expression analysis, a clustering algorithm is applied to scRNA-seq data, and then a differential expression test is conducted in order to identify genes that are differentially expressed between the clusters. However, this procedure constitutes "double dipping", as it first clusters the data to identify cell types, and then uses those same clusters to identify cell-type marker genes. This leads to an inflated FDR for DE genes. #cite(<ClusterDE>) propose ClusterDE, a post-clustering DE method that controls the FDR of DE genes. ClusterDE generates a synthetic null dataset that preserves the structure of the real data, computes differences between this null dataset and the real data, then performs FDR control on the results. Simulations and real data analysis demonstrate that ClusterDE controls the FDR and identifies cell-type marker genes as top DE genes, successfully distinguishing them from housekeeping genes.

],
  font: ("arial",),
  fontsize: 11pt,
  toc: true,
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)


== Introduction
<introduction>
=== Cell-type annotation
<cell-type-annotation>
==== Motivation
<motivation>
Understanding which types of cells are in a data sample allows an analyst to better make use of existing knowledge about those cells. "Cell annotation" is the process of labeling cells in a sample of data. In this paper, the focus is on annotating the "cell type" of each cell: a cellular phenotype that is robust across datasets #cite(<HeumosSchaarLance2023>);. For example, plasma B cells are one type of white blood cell that are involved in the human body’s immune response by secreting antibodies #cite(<HeumosSchaarLance2023>);. T cells are another type of white blood cell that are also involved in the immune response #cite(<GlossaryLymphocyte>);. They produce cytokines, which are signaling proteins that activate other parts of the human immune system. A scientist interested in a patient’s immune response may be interested in the counts of B cells and T cells \(and their subtypes): for example, in order to better understand the roles of each cell, or how they affect patient outcomes. Cell-type annotation is required in order to obtain this information from e.g.~a blood sample.

==== Cell-type markers
<cell-type-markers>
==== \(TODO: other methods of annotation)
<todo-other-methods-of-annotation>
=== Differential expression testing
<differential-expression-testing>
Differential expression testing is the primary method by which scientists identify marker genes. If genes are

\(define validity)

\(define FDR)

\(mention FDR control like Benjamini-Hochberg)

TODO: mention Scanpy, Seurat, and their default methods

To identify these cell types, they identify a set of cell-type marker genes.

To identify these cell-type marker genes, they perform differential expression testing.

Naive differential expression testing is susceptible to false discoveries caused by double-dipping.

Biologists like to identify the cell types in their scRNA-seq samples.

To identify these cell types, they identify a set of cell-type marker genes.

To identify these cell-type marker genes, they perform differential expression testing.

Naive differential expression testing is susceptible to false discoveries caused by double-dipping.Biologists like to identify the cell types in their scRNA-seq samples.

To identify these cell types, they identify a set of cell-type marker genes.

To identify these cell-type marker genes, they perform differential expression testing.

Naive differential expression testing is susceptible to false discoveries caused by double-dipping.

== The double-dipping issue
<the-double-dipping-issue>
== Differential expression methods that address the issue
<differential-expression-methods-that-address-the-issue>
=== Count splitting
<count-splitting>
=== TN Test
<tn-test>
== ClusterDE
<clusterde>
=== Summary of steps
<summary-of-steps>
The ClusterDE method consists of four basic steps.

+ Generate a synthetic null dataset that consists of a single cluster but otherwise mimics the real data.

+ Separately for each dataset, cluster the cells into two groups.

+ Separately for each dataset, perform differential expression testing between the two groups from step 2.

+ Combine the results to determine which genes to output as discoveries \(DE genes). #cite(<RafiGreenland2020>)

TODO: similarity to the S-value.

#figure([
#box(width: 1304.9096989966556pt, image("figures/ClusterDE_illustration.png"))
], caption: figure.caption(
position: bottom, 
[
A visual overview of the ClusterDE method.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
numbering: "1", 
)
<fig-ClusterDE-illustration>


==== Step 1
<step-1>
===== Idea: negative control
<idea-negative-control>
==== Step 2
<step-2>
==== Step 3: choice of tests
<step-3-choice-of-tests>
==== Step 4: false discovery rate control using Clipper
<step-4-false-discovery-rate-control-using-clipper>
In Step 4, ClusterDE uses the Clipper method to choose which discoveries from step 3 to output as true DE genes.

===== Intuition
<intuition>
Given that the negative control generated in step 1 accomplished its goal, the two datasets should be similar, and therefore the p-values \(and DE scores) outputted by each test should be similar. This means that, when a test on a given gene has a very low p-value, but this p-value is similar across both datasets, it is reasonable to believe that this low p-value occurred due to noise. However, when a p-value is much lower in the real data than in the synthetic null data, this indicates that the gene is truly differentially expressed between the two clusters.

== Practical notes on ClusterDE usage
<practical-notes-on-clusterde-usage>
- Only 2 clusters

#block[
```r
# library(flextable)
# library(tinytable)
# library(readr)
```

]
#block[
```r
# table_e1 = read_csv("../seminar_paper-bacsc/python/reproduce_results/table_e1.csv")
# flextable(table_e1)
# 
# table_e1_synthetic = read_csv("../seminar_paper-bacsc/python/synthetic_null_generation/table_e1-synthetic_data.csv")
# flextable(table_e1_synthetic)
```

]
== Data analysis
<data-analysis>
=== BacSC data
<bacsc-data>
=== Synthetic null data generation
<synthetic-null-data-generation>
=== Schäfer-Strimmer
<schäfer-strimmer>
=== Results
<results>
== Simulation study
<simulation-study>
#cite(<HeumosSchaarLance2023>) #cite(<BenjaminiHochberg1995>)

== Appendix
<appendix>



#set bibliography(style: "apa")

#bibliography("refs.bib")

