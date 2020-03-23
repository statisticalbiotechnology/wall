# Disposition

- Data folder
	- Contains most of the data needed for the code, some files are too big to store on github
- Src folder
	- Contains the code base, mostly with ipynb files for developing scripts but also py files for already completed scripts
	- The sunburst folder contains all the files needed to view the sunburst which can be run through pythons http.server module

```bash
#installing the module
python3 -m pip install http.server
#running the server
python3 -m http.server
#running the server with custom port e.g. 8001
python -m http.server 8001
#Accessing the server from the browser (assuming you are in the sunburst folder)
localhost:8001/sunburst.html
```
- Exp folder
	- Contains the experiments. The first version contains the csv files for the intclust qvalues and the receptor qvalues
- Porch
	- Porch is included especially for the use of the q-value function. As of the first version of wall, porch is not used ext

# Results
Results are derived from statistical tests comparing each cluster to all the other clusters in the context of pathway activity


## What is known about the clusters from the metabric paper
| Cluster | p53 Mutation Frequency| Driver Genes | MATH score |Outcome|
| ------------- |:-------------:| :---:|:---:|:---:|
| 1|29%| GATA3, fewer alterations in Akt|-|4|
| 2|24.1%|CCND1, PAK1 |Low|4|
| 3|10.0% |Clonal PIK3CA, clonal inactivating MAP3K1, Subclonal inactivating MAP3K |-|1|
| 4ER+|21.1% |Clonal PIK3CA |-|1.5|
| 4ER-|50.5% | -|Low|4|
| 5|64.2%| ERBB2 (HER2)|-|8|
| 6|40.7% | ZNF703|-|3|
| 7|14.0% | MAPK|Low|1|
| 8|4.4% | GATA3|Low|2|
| 9|44.7% | DNA Damage Response, Activating PIK3CA|-|4|
| 10|84.6% | DNA Damage Response, Cell Cycle regulation, Ubiquitination|Highest|6|





## What the 30 most differentially activated pathways tell us

| Cluster|Pathways|
| ------- |:----:|
| 1|Epigenetic changes, senescence|
| 2|Insulin receptor, TRAIL Signaling|
| 3|Mitosis/cell cycle, alpha-linoleic(omega-3/omega-6) degradation|
| 4ER+|HDR, D loops|
| 4ER-|Akt Signaling, Cell-Cell Communication|
| 5|Hormones, especially amines (epinephrine/melatonin)|
| 6|FGFR, Wnt, IP3|
| 7|NMDA, Ca2+, ion channels, CaM|
| 8|Noncanonical activation of NOTCH3|
| 9|Retroviral genome, riboflavin metabolism |
| 10| Cell Cycle, CDK|


# Sunburst plots

## Sunburst Plots Of Intclust vs Other Intclusts
| Cluster | Pathways |
|:---:|:---:|
| 1| Cell Cycle, Base Excision Repair, DNA Double Strand Break Repair,  RNA Pol I transcription, Gene Silencing by RNA, Senescence, Rho GTPase Signaling, ESR Signaling, TCF complex|
| 2| mTORC1, Iron Uptake and Transport, Insulin Receptor|
| 3| Cell Cycle, DNA Replication, DNA Double Strand Break Repair, Mismatch Repair, DNA Damage Bypass, SUMOylation, Megakaryocyte Production|
| 4ER+| Cell Cycle, DNA Double Strand Break Repair, Rho GTPase Signaling|
| 4ER-| SLC Mediated Transport, Immune System, Neuronal System, Metabolism of Vitamins and Cofactors, DAG/IP3 Signaling, CaM Pathway|
| 5| Cell Cycle, Metabolism of Amine Derived Hormones, Tryptophan Catabolism, Rho gtpases|
| 6| IP Metabolism, Wnt|
| 7| CaM Pathway, NMDA Receptors, Cilium Assembly, Metabolism of Vitamins/Cofactors, Metabolism of Nucleotides, Ion Channel Transport|
| 8| Immune System, Transport of Small Molecules, Sphingolipid Metabolism, GPCR Signaling, DAG/IP3, Negative Regulation of Akt Network, Hemostasis|
| 9| Cell Cycle, KSHRP/TTP/TRF1, Riboflavin Metabolism, Glutamine/Glutamate Metabolism|
| 10| Cell Cycle, DNA Replication, Rho GTPases, Sphingolipid Metabolism, Nucleotide Metabolism|

## Sunburst Plots Of Intclust vs Benign Tissue
| Cluster | Pathways |
|:---:|:---:|
| 1| Cell Cycle, Signal Transduction, DNA repair, transport of small molecules, cell-cell communications, DNA replication, muscle contraction, neuronal system|
| 2| Cell Cycle, Signal Transduction, DNA repair, transport of small molecules, cell-cell communications, DNA replication, muscle contraction, neuronal system, Hemostasis|
| 3| Signal transduction, cell cycle, cell-cell communication, neuronal system|
| 4ER-| Cell Cycle, signal transduction, DNA replication|
| 4ER+| Cell Cycle, signal transduction, DNA replication, transport of small molecules, neuronal system, cell-cell communication|
| 5| Cell cycle, transport of small molecules, signal transduction, DNA replication |
| 6| Cell Cycle, Signal transduction, transport of small molecules, hemostasis, DNA repair, cell-cell communication, dna replication, neuronal system, muscle contraction |
| 7| Signal trnasduction, transport of small molecules, neuronal system, cell cycle, cell-cell communication, hemostasis, muscle contraction |
| 8| Signal transduction, hemostasis, transport of small molecules, muscle contraction, cell-cell communication, neuronal system|
| 9| Cell Cycle, Signal Transduction, DNA repair, transport of small molecules, cell-cell communications, DNA replication, muscle contraction, neuronal system|
| 10| Cell Cycle, Signal Transduction, DNA repair, transport of small molecules, cell-cell communications, DNA replication, muscle contraction, neuronal system, immune system|





## Sunburst Plots For Receptors/Receptor Signatures

| Cluster | Pathways |
|:---:|:---:|
| ER Receptor| Immune System, Cell Cycle, Neuronal System, Signal transduction, SLC mediated transport, Keratinization|
| PR Receptor| Immune System, Cell Cycle, NMDA receptor, G protein mediated receptor, DAG, IP3, fatty acid metabolism|
| HER2 Receptor| Cell Cycle, visual phototransduction, metabolism of amine derived hormones, tryptophan metabolism|
| Triple Negative signature| Neuronal system, small molecule transport, lipid metabolism, IP metabolism, metabolism of cofactors, G protein modulated events, CaM pathway|
| ER-/PR-/HER2+| metabolism of amine derived hormones, tryptophan metabolism, Cell Cycle, G1 phase|



## What the 30 most differently activated GSEA pathways tell us

| Cluster | Pathways |
|:----:|:----:|
| 1 | Cell Cycle, Chromosome Maintenance, Mitochondrial Translation |
| 2 | Collagen Formation, Chromosome Maintenance, Telomeres, Cell Cycle |
| 3 | Cell Cycle, Chromosome Stability and Maintenance|
| 4ER+ | Cell Cycle, Mitochondrial Translation, HDR, SUMOylation|
| 4ER- | Mitochondrial Translation, Signal Transduction, Immune System, Intraflaggelar Transport |
| 5 | Cell Cycle, Immune System, Metabolism (of Proteins) |
| 6 | Cell Cycle, Complements, Chemokines, PD-1 signaling(Nobel 2018), Immunoregulatory Functions |
| 7 | Cell Cycle, Immunoregulatory Between Lymphoid/Non-Lymphoid Cell |
| 8 | Immune System, Parts of Cell Cycle |
| 9 | Signal Transduction |
| 10 | Cell Cycle, Metabolism of ncRNA, DNA Replication |


Cluster 1, 4ER+  and 10 had more than 30 pathways with fdr == 0. 4ER+ had 72 for example





### Comparing wall and GSEA significant pathways


![Comparing wall and GSEA significant pathways](/src/gsea_wall_comparison.png)
