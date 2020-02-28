# wall


## Disposition
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
	- contains the experiments. The first version contains the csv files for the intclust qvalues and the receptor qvalues
- Porch
	- porch is included especially for the use of the q-value function. As of the first version of wall, porch is not used ext

# Results
Results are derived from statistical tests comparing each cluster to all the other clusters in the context of pathway activity


## What is known about the clusters from the metabric paper
| cluster  | p53 mutation frequency| driver gene|MATH score|Outcome|
| ------------- |:-------------:| :---:|:---:|:---:|
| 1|29%| GATA3, fewer alterations in Akt|-||
| 2|24.1%|CCND1, PAK1 |Low||
| 3|10.0% |Clonal PIK3CA, clonal inactivating MAP3K1, Subclonal inactivating MAP3K |-||
|4ER+|21.1% |Clonal PIK3CA |-|| 
|4ER-|50.5% | -|Low||
| 5|64.2%| ERBB2 (HER2)|-||
| 6|40.7% | ZNF703|-||
| 7|14.0% | MAPK|Low||
| 8|4.4% | GATA3|Low||
| 9|44.7% | DNA Damage Response, Activating PIK3CA|-||
| 10|84.6% | DNA Damage Response, Cell Cycle regulation, Ubiquitination|Highest||





## What the most differentially activated pathways tell us

|Cluster|Pathways|
| ------- |:----:|
|1|Epigenetic changes, senescence|
|2|Insulin receptor, TRAIL Signaling|
|3|Mitosis/cell cycle, alpha-linoleic(omega-3/omega-6) degradation|
|4ER+|HDR, D loops|
|4ER-|Akt Signaling, Cell-Cell Communication|
|5|Hormones, especially amines (epinephrine/melatonin)|
|6|FGFR, Wnt, IP3|
|7|NMDA, Ca2+, ion channels, CaM|
|8|Noncanonical activation of NOTCH3|
|9|Retroviral genome, riboflavin metabolism |
|10| Cell Cycle, CDK|
