# Single cell pathway analysis
## Getting started

The entire pipeline can be executed via a shellscript. Depending on how much of the example you want to run, you can use part of the shell script.
However, to make the shell script executable, run:

```bash
#locate the src folder
cd src
#make the file executable
chmod +x single_cell_pipeline
#Run the file:
./single_cell_pipeline.sh

```
Due to difficulties of memory management in linux, the differential expression results are provided in /data/lfc_df.csv and the shell script contains an if statement that checks for this file. If it exists, it skips the differential expression, if it doesnt exist, differential expression will be run.
