echo "Welcome to the single cell pathway analysis/Sunburst plotting pipeline!";
echo "To start, we will download the datafiles and put them together in a dataframe we will call full_df.csv";
echo "There will also be a pheno_df file containing informartion about which cells are dress and which are hv";
echo "but first, let's install the required modules";
python3 -m pip install -r requirements.txt;
Rscript requirements.R;
python3 file_preparation.py;
FILE=../data/lfc_df.csv;
if [ -e "$FILE" ]; then
	echo "differential expression file already exist, proceeding with pathway analysis"
	python3 prepare_idea.py
	Rscript iDEA.R;
else
	echo "differential expresssion file does not exist, proceeding with differential expression"
	Rscript differential_expression.R;
	python3_prepare_idea.py
	Rscript iDEA.R
fi
echo "Lets gather the results and convert them to a csv containing the -log10 q-values";


python3 results.py;
echo "Lets make the sunburst plots";

python3 run_sunburst.py;
