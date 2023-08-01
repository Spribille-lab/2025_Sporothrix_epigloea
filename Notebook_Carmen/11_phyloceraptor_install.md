##SET-UP Phyloceraptor in Narval

module load apptainer/1.1
module load python/3.7
virtualenv --no-download ~/ENV
source ~/ENV/bin/activate
pip install --no-index --upgrade pip
pip install snakemake
pip install pandas
git clone --recursive https://github.com/reslp/phylociraptor.git

cd ./phylociraptor
#-t serial=1 is to pass directly to the current processor in a single thread
./phylociraptor setup --verbose -t serial=1

if things fail then delete the logs but not the folder, delete the results and delete the snakemake logs

pull the .sif files



module load apptainer/1.1
module load python/3.7
source ~/ENV/bin/activate
cd $SCRATCH/phylociraptor/

# Setup the pipeline:
./phylociraptor setup -t serial=2
# Identify orthologous genes for all the genomes:
./phylociraptor orthology -t slurm -c data/cluster-config-SLURM.yaml --verbose
# Filter orthologs using according to settings in the config.yaml file:
./phylociraptor filter-orthology -t slurm -c data/cluster-config-SLURM.yaml --verbose
# Create alignments and trim them:
./phylociraptor align -t slurm -c data/cluster-config-SLURM.yaml --verbose
# Filter alignments according to settings in the config.yaml file:
./phylociraptor filter-align -t slurm -c data/cluster-config-SLURM.yaml --verbose
# Optionally you can run extensive model testing for individual alignments. This is done using iqtree. In case you run this step, the next step will use these models. Otherwise phylociraptor will use models specified in the config file.
./phylociraptor modeltest -t slurm -c data/cluster-config-SLURM.yaml --verbose
# Reconstruct phylogenies:
./phylociraptor njtree -t slurm -c data/cluster-config-SLURM.yaml --verbose
./phylociraptor mltree -t slurm -c data/cluster-config-SLURM.yaml --verbose
./phylociraptor speciestree -t slurm -c data/cluster-config-SLURM.yaml --verbose
# Create a report of the run:
# this wont work in cedar.  Use the script below.
./phylociraptor report --verbose --figure

# check status of the pipeline
./phylociraptor  check