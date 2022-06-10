#!/bin/bash
#SBATCH --job-name=spacia_simulation                                           # job name
#SBATCH --partition=256GB                                         # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                                    # number of nodes requested by user
#SBATCH --cpus-per-task=48                                                  # number of total tasks
#SBATCH --time=2-00:00:00                                            # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=./spacia_simulation_%j                                 # redirect both standard output and erro output to the same file
#SBATCH --error=./spacia_simulation_%j

source /home2/s190548/.bashrc
module purge
module load shared
module load R/4.0.2-gccmkl
module load python/3.6.1-2-anaconda
conda activate p3s

# python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py 50 5
# no noise
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 3 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/ &
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 4 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/ &
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 5 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/ &
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 6 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/ &
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 7 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 8 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/ &
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 9 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/ &
python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/simulation.py -n 10 -o /project/shared/xiao_wang/projects/cell2cell_inter/data/shijia_simulation/
wait 
