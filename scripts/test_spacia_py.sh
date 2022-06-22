#!/bin/bash

# Need p3s
# Make sure your python have scipy, scikitlean, and pandas.
module purge
module load shared
module load R/4.0.2-gccmkl

python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/spacia.py \
    /project/shared/xiao_wang/projects/cell2cell_inter/data/Counts.txt \
    /project/shared/xiao_wang/projects/cell2cell_inter/data/Spot_metadata.txt \
    -rf FGFR1,FGFR2 \
    -o /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spaciapy_test \
    -rc C_1 -sc C_2 \
    -sf FGF1,FGF2,FGF7,FGF9,FGF10,FGF11,FGF12,FGF13,FGF14,FGF17,FGF18

# python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/spacia.py \
#     /project/shared/xiao_wang/projects/cell2cell_inter/data/Counts.txt \
#     /project/shared/xiao_wang/projects/cell2cell_inter/data/Spot_metadata.txt \
#     -rf FGFR1,FGFR2 \
#     -o /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spaciapy_test_cf \
#     -cf /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spacia_test_input/input_cells.csv \
#     -sf FGF1,FGF2,FGF7,FGF9,FGF10,FGF11,FGF12,FGF13,FGF14,FGF17,FGF18

# python /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/spacia.py \
#     /project/shared/xiao_wang/projects/cell2cell_inter/data/Counts.txt \
#     /project/shared/xiao_wang/projects/cell2cell_inter/data/Spot_metadata.txt \
#     -rf /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spacia_test_input/example_receiver_pathways.csv \
#     -o /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spaciapy_test_cf_rf \
#     -cf /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spacia_test_input/input_cells.csv \
#     -sf FGF1,FGF2,FGF7,FGF9,FGF10,FGF11,FGF12,FGF13,FGF14,FGF17,FGF18
