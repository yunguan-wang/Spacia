
suppressPackageStartupMessages(library(distances))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
########  read arguments  ########
parser <- OptionParser(
    description = "Spacia plotting master script.")
parser <- add_option(
    parser, c("-i", "--spacia_results"),
    help="Spacia result path containing results.")

parser <- add_option(
    parser, c("-p", "--spacia_path"), default='../spacia',type='character',
    help="Spacia path containing all individual dependant scripts")

parser <- add_option(
    parser, c("-j", "--job_id"), default=NULL,
    help="Spacia job id, denotes a receiving pathway.")


########  define input  ########
# args = parse_args(parser)
args = parse_args(
    parser, 
    c(
        '-i', '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spaciapy_test',
        '-j','FGFR1_correlated_genes')
    )
if (is.null(args$spacia_results)) {
    stop("Spacia results must be provided. See script usage (--help)")
}

results_path = args$spacia_results
spacia_path = args$spacia_path
receiving_pathway = args$job_id

model_input_path = file.path(results_path, 'model_input')
exp_sender_fn = file.path(model_input_path, 'exp_sender.json')
metadata_fn = file.path(model_input_path, 'metadata.csv')

exp_sender_fn = file.path(
    model_input_path, paste(receiving_pathway,'exp_receiver.csv', sep='_')
)
beta_fn = file.path(
    results_path, receiving_pathway, paste(receiving_pathway,'beta.txt', sep='_')
)

########  Load all required plotting functions  ########
source(file.path(spacia_path, "plot_density_ridgeline.R"))
source(file.path(spacia_path, "plot_interaction_circos.R"))
source(file.path(spacia_path, "plot_interaction_direction.R"))
source(file.path(spacia_path, "plot_interaction_proportion.R"))
source(file.path(spacia_path, "plot_interaction_river.R"))


setwd('/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts')

