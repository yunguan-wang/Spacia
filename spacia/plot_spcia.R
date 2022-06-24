
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

parser <- add_option(
    parser, c("-m", "--image_file"), default=NULL,
    help="Matching image name.")

parser <- add_option(
    parser, c("-t", "--theta"), default=0.9,
    type="double",
    help="Matching image name.")

parser <- add_option(
    parser, c("-r", "--sankey_rp_topN"), default=5,
    type="double",
    help="Number of top receiver pathways to be shown in the Sankey plot.")

parser <- add_option(
    parser, c("-s", "--sankey_sp_topN"), default=10,
    type="double",
    help="Number of top sender pathways to be shown in the Sankey plot.")

parser <- add_option(
    parser, c("--sankey_short"), default=F, action='store_true',
    type="logical",
    help="Number of top sender pathways to be shown in the Sankey plot.")


########  define input  ########
# args = parse_args(parser)
args = parse_args(
    parser, 
    c(
        '-i', '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spaciapy_test',
        '-j', 'FGFR1_correlated_genes',
        '-m', 'LN_visium_rescaled.jpg',
        '-r', '2',
        '-s', '10',
        '--sankey_short')
    )
if (is.null(args$spacia_results)) {
    stop("Spacia results must be provided. See script usage (--help)")
}

results_path = args$spacia_results
spacia_path = args$spacia_path
receiving_pathway = args$job_id
img_fn = args$image_file
f = args$theta
sankey_rn = args$sankey_rp_topN
sankey_sn = args$sankey_sp_topN
sankey_shortname = args$sankey_short
job_id = args$job_id

model_input_path = file.path(results_path, 'model_input')
exp_sender_fn = file.path(model_input_path, 'exp_sender.json')
metadata_fn = file.path(model_input_path, 'metadata.txt')
interactions_fn = file.path(results_path, 'Interactions.csv')
b_and_fdr_fn = file.path(results_path,'B_and_FDR.csv')
pathway_betas_fn = file.path(results_path,'Pathway_betas.csv')
img_fn = file.path(results_path,img_fn)

exp_receiver_fn = file.path(
    model_input_path, paste(receiving_pathway,'exp_receiver.csv', sep='_')
)
beta_fn = file.path(
    results_path, receiving_pathway, paste(receiving_pathway,'beta.txt', sep='_')
)

########  Load all required plotting functions  ########
source(file.path(spacia_path, "plot_interaction.R"))
source(file.path(spacia_path, "plot_density_ridgeline.R"))
source(file.path(spacia_path, "plot_interaction_circos.R"))
source(file.path(spacia_path, "plot_interaction_direction.R"))
source(file.path(spacia_path, "plot_interaction_proportion.R"))
source(file.path(spacia_path, "plot_interaction_river.R"))
source(file.path(spacia_path, "plot_b.R"))
source(file.path(spacia_path, "plot_sankey.R"))

# Plot posterior b
pdf(
    file.path(
        results_path,
        paste('B_', 'Theta_', as.character(theta), '.pdf', sep='')),
    width = 8, height = 8
    )

plot_b(b_and_fdr_fn, theta)
dev.off()

# Plot pathway interactions sankey
pdf(file.path(results_path,'Pathway Sankey.pdf'),width = 10, height = 6)
plot_sankey(pathway_betas_fn, sankey_rn, sankey_sn, sankey_shortname)
dev.off()

# plot_interaction(metadata_fn,exp_receiver_fn, exp_sender_fn, beta_fn, interactions_fn, 'all', 0.2)

interactions = read.table(interactions_fn, sep=',', header=T)
job_interactions = interactions[(interactions[,1] == job_id) & (interactions$Primary_instance_score>=theta),]
job_interactions = job_interactions[,2:4]
colnames(job_interactions) = c('X','Y','weight')
coordinates = read.table(metadata_fn, header=T, sep='\t')
coordinates$X = coordinates$X * 0.3
coordinates$Y = coordinates$Y * 0.3

plot_interaction_direction(img_fn, job_interactions, coordinates,position_scale = 1, point_scale = 2,
                           sender_col='darkred', reciever_col='steelblue', interaction_col='darkgreen',
                           image_alpha=0.5, cell_alpha=1, interaction_alpha=5, interaction_size=0.1)


