# Standard library imports
import os
from dataclasses import dataclass
from typing import List, Optional, Tuple

# Third-party library imports
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.path import Path
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.spatial import ConvexHull

def process_spacia_data(receiving_cell: str, sending_cell: str, directory: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process Spacia data files from a given directory, combining beta and pip files.

    This function walks through the specified directory, processes 'betas.csv' and 'pip.csv' files,
    and combines them into two separate DataFrames.

    Parameters:
    receiving_cell (str): Identifier for the receiving cell type
    sending_cell (str): Identifier for the sending cell type
    directory (str): Path to the directory containing the data files

    Returns:
    Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
        - combined_betas: DataFrame with combined data from all 'betas.csv' files
        - combined_pip: DataFrame with combined data from all 'pip.csv' files
    """
    def read_csv_files(file_paths: List[str]) -> pd.DataFrame:
        """Read and combine multiple CSV files into a single DataFrame."""
        return pd.concat([pd.read_csv(file) for file in file_paths], ignore_index=True)

    beta_files = []
    pip_files = []

    # Walk through the directory and collect file paths
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if file.endswith('betas.csv'):
                beta_files.append(file_path)
            elif file.endswith('pip.csv'):
                pip_files.append(file_path)

    # Process beta files
    if beta_files:
        combined_betas = read_csv_files(beta_files)
        combined_betas['sending_cell'] = sending_cell
        combined_betas['receiving_cell'] = receiving_cell
    else:
        combined_betas = pd.DataFrame()

    # Process pip files
    combined_pip = read_csv_files(pip_files) if pip_files else pd.DataFrame()

    return combined_betas, combined_pip

# Example usage:
# betas_df, pip_df = process_spacia_data("TumorCells", "EndothelialCells", "/path/to/data/directory")



def calculate_weighted_activation_scores(expression_data: pd.DataFrame, 
                                         cell_metadata: pd.DataFrame, 
                                         beta_values: pd.DataFrame, 
                                         pi_scores: pd.DataFrame, 
                                         proximity_dict: dict, 
                                         receiver_type: str, 
                                         sender_type: str, 
                                         sending_genes: list, 
                                         receiving_genes: list,
                                         score_name: str) -> pd.DataFrame:
    """
    Calculate weighted activation scores for cell-cell interactions based on gene expression,
    interaction strengths, and spatial proximity.

    Parameters:
    expression_data (pd.DataFrame): Gene expression data (cells as rows, genes as columns)
    cell_metadata (pd.DataFrame): Cell metadata including cell type
    beta_values (pd.DataFrame): Beta values for gene interactions
    pi_scores (pd.DataFrame): Primary instance scores
    proximity_dict (dict): Dictionary of receiver cells and their nearby sender cells
    receiver_type (str): Cell type of the receiver cells
    sender_type (str): Cell type of the sender cells
    sending_genes (list): List of sending genes to consider
    receiving_genes (list): List of receiving genes to consider
    score_name (str): Name of the activation score

    Returns:
    pd.DataFrame: Final data with weighted activation scores for each receiver cell
    """

    def calculate_activation_scores(proximity_dict: dict, 
                                    beta_filtered: pd.DataFrame, 
                                    exp_data: pd.DataFrame) -> pd.DataFrame:
        """Calculate activation scores for cells using the pre-computed proximity dictionary."""
        activation_scores = []
        unique_sender_cells = set(cell for cells in proximity_dict.values() for cell in cells)
        
        for sender_cell in unique_sender_cells:
            cell_expression = exp_data.loc[sender_cell]
            merged_data = beta_filtered.merge(cell_expression, left_on='sending_gene', right_index=True, how='left')
            merged_data['expression'] = merged_data.iloc[:, -1]
            merged_data['score'] = np.maximum(0, merged_data['avg_beta_sampled']) * merged_data['expression']
            
            activation_scores.append({
                'sending_cell': sender_cell,
                'activation_score': merged_data['score'].sum()
            })
        
        return pd.DataFrame(activation_scores)

    # Filter beta values
    beta_filtered = beta_values[
        (beta_values['b'] < 0) & 
        (beta_values['sending_gene'].isin(sending_genes)) & 
        (beta_values['receiving_gene'].isin(receiving_genes)) &
        (beta_values['sending_cell'] == sender_type) &
        (beta_values['receiving_cell'] == receiver_type)
    ]

    # Calculate activation scores
    activation_scores = calculate_activation_scores(proximity_dict, beta_filtered, expression_data)

    # Summarize PI scores
    pi_summary = pi_scores.groupby(['receiving_cell', 'sending_cell'])['avg_primary_instance_score'].mean().reset_index()

    # Merge activation scores with PI scores
    merged_scores = pi_summary.merge(activation_scores, on='sending_cell', how='left')
    merged_scores['weighted_score'] = merged_scores['avg_primary_instance_score'] * merged_scores['activation_score']

    # Calculate final scores for each receiving cell
    final_scores = []
    for receiver, senders in proximity_dict.items():
        cell_scores = merged_scores[
            (merged_scores['receiving_cell'].astype(str) == str(receiver)) & 
            (merged_scores['sending_cell'].isin([str(cell) for cell in senders]))
        ]
        if not cell_scores.empty:
            final_scores.append({
                'receiving_cell': receiver,
                score_name: cell_scores['weighted_score'].sum()
            })

    final_scores_df = pd.DataFrame(final_scores)

    # Merge with receiver cell metadata
    receiver_cells = cell_metadata[cell_metadata['cell_type'] == receiver_type].reset_index()
    final_data = receiver_cells.merge(final_scores_df, 
                                      left_on='index', 
                                      right_on='receiving_cell',
                                      how='inner')

    return final_data

# Example usage:
# result = calculate_weighted_activation_scores(
#     expression_data=expression_data,
#     cell_metadata=cell_metadata,
#     beta_values=beta_values,
#     pi_scores=pi_scores,
#     proximity_dict=proximity_dict,
#     receiver_type="Tumor Cells",
#     sender_type="Endothelial Cells",
#     sending_genes=["HGF", "WNT5A", "FGF2", "IL6", "CXCL8", "FGF1", "TGFB1", "TGFB2"],
#     receiving_genes=["JAK1", "AKT2", "SMO", "CTNNB1", "SMAD2", "NFKB2"],
#     score_name="Endothelial_Activation_Score"    
# )




@dataclass 
class SpatialScore:
    """
    Class to hold spatial scoring results and perform density estimation
    
    Attributes:
        score_data (pd.DataFrame): DataFrame containing cell coordinates and scores
        score_name (str): Name of the score (e.g. "EMT", "Activation")
        mask (Optional[np.ndarray]): Custom binary mask for the tissue area
        xi (Optional[np.ndarray]): Optional x coordinates of the grid
        yi (Optional[np.ndarray]): Optional y coordinates of the grid
    """
    score_data: pd.DataFrame
    score_name: str
    mask: Optional[np.ndarray] = None
    xi: Optional[np.ndarray] = None
    yi: Optional[np.ndarray] = None
    
    def __post_init__(self):
        self.x_coords = self.score_data['X'].values
        self.y_coords = self.score_data['Y'].values
        self.scores = self.score_data[self.score_name].values
        
        # Initialize as None, compute only when needed
        self._kde = None
        self._density = None
        self._hull_mask = None
        
        # Validate xi and yi if provided
        if (self.xi is not None) != (self.yi is not None):
            raise ValueError("Both xi and yi must be provided together")
        
        if self.xi is not None:
            if self.xi.shape != self.yi.shape:
                raise ValueError("xi and yi must have the same shape")
            self._xi = self.xi
            self._yi = self.yi
        else:
            self._xi = None
            self._yi = None

    def compute_kde(self, bandwidth: float = 0.1, grid_size: int = 300):
        """
        Compute kernel density estimation
        
        Parameters:
            bandwidth: KDE bandwidth parameter
            grid_size: Size of the grid (used only if xi/yi not provided at init)
        """
        # Compute KDE
        self._kde = gaussian_kde(
            np.vstack([self.x_coords, self.y_coords]),
            weights=self.scores,
            bw_method=bandwidth
        )
        
        # Create grid if not provided
        if self._xi is None:
            x_min, x_max = self.x_coords.min(), self.x_coords.max()
            y_min, y_max = self.y_coords.min(), self.y_coords.max()
            self._xi, self._yi = np.mgrid[x_min:x_max:grid_size*1j, 
                                         y_min:y_max:grid_size*1j]
        
        # Compute density
        positions = np.vstack([self._xi.ravel(), self._yi.ravel()])
        self._density = self._kde(positions).reshape(self._xi.shape)
        
        # Generate automatic hull mask if no custom mask provided
        if self.mask is None:
            self._generate_hull_mask()
        else:
            if self.mask.shape != self._xi.shape:
                raise ValueError("Provided mask shape must match grid shape")
            self._hull_mask = self.mask
            
        return self

    def _generate_hull_mask(self):
        """Generate convex hull mask for the tissue area"""
        coords = np.column_stack([self.x_coords, self.y_coords])
        hull = ConvexHull(coords)
        hull_points = coords[hull.vertices]
        
        # Create mask using hull
        points = np.vstack([self._xi.ravel(), self._yi.ravel()]).T
        hull_path = Path(hull_points)
        self._hull_mask = hull_path.contains_points(points).reshape(self._xi.shape)
    
    
    def get_masked_density(self) -> np.ndarray:
        """Return density values masked by tissue area"""
        if self._density is None:
            self.compute_kde()
        return np.where(self._hull_mask, self._density, np.nan)
    
    def plot(self, cmap: str = 'viridis', title: Optional[str] = None):
        """Plot the spatial distribution of scores"""
        if self._density is None:
            self.compute_kde()
            
        plt.figure(figsize=(10, 8))
        plt.contourf(self._xi, self._yi, self.get_masked_density(), 
                    cmap=cmap, alpha=0.6)
        
        if title:
            plt.title(title)
        plt.colorbar(label=f'KDE of {self.score_name}')
        plt.axis('off')
        return plt.gca()


# Helper function to create common grid for multiple scores
def create_common_grid(score_objects: List[SpatialScore], 
                      grid_size: int = 300) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a common grid for multiple spatial scores
    
    Parameters:
        score_objects: List of SpatialScore objects
        grid_size: Size of the grid
        
    Returns:
        xi, yi: Grid coordinates
    """
    # Find common bounds
    x_min = min(obj.x_coords.min() for obj in score_objects)
    x_max = max(obj.x_coords.max() for obj in score_objects)
    y_min = min(obj.y_coords.min() for obj in score_objects)
    y_max = max(obj.y_coords.max() for obj in score_objects)
    
    # Create common grid
    xi, yi = np.mgrid[x_min:x_max:grid_size*1j, 
                      y_min:y_max:grid_size*1j]
    
    return xi, yi

def compute_correlation(score1: SpatialScore, 
                       score2: SpatialScore,
                       method: str = 'spearman') -> Tuple[float, float]:
    """
    Compute correlation between two spatial scores
    
    Parameters:
        score1, score2: SpatialScore objects
        method: 'spearman' or 'pearson'
        
    Returns:
        correlation coefficient and p-value
    """
    if score1._density is None:
        score1.compute_kde()
    if score2._density is None:
        score2.compute_kde()
        
    # Get masked values
    masked1 = score1.get_masked_density().ravel()
    masked2 = score2.get_masked_density().ravel()
    
    # Remove NaN values
    valid_idx = ~(np.isnan(masked1) | np.isnan(masked2))
    masked1 = masked1[valid_idx]
    masked2 = masked2[valid_idx]
    
    if method == 'spearman':
        return stats.spearmanr(masked1, masked2)
    elif method == 'pearson':
        return stats.pearsonr(masked1, masked2)
    else:
        raise ValueError("Method must be 'spearman' or 'pearson'")

def create_common_mask(score_objects: List[SpatialScore], xi: np.ndarray, yi: np.ndarray) -> np.ndarray:
    """
    Create a common mask using all cell coordinates from all score objects
    
    Parameters:
        score_objects: List of SpatialScore objects
        xi, yi: Grid coordinates
        
    Returns:
        mask: Common binary mask
    """
    # Combine all cell coordinates
    all_x_coords = np.concatenate([obj.x_coords for obj in score_objects])
    all_y_coords = np.concatenate([obj.y_coords for obj in score_objects])
    
    # Create convex hull from all coordinates
    coords = np.column_stack([all_x_coords, all_y_coords])
    hull = ConvexHull(coords)
    hull_points = coords[hull.vertices]
    
    # Create mask using hull
    points = np.vstack([xi.ravel(), yi.ravel()]).T
    hull_path = Path(hull_points)
    return hull_path.contains_points(points).reshape(xi.shape)

# Example usage:
# """
# 
# endo_score = pd.DataFrame({
#     'index': ['cell_2', 'cell_8', ...],  # object (str)
#     'X': [23.5, 45.2, ...],              # float64
#     'Y': [-12.3, 8.7, ...],              # float64
#     'cell_type': ['Tumor_Cells', 'Tumor_Cells', ...],  # object (str)
#     'receiving_cell': ['cell_5', 'cell_8', ...],  # object (str)
#     'Endothelial_Activation_Score': [0.85, 0.92, ...]  # float64
# })
# 
# fibro_score = pd.DataFrame({
#     'index': ['cell_4', 'cell_10', ...],           # object (str)
#     'X': [23.5, 45.2, ...],                       # float64
#     'Y': [-12.3, 8.7, ...],                       # float64
#     'cell_type': ['Tumor_Cells', 'Tumor_Cells',...],   # object (str)
#     'receiving_cell': ['cell_45', 'cell_23', ...],# object (str)
#     'Fibroblast_Activation_Score': [0.80, 0.32, ...]  # float64
# }
#
# # Create objects without computing KDE
# endo_spatialscore = SpatialScore(endo_score, 'Endothelial_Activation_Score')
# fibro_spatialscore = SpatialScore(fibro_score, 'Fibroblast_Activation_Score')
#
# # Create common grid
# xi, yi = create_common_grid([endo_spatialscore, fibro_spatialscore])
# common_mask = create_common_mask([endo_spatialscore, fibro_spatialscore, tumor_spatialscore], xi, yi)
# 
# # Create new objects with common grid
# endo_spatialscore = SpatialScore(endo_score, 'Endothelial_Activation_Score' xi=xi, yi=yi)
# fibro_spatialscore = SpatialScore(fibro_score, 'Fibroblast_Activation_Score' xi=xi, yi=yi)
# 
# # Compute KDE and correlations
# endo_spatialscore.compute_kde()
# fibro_spatialscore.compute_kde()

# corr, pval = compute_correlation(endo_spatialscore, fibro_spatialscore)
# """


