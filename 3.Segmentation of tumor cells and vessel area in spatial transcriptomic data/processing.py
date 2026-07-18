from skimage.filters import threshold_otsu
from scipy.ndimage import binary_erosion
from scipy.sparse import csr_matrix

import numpy as np
import pandas as pd
import scanpy as sc


def find_cell_areas(adata, cell_type='Cancer_Epithelial_filtered'):
    # Extract coordinates and cell-type specific values
    cols = adata.obs['array_col']
    rows = adata.obs['array_row']
    values = adata.obs[cell_type]
    
    # Determine matrix bounds
    min_row, max_row = rows.min(), rows.max()
    min_col, max_col = cols.min(), cols.max()
    
    n_rows = int(max_row - min_row + 1)
    n_cols = int(max_col - min_col + 1)
    
    # Initialize a matrix with NaN values (positions with no data remain NaN)
    matrix = np.full((n_rows, n_cols), np.nan)
    
    row_index = range(int(min_row), int(max_row) + 1)
    col_index = range(int(min_col), int(max_col) + 1)
    
    # Create a DataFrame using the initialized matrix
    df = pd.DataFrame(matrix, index=row_index, columns=col_index)
    
    # Fill the DataFrame with observed values
    for r, c, v in zip(rows, cols, values):
        if np.isnan(df.loc[r, c]):
            df.loc[r, c] = v
    df = df.fillna(0)
    
    # Calculate the Otsu threshold using the cell-type values (drop NaN first)
    cell_values = adata.obs[cell_type].dropna().values
    
    if len(cell_values) == 0 or np.all(np.isnan(cell_values)):
        # No valid values: skip thresholding, set all cells as 'others'
        import warnings
        warnings.warn(f"No valid (non-NaN) values found for '{cell_type}'. Skipping threshold computation.")
        
        adata.obs[f"edge_{cell_type}"] = False
        adata.obs[f"is_{cell_type}"] = False
        adata.obs[f"area_{cell_type}"] = 'others'
        return adata
    
    thresh = threshold_otsu(np.array(cell_values, dtype=float))
    
    # Create a binary mask and determine the edge using erosion
    binary = df > thresh
    eroded = binary_erosion(binary.values)
    edge_mask = binary & ~eroded
    
    # Convert the edge mask to a DataFrame and rename the edge column
    edge_series = edge_mask.stack().rename(f"edge_{cell_type}").reset_index()
    edge_series.columns = ["array_row", "array_col", f"edge_{cell_type}"]
    
    # Merge the edge information back to adata.obs
    adata_obs_updated = adata.obs.merge(edge_series, on=["array_row", "array_col"], how="left")
    # Fill NaN values in the edge column with False
    adata_obs_updated[f"edge_{cell_type}"] = adata_obs_updated[f"edge_{cell_type}"].fillna(False)

    adata_obs_updated.index = adata.obs_names

    # Update adata.obs with the merged DataFrame
    adata.obs = adata_obs_updated
    
    # Create a binary column indicating if a cell is above the threshold
    # NaN values in cell_type column are treated as below threshold (False)
    adata.obs[f"is_{cell_type}"] = np.where(
        adata.obs[cell_type].notna() & (adata.obs[cell_type] > thresh), True, False
    )
    
    # Initialize an area column with a default value 'others'
    adata.obs[f"area_{cell_type}"] = 'others'
    # Label cells that are cancer (not on the edge)
    adata.obs.loc[(~adata.obs[f"edge_{cell_type}"]) & (adata.obs[f"is_{cell_type}"]), f"area_{cell_type}"] = 'cell_area'
    # Label cells that are on the edge
    adata.obs.loc[(adata.obs[f"edge_{cell_type}"]) & (adata.obs[f"is_{cell_type}"]), f"area_{cell_type}"] = 'edge'
    
    return adata



def select_specific_cell(adata, cell_type='Cancer_Epithelial_filtered'):
    """
    Filters the AnnData object to include only cells of a specific type and applies
    cell weights to each row of the expression matrix.
    
    Parameters:
    -----------
    adata : AnnData
        The input AnnData object.
    cell_type : str
        The cell type name. The following columns should exist in adata.obs:
            - 'is_{cell_type}': Boolean column for filtering cells.
            - '{cell_type}': Column containing cell weight values.
            - 'edge_{cell_type}': Boolean column indicating whether to apply the weight (True uses the weight, False uses 1).
    
    Returns:
    --------
    AnnData
        A filtered AnnData object with the cell weights applied to the expression matrix.
    """
    # Create a copy of the AnnData object to avoid modifying the original data
    adata = adata.copy()
    
    # Filter cells based on the 'is_{cell_type}' column in adata.obs
    adata_sub = adata[adata.obs[f'is_{cell_type}']].copy()
    
    # Retrieve cell weight from the '{cell_type}' column
    cell_weight = adata_sub.obs[f'{cell_type}']
    # Use cell_weight if 'edge_{cell_type}' is True, otherwise use 1
    cell_weight_ = np.where(adata_sub.obs[f'edge_{cell_type}'], cell_weight, 1)
    
    # Multiply each row of the expression matrix by the corresponding cell weight
    if hasattr(adata_sub.X, "multiply"):
        weighted_X = adata_sub.X.multiply(cell_weight_.reshape(-1, 1))
    else:
        weighted_X = csr_matrix(adata_sub.X * cell_weight_.reshape(-1, 1))
    
    # Convert COO matrices to CSR to avoid future errors
    if hasattr(weighted_X, 'format') and weighted_X.format == 'coo':
        weighted_X = weighted_X.tocsr()
    
    # Update the AnnData object with the weighted expression matrix
    adata_sub.X = weighted_X
    return adata_sub
