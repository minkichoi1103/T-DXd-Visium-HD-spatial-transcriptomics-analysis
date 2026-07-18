import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from anndata import AnnData as an
import warnings
from pandas.errors import PerformanceWarning

warnings.simplefilter("ignore", PerformanceWarning)
save_path = "../data/h5ad_55um_merged"
if not os.path.exists(save_path):
    os.makedirs(save_path, exist_ok=True)

root_path = "../data/h5ad_segmented_mf/"
dataset_list = sorted([i for i in os.listdir(root_path) if i.endswith(".h5ad") and not ('B05' in i)])

celltype_colnames = ['B_cells_filtered', 'CAFs_filtered',
                     'Cancer_Epithelial_filtered', 'Endothelial_filtered',
                     'Myeloid_filtered', 'Normal_Epithelial_filtered', 'PVL_filtered',
                     'Plasmablasts_filtered', 'T_cells_filtered']

for idx, dataset in enumerate(dataset_list):
    print(f"Dataset Name: {dataset}")
    adata = sc.read_h5ad(os.path.join(root_path, dataset))
    
    print(f"Running Gridding Process")
    pixel_size = adata.uns['spatial']['Visium_HD']['scalefactors']['microns_per_pixel']
    ## Grid-based aggregation of image-based ST: divide coordinates by x_bins and y_bins and aggregate
    spatial_coord = adata.obsm['spatial']
    # Find the x and y coordinate arrays
    x_coord = spatial_coord[:,0] * pixel_size
    y_coord = spatial_coord[:,1] * pixel_size
    # Find the coordinates that equally divides the x and y axis into x_bins and y_bins
    x_div_arr = np.linspace(np.min(x_coord), np.max(x_coord), num=int((np.max(x_coord)-np.min(x_coord))//55), endpoint=False)[1:]
    y_div_arr = np.linspace(np.min(y_coord), np.max(y_coord), num=int((np.max(y_coord)-np.min(y_coord))//55), endpoint=False)[1:]
    # Assigning the grid column and row number to each transcript based on the coordinates by x_div_arr and y_div_arr
    adata.obs['grid_array_col'] = np.searchsorted(x_div_arr, x_coord, side='right')
    adata.obs['grid_array_row'] = np.searchsorted(y_div_arr, y_coord, side='right')
    
    # Calculate means between boundaries and the first/last bin lines
    x_boundary_midpoints = [
        (np.min(x_coord) + x_div_arr[0]) / 2,  # Between boundary and first bin
        (np.max(x_coord) + x_div_arr[-1]) / 2  # Between last bin and boundary
    ]
    y_boundary_midpoints = [
        (np.min(y_coord) + y_div_arr[0]) / 2,  # Between boundary and first bin
        (np.max(y_coord) + y_div_arr[-1]) / 2  # Between last bin and boundary
    ]

    # Calculate midpoints for inside the grid
    array_col_center = (x_div_arr[:-1] + x_div_arr[1:]) / 2
    array_row_center = (y_div_arr[:-1] + y_div_arr[1:]) / 2
    # Combine midpoints and boundaries into a single array
    array_col_with_boundaries = np.concatenate(([x_boundary_midpoints[0]], array_col_center, [x_boundary_midpoints[1]]))
    array_row_with_boundaries = np.concatenate(([y_boundary_midpoints[0]], array_row_center, [y_boundary_midpoints[1]]))
    array_col_with_boundaries, array_row_with_boundaries = array_col_with_boundaries / pixel_size, array_row_with_boundaries / pixel_size
    
    array_col_index = np.arange(len(array_col_with_boundaries))
    array_row_index = np.arange(len(array_row_with_boundaries))

    col_idx, row_idx = np.meshgrid(array_col_index, array_row_index, indexing="ij")
    col_coords, row_coords = np.meshgrid(array_col_with_boundaries, array_row_with_boundaries, indexing="ij")
    
    # Flatten the meshgrid arrays for the DataFrame
    flat_col_idx = col_idx.flatten()
    flat_row_idx = row_idx.flatten()
    flat_col_coords = col_coords.flatten()
    flat_row_coords = row_coords.flatten()
    
    # Create the DataFrame
    df_meshgrid = pd.DataFrame({
        "Array_Col_Index": flat_col_idx,
        "Array_Row_Index": flat_row_idx,
        "Array_Col_Coordinate": flat_col_coords,
        "Array_Row_Coordinate": flat_row_coords
    })
    df_meshgrid.index = df_meshgrid['Array_Col_Index'].astype(str) + '_' + df_meshgrid['Array_Row_Index'].astype(str)
    
    ## Normalize the transcript number in each grid by total count in the cell
    tx_by_cell_grid = pd.concat([adata.obs.loc[:,['grid_array_col','grid_array_row']], 
                                 pd.DataFrame(adata.X.toarray(), 
                                             index=adata.obs_names, columns=adata.var_names)], axis=1)
    # Generate normalization count matrix by grid
    grid_tx_count = tx_by_cell_grid.groupby(['grid_array_col','grid_array_row']).sum()
    
    # Saving grid barcode and gene symbol names
    var_names = grid_tx_count.columns
    grid_metadata = grid_tx_count.index.to_frame(name=['array_col','array_row'])
    grid_metadata.index = grid_metadata['array_col'].astype(str) + '_' + grid_metadata['array_row'].astype(str)
    # Log transformation of grid based count
    grid_tx_count = (sparse.csr_matrix(grid_tx_count, dtype=np.float32))

    grid_celltype = adata.obs.loc[:,['grid_array_col','grid_array_row']+celltype_colnames]
    grid_celltype = grid_celltype.pivot_table(index=['grid_array_col','grid_array_row'], aggfunc=['sum']).fillna(0)
    grid_celltype.columns = grid_celltype.columns.to_frame()[1].values
    # Assign index names to the dataframe
    grid_index = grid_celltype.index.to_frame()
    grid_celltype.index = grid_index['grid_array_col'].astype(str) + '_' + grid_index['grid_array_row'].astype(str)
    # Modify metadata to contain cell type information in each grid
    grid_metadata = grid_metadata.join(grid_celltype, how='left').fillna(0)

    ## Generating grid-based image-based ST anndata
    sp_adata_grid = an(X = grid_tx_count, obs=grid_metadata)
    sp_adata_grid.var_names = var_names

    sp_adata_grid.obs = sp_adata_grid.obs.join(df_meshgrid)
    sp_adata_grid.uns['spatial'] = adata.uns['spatial'].copy()
    sp_adata_grid.obsm['spatial'] = sp_adata_grid.obs[['Array_Col_Coordinate','Array_Row_Coordinate']].to_numpy()
    
    df_core = adata.obs[['grid_array_col','grid_array_row','core_name','orig_cluster','subcluster']]
    df_core['grid_name'] = df_core['grid_array_col'].astype(str) + '_' + df_core['grid_array_row'].astype(str)
    df_core.drop_duplicates(subset = ['grid_array_col','grid_array_row'], inplace=True)
    df_core.set_index('grid_name', inplace=True, drop=True)
    sp_adata_grid.obs = sp_adata_grid.obs.join(df_core)
    sp_adata_grid.write_h5ad(os.path.join(save_path, f'sp_grid_{dataset}.h5ad'))
    print("-"*30)