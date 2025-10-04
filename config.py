import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
from skimage import io
from scipy.io import mmread

# Genes.csv - List of genes in the feature matrix
# FeatureMatrix.mtx - Feature matrix, stored as a sparse matrix
# ClusterGeneExpression.csv - Gene expression per cell-type cluster
# SpotClusterMembership.csv - Cell type proportions per spot
# SpotPositions.csv - Spot positions and radiuses
# Images/scalefactors_json.json - Scale factors between spot positions
# Images/tissue_hires_image.png - H&E image

data_dir = "data/"
genes_file = data_dir + "Genes.csv"
matrix_file = data_dir + "FeatureMatrix.mtx"
clusterMembership_file = data_dir + "SpotClusterMembership.csv"
positions_file = data_dir + "SpotPositions.csv"

imgs_path = data_dir + "images/"
tissue_img = imgs_path + "tissue_hires_images.png"

def load_data():
    genes = pd.read_csv(genes_file, header=0, names=['gene', 'feature_type'])
    spotPosition = pd.read_csv(positions_file)
    clusterMembership = pd.read_csv(clusterMembership_file)
    spotPosition['y'] = 17244 - spotPosition['y']
    data = pd.merge(spotPosition, clusterMembership, on='barcode')
    
    expression_matrix = mmread(matrix_file).T.tocsr()
    print(f"Number of genes: {len(genes)}")
    print(f"Number of spots: {len(spotPosition)}")

    adata = ad.AnnData(
        X = expression_matrix
    )
    adata.var_names = genes['gene'].values
    adata.obsm['spatial'] = spotPosition[['x', 'y']].values
    adata.obs = clusterMembership.set_index(spotPosition.index)
    
    print(f"data loaded: {adata}")

    return adata

def process_data(adata):
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    return adata

def draw_visual(adata):
    adata.obs['dominant_cluster'] = adata.obs[['X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9']].idxmax(axis=1)
    sc.pl.spatial(adata, color='dominant_cluster', spot_size=50)

    top_genes = adata.var_names[adata.var['highly_variable']][:10].tolist()
    sc.pl.spatial(adata, color=top_genes, ncols=3)