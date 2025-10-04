import config

# Quick spatial plot
import matplotlib.pyplot as plt

def main():
    adata = config.load_data()# Check what's in the cluster columns
    print(adata.obs.head())
    print(adata.obs.describe())

    # Check expression matrix statistics
    print(f"Expression range: {adata.X.min()} to {adata.X.max()}")
    print(f"Sparsity: {(adata.X == 0).sum() / adata.X.size * 100:.2f}%")

    adata = config.process_data(adata)
    config.draw_visual(adata)

    plt.scatter(adata.obsm['spatial'][:, 0], adata.obsm['spatial'][:, 1], s=50)
    plt.title("Tissue spatial layout")
    plt.show()

if __name__ == "__main__":
    main()