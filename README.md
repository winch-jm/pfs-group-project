# **LeiDAN: A Go Implementation of Leiden Clustering and Weighted KNN for Drug Association Networks**

This repository provides a **pure Go implementation of the Leiden community detection algorithm** together with a **weighted K-nearest neighbors (KNN)** graph builder designed for large-scale transcriptomic data. The full pipeline has been applied to **LINCS L1000 Level 5 signatures** to reproduce the **Drug Association Networks (DANs)** described in:  
> *Musa et al., “A Drug-Centric View of Drug Development: Drug Association Networks,” Scientific Reports (2019).*  
> https://www.nature.com/articles/s41598-019-44291-3

The project includes both the high-performance backend for constructing and clustering similarity graphs **and** an interactive **RShiny dashboard** that visualizes community evolution, enrichment, and metadata structure across Leiden levels.

Group members:
1. [Jeff Winchell](https://github.com/winch-jm/)
2. [Ajay Prabhakar](https://github.com/ajayprab20/)
3. [Nethan Ramachandran](https://github.com/nethanr/)

[Link to Data](https://drive.google.com/drive/u/0/folders/1Ge4OSOuxCWXZXDOTZMt1taMk5nvEdpj9)
---

## 🚀 **Features**

### **Algorithmic + Graph Features**
- **Pure Go Leiden clustering**
  - Fast local moving  
  - Refinement step  
  - Multi-level aggregation  
  - Modularity-based optimization  

- **Weighted KNN graph construction**
  - Cosine-similarity–based  
  - Efficient for large, sparse data  
  - Integrates naturally with CSR structure  

- **Custom CSR Graph**
  - Memory-efficient adjacency matrix  
  - Fast neighbor iteration  
  - Suitable for large LINCS-scale graphs  

### **Data-Driven Features**
- **Reconstruction of Drug Association Networks (DANs)**  
- **Integration of ATC, MoA, and TAS metadata**  
- **Discovery of functional drug similarity beyond predefined labels**

---

## 📊 **Interactive RShiny Dashboard**

This project includes a fully implemented **RShiny dashboard** for visual exploration of the clustering results.  
Key capabilities:

- **Community label evolution viewer**  
  - Track how individual nodes move across Leiden levels  
  - Visualize stability and transitions between communities  

- **TAS, ATC, MoA, and metadata enrichment panels**  
  - Automatically computes enrichment for final communities  
  - Provides interpretable summaries and visualizations  

- **Multiple quality functions**  
  - Compare Leiden vs Louvain  
  - Explore different resolution (γ) settings  

- **Optimization utilities**
  - Visualize how clustering quality varies across γ  
  - Aid in selecting resolution parameters
  - 
---

## 📁 **Repository Structure**
```
├── code       <- go code
├── dashboard  <- rshiny dashboard
├── data       <- lincs l1000 subset
└── notebooks  <- python notebooks used to preprocess data
```
## **How To Run**
1. Download mcf7_subset_pca.csv file and place in /data sub-directory in folder
2. Run command listed below:
```
./code dataset prefix outputDir gridSearch k(optional) gamma(optional)

```
- dataset:  input dataset in .csv format
- prefix: tag you want for graphs and images
- outputDir: filepath that the algorithm results go to
- gridSearch: true or false, if user wants to find optimal hyperparameters
- k: neighbors in initial graph
- gamma: float hyperparameter
Dow
