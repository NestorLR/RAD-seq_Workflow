################################################################################
# Population genomics of a cosmopolitan weed provides insights into its 
# local adaptation and recent demographic history
# Year: 2026
# FEEMs analysis
################################################################################

#####----------------------------------------------------------------------#####

library(dggridR)
library(sf)

#####----------------------------------------------------------------------#####

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from pandas_plink import read_plink
from sklearn.impute import SimpleImputer

from feems import SpatialGraph, Viz
from feems.utils import prepare_graph_inputs

bim, fam, G = read_plink("Conyza")
genotypes = SimpleImputer(strategy="mean").fit_transform(np.array(G).T)

coord = np.loadtxt("Conyza.coord")
outer = np.loadtxt("Conyza.outer")
grid_path = "world.shp"

outer, edges, grid, _ = prepare_graph_inputs(
  coord=coord,
  ggrid=grid_path,
  translated=False,
  buffer=0,
  outer=outer
)

sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)
sp_graph.fit(lamb=2.0, lamb_q=10.0, optimize_q="n-dim")

projection = ccrs.AzimuthalEquidistant(central_longitude=-100)
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)

v = Viz(
  ax,
  sp_graph,
  projection=projection,
  edge_width=0.5,
  edge_alpha=1,
  edge_zorder=100,
  sample_pt_size=20,
  obs_node_size=7.5,
  sample_pt_color="black",
  cbar_font_size=10
)

v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()

plt.show()



#####---------------------------------------------------------------------------#####
