#!/usr/bin/env python
# coding: utf-8

import plotly.express as px
import pandas as pd, os


output_dir = "output_dir/"
dry_lab_dir = "dry_lab/"

to_plot = pd.read_csv(os.path.join(dry_lab_dir, "SciData_gt_demux.csv"))

x_min = to_plot.loc[to_plot['cell_counts'] > 800, 'perc_het_consistent'].min()
y_min = to_plot.loc[to_plot['cell_counts'] > 800, 'perc_hom_consistent'].min()

to_plot['color'] = to_plot['color'].astype('str').astype('category')
to_plot['color'] = to_plot['color'].map({'Expected Mismatch': 'Mismatch', 'Expected Match': 'Match', 
                                         'Removed Donor (gt issue)': 'Mismatch', 'Expected Match (QC fix)': 'Match',})
to_plot = to_plot.loc[to_plot['cell_counts'] > 800, :]
fig = px.scatter(
    to_plot, x="perc_het_consistent", y="perc_hom_consistent", 
    color="color",
    color_discrete_map = {'Mismatch': 'blue', 
                          'Match': 'green',
                         },
    hover_data=['pool_name', 'donor_name', 'SampleID']
)

fig.update(layout_yaxis_range = [y_min,1])
fig.update(layout_xaxis_range = [x_min,1])
fig.update_layout(
    margin=dict(l=20, r=20, t=20, b=20),
    paper_bgcolor="LightSteelBlue",
    minreducedwidth=1000,
    minreducedheight=1300,
)
fig.update_layout(width=int(1300))
fig.write_image(os.path.join(output_dir, "SciData_fig4h.pdf"))
