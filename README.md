# AMP-PD_SciData

This relates to the publication: **A multi-region single nucleus transcriptomic atlas of Parkinson’s disease**

doi: TBA

## Creating virtual environment

Using [anaconda](https://conda.io/projects/conda/en/latest/index.html), and the *env.yaml* file provided here the following scripts can be run.

```bash
conda env create -f env.yml
```

### pd-metadata-overview.R

This script when run as:

```bash
Rscript pd-metadata-overview.R
```

produces the following output files:

- output_dir/plot_age_by_diagnosis.pdf
- output_dir/plot_pd_phenotypes_correlation_spearman.pdf
- output_dir/plot_age_by_sex.pdf
- output_dir/plot_race.pdf
- output_dir/plot_avgCellCountPerDiss_byBrainRegion.pdf
- output_dir/plot_sampleCount_byBrainRegion.pdf
- output_dir/plot_braak_by_brainbank.pdf
- output_dir/plot_sampleCount_byBrainRegion_Upset.pdf
- output_dir/plot_brainbank.pdf
- output_dir/plot_sex_by_brainbank.pdf
- output_dir/plot_cellCount2.pdf
- output_dir/plot_sex.pdf
- output_dir/plot_diagnosis_by_brainbank.pdf
- output_dir/plot_totalCellCount_byBrainRegion.pdf
- output_dir/plot_diagnosis.pdf

> This script is designed to run as is when the whole git repo is downloaded.
> Otherwise please fix **line 12** to reflect your working directory.

### qtltools_plot_dash_SciData.py

This script when run as:

```bash
python qtltools_plot_dash_SciData.py
```

produces the following outputs:

- output_dir/SciData_fig4h.pdf

### SciDatafigure_figs3abc.py

This script when run as:

```bash
python SciDatafigure_figs3abc.py
```

produces the following outputs:

- output_dir/SciData_fig3a.pdf  
- output_dir/SciData_fig3b_1.pdf  
- output_dir/SciData_fig3b_2.pdf  
- output_dir/SciData_fig3c.pdf
