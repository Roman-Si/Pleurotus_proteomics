# Pleurotus_proteomics
This repository contains code and output files for the proteomics analysis of *Pleurotus pulmonarius* LGAM 28684 grown on three different substrates: xylose, corn stover, and beechwood.


The analysis is in [Quarto](https://quarto.org/) documents (`.qmd`), with rendered HTML outputs available in the [`docs/`](docs/) directory. These HTML files can be previewed in-browser using [htmlpreview.github.io](https://htmlpreview.github.io/), avoiding the need to download them.

The analysis is in Quarto documents (.qmd), with rendered HTML outputs available in the docs/ directory. These HTML files can be previewed in-browser using htmlpreview.github.io, avoiding the need to download them.

### Rendered reports

* Quantitative proteomics analysis:
  [prolfqua.html](https://htmlpreview.github.io/?https://github.com/Roman-Si/Pleurotus_proteomics/blob/main/docs/prolfqua.html)

* Exploratory analysis of batch effects:
  [batches\_qc.html](https://htmlpreview.github.io/?https://github.com/Roman-Si/Pleurotus_proteomics/blob/main/docs/batches_qc.html)

---
  
## Packages used

* For protein functional annotation: [ProtAnnoScripts](https://github.com/Roman-Si/ProtAnnoScripts)
* For quantitative proteomic analysis: [Rscripts\_for\_proteomics](https://github.com/Roman-Si/Rscripts_for_proteomics)
* All plots were created in R using `ggplot2` (scripts in [`code/`](code/)) and further edited in [Inkscape](https://inkscape.org/).

---

## To-do

* Add `download_data.sh` to retrieve raw and processed data from PRIDE once published

---
