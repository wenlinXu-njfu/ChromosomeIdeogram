[![PyPI - Version](https://img.shields.io/pypi/v/chromosome-ideogram)](https://pypi.org/project/chromosome-ideogram)

# ChromosomeIdeogram
A python package for drawing genome-wide data on ideograms.

## Install
```shell
pip install chromosome-ideogram --upgrade
```

## Usage example
```python
import matplotlib.pyplot as plt
from matplotlib import gridspec
from chromosome_ideogram import ChromosomeIdeogram

plt.rcParams['font.size'] = 12
plt.rcParams['pdf.fonttype'] = 42

fig = plt.figure(figsize=(18, 15))
gs = gridspec.GridSpec(nrows=12, ncols=4)
axes1 = fig.add_subplot(gs[:10, :])
axes2 = fig.add_subplot(gs[10:, :3])  # Used to display the detailed karyotype of a certain chromosome (such as chromosome 14)

ci1 = ChromosomeIdeogram(
    karyotype='example/karyotype.xls',
    radius=0.4,  # Radius of chromosome
    space=4  # Chromosome center coordinate interval (radius x space)
)
ci1.draw_chromh(axes=axes1, centromere_angle=60, facecolor='lightgrey', edgecolor=None, linewidth=0)
ci1.annotations(axes=axes1, annotations='example/satellite_DNA.bed', height=None)
# show rDNA
ci1.mark_loci(axes=axes1, markers='example/hapA_rDNA.bed', pos='bottom', size=None)  # Draw the rDNA beneath haplotype A chromosomes
ci1.mark_loci(axes=axes1, markers='example/hapB_rDNA.bed', pos='top', size=None)  # Draw the rDNA above haplotype B chromosomes

# Create a second karyotype ideogram object (containing only chromosome 14)
ci2 = ChromosomeIdeogram(
    karyotype='example/Chr14_karyotype.xls',
    radius=0.4,
    space=3
)
ci2.draw_chromh(axes=axes2, centromere_angle=60, facecolor='white', edgecolor='#4d4d4d', linewidth=1)
ci2.annotations(axes=axes2, annotations='example/satellite_DNA.bed', height=None)
ci2.mark_loci(axes=axes2, markers='example/hapA_rDNA.bed', pos='bottom', size=None)
ci2.mark_loci(axes=axes2, markers='example/hapB_rDNA.bed', pos='top', size=None)
ci2.draw_synteny(axes=axes2, synteny_file='example/Chr14_SV.bed') # show synteny

axes1.legend(handles=ci1.handles, loc=(0.85, 0.85))  # Set axes1 legend
axes2.legend(handles=ci2.handles, loc=(1.05, 0.1), ncol=2)  # You can also set axes2 legend

plt.subplots_adjust(hspace=9)  # Adjust the vertical distance between the two subgraphs

plt.savefig('example/satellite_DNA.pdf', bbox_inches='tight')
plt.savefig('example/satellite_DNA.png', bbox_inches='tight')
```
![image](https://github.com/wenlinXu-njfu/ChromosomeIdeogram/blob/main/example/satellite_DNA.png)
