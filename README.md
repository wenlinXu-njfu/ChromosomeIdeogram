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

fig = plt.figure(figsize=(10, 12))
gs = gridspec.GridSpec(nrows=10, ncols=1)
axes1 = fig.add_subplot(gs[:8, :])
axes2 = fig.add_subplot(gs[8:, :])  # Used to display the detailed karyotype of a certain chromosome (such as chromosome 14)

ci1 = ChromosomeIdeogram(
    karyotype='example/karyotype.xls',
    radius=0.4,  # Radius of chromosome (it is not recommended to make any changes)
    space=4
)
ci1.draw_chromh(axes=axes1, centromere_angle=60, facecolor='lightgrey', edgecolor=None, linewidth=0)
ci1.annotations(axes=axes1, annotations='example/satellite_DNA.bed')

# Create a second karyotype diagram object (containing only chromosome 14)
ci2 = ChromosomeIdeogram(
    karyotype='example/Chr14_karyotype.xls',
    radius=0.4,
    space=4
)
ci2.draw_chromh(axes=axes2, centromere_angle=60, facecolor='white', edgecolor='#4d4d4d', linewidth=1)
ci2.annotations(axes=axes2, annotations='example/Chr14_satellite_DNA.bed')
ci2.draw_synteny(axes=axes2, synteny_file='example/Chr14_SV.bed') # show synteny

axes1.legend(loc=(0.85, 0.85))  # Set axes1 legend, you can also set axes2 legend

plt.subplots_adjust(hspace=6)  # Adjust the vertical distance between the two subgraphs

plt.savefig('example/satellite_DNA.pdf', bbox_inches='tight')
plt.savefig('example/satellite_DNA.png', bbox_inches='tight')
```
![image](https://github.com/wenlinXu-njfu/ChromosomeIdeogram/blob/main/example/satellite_DNA.png?raw=true)
