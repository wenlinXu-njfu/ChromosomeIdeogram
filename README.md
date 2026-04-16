[![PyPI - Version](https://img.shields.io/pypi/v/chromosome-ideogram)](https://pypi.org/project/chromosome-ideogram)

# ChromosomeIdeogram
A python package for drawing genome-wide data on ideograms.

## Install
```shell
pip install chromosome-ideogram
```

## Usage example
```python
import matplotlib.pyplot as plt
from chromosome_ideogram import ChromosomeIdeogram

plt.rcParams['font.size'] = 8

fig = plt.figure(figsize=(10, 7))
axes = fig.add_subplot(111)

ci = ChromosomeIdeogram(
    karyotype='example/karyotype.xls',
    radius=0.4,  # Radius of chromosome (it is not recommended to make any changes)
    space=3  # Chromosome center coordinate interval (radius x space, also not recommended for modification)
)
ci.draw_chromh(axes=axes, centromere_angle=60, color='lightgrey')
ci.annotations(axes=axes, annotations='example/satellite_DNA.bed')

plt.savefig('example/satellite_DNA.pdf', bbox_inches='tight')
plt.savefig('example/satellite_DNA.png', bbox_inches='tight')
```
![image](https://github.com/wenlinXu-njfu/ChromosomeIdeogram/blob/main/example/satellite_DNA.png?raw=true)
