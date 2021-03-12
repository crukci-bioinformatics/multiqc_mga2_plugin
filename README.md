# [<img src="https://raw.githubusercontent.com/ewels/MultiQC/master/docs/images/MultiQC_logo.png" width="300" height="80" title="MultiQC">](https://multiqc.info/)

# Multi Genome Alignment 2 Plugin for MultiQC

For more information about MultiQC, see [https://multiqc.info](https://multiqc.info)

For more information about Multi Genome Alignment 2, see [https://github.com/crukci-bioinformatics/MGA2](https://github.com/crukci-bioinformatics/MGA2)

## Multi Genome Alignment 2 Plugin

The plugin will add a section for Multi Genome Alignment reports to the
MultiQC report. There is a plot of how reads are assigned to reference
species and tabular data for each data set.

### Requirements

The MGA2 plugin required Python 3.6 minimum and MultiQC version 1.9 or newer.

### Installing

The easiest way to install the plugin is to add it to your Python virtual
environment using _pip_.

```
python3 -m venv <path to virtual environment>
<path to virtual environment>/bin/activate
pip install git+https://github.com/crukci-bioinformatics/multiqc_mga2_plugin.git
```

The _pip_ command above will also install MultiQC if it is not already present.

### Installing for Development

When doing development work, the plugin will be cloned from GitHub onto your
local file system. You will need to create a virtual environment and add MultiQC
to it:

```
python3 -m venv <path to virtual environment>
<path to virtual environment>/bin/activate
pip install multiqc
```

Then, from the check out of the plugin, add a development installation of the
plugin to the environment:

```
python3 setup.py develop
```
