# [<img src="https://raw.githubusercontent.com/ewels/MultiQC/master/docs/images/MultiQC_logo.png" width="300" height="80" title="MultiQC">](https://multiqc.info/)

# Multi Genome Alignment Plugin for MultiQC

For more information about MultiQC, see [https://multiqc.info](https://multiqc.info)

For more information about Multi Genome Alignment, see [https://github.com/crukci-bioinformatics/MGA](https://github.com/crukci-bioinformatics/MGA)

## Multi Genome Alignment Plugin

The plugin will by default look for MGA summary files called "`*.mga.xml`". This is
not automatically how MGA names its files, but searching for "`*.xml`" is too greedy
and swallows files that might be read by other plugins.

This can be changed with a [MultiQC YAML config file](https://multiqc.info/docs/#configuring-multiqc)
setting the name of the files the MGA plugin will claim for itself:

```YAML
sp:
    mga:
        fn: '*.mga.xml'
```

Change the `fn` pattern to a glob that matches your MGA files.
