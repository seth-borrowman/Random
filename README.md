# Random
- qPCR_Template.qmd
  - Takes a .xls file exported from QuantStudio as input. Outputs a PDF plotting a standard curve with statistics and predicted copy \#'s for unknown samples as well as a .csv containing data for all replicates.
  - Download, install Quarto from here: https://quarto.org/docs/get-started/
  - In R Studio: Update header with title, name, etc. Update input_file with a .xls file path (eg. "D:/filename.xls"). Install any missing packages. Click "Render" (Ctrl+Shif+K on Windows).
- fastaRegexRemove.py
  - Can quickly trim sample names in a fasta by removing part of the name matching a regular expression.
  - Takes an input fasta (-i), a regular expression (-r), a string to replace the regex (-n, default=''), and a name for the output fasta (-o).
  - I recommend [this site](https://regex-generator.olafneumann.org/) to easily generate regex.
