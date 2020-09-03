---
title: "PO activity calculator"
---

![](logo.png)

> Automated phenoloxidase activity calculation from absorbance assay data

This shiny app takes absorbance assay data for 96-well plate runs in an xlsx file and protein mass and sample IDs in a csv file and automates calculation of the maximum velocity and phenoloxidase activity standardised by protein mass and time (in minutes).  

Outputs a csv file with results and an html report containing a plot of absorbance vs time for each reaction, showing any discarded datapoints, the fitted curve and the timepoint used for the estimated max. velocity.  

The deployed app is available [here](https://sjp-analytics.shinyapps.io/POactivity/).  

## Developing

To modify this app to work with similar assays/data:

```shell
git clone https://github.com/lizardburns/po_activity_app
```

## Features

This app analyses absorbance data in order to calculate phenoloxidase activity per mg of protein per minute.

* It takes two files supplied by the user:
  - an absorbance assay data file (development to date limited to files generated by a Tecan 96-well plate reader) in xlsx format
  - a csv file with 3 fields ("well", "sampleID", "mg_protein") which is joined with velocities calculated from absorbance data using well ids (e.g. "A1", "A2", ..., "H12")
* Accomodates absorbance data in two formats (see example_data):
  - datatables for 475nm, 600nm and 490nm arranged sequentially in first worksheet
  - 490nm readings only where datatable is transposed
* Finds data in worksheet based on presence of "Cycle Nr." variable. Looks for string "cycle" (case insensitive)
* Calculates "Vmax" for each reaction as follows:
  - Scans time series for two consecutive increases in absorbance (i.e. 3 datapoints enclosing positive changes in absorbance). Throws out data preceding such a consistent ascent to facilitate fitting of a smoothed line (loess fitting). Ceases scan if condition not met before approx. two-thirds of total assay time has elapsed, to avoid use of very late activity unlikely to represent true enzyme activity.
  - Generates `NA` result if no suitable datapoints to fit curve.
  - Performs local polynomial regression fitting (`loess` function with span argument fixed at 0.6) to estimate smooth fitted line and approximates derivative to get max change in absorbance per second.
* Converts max velocity to change in absorbance per minute (Dabs / min, for each wavelength processed on run). Generates warning if either sample replicate is twice the other prior to calculating mean.
* In the case of all three wavelengths on a run, calculates difference at 475nm and 600nm wavelengths and performs remaining operations on these differences as well as the change in absorbance per minute at 490nm.
* Checks for presence of "blank" (case insensitive) among sampleIDs and prompts user for blank values if not found.
* Deducts blank value and standardises per milligram of protein (using user supplied protein mass).
* Makes two files available for download:
  - csv file with results
  - html report (rendered rmarkdown document) with a plot of absorbance against time for each reaction, which should be scanned by users to verify accuracy of curve-fitting and identify any reactions where it may be better to calculate the maximum velocity manually. Fitted line shown for span of timepoints used in estimation and vertical line on plot marks timepoint where maximum velocity recorded.
* example input data files are available to present expected input file format.

## Contributing

If you'd like to contribute, please fork the repository and use a feature branch. Pull requests are warmly welcome.  

## Links

- Project homepage: https://github.com/lizardburns/po_activity_app
- Related projects:
  - Seeking Survivors: https://www.seekingsurvivors.org

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/) (CC-BY-NC).  
