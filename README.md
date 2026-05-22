# scutKinDev

`scutKinDev` is an R data package for working with *Scutigera*
kinematics from Physlets/Tracker `.trk` files, with a focus on
extracting top-view trajectories and generating exploratory plots from
tracked coordinates.

## Purpose

Tracker stores analysis data in `.trk` files, and Tracker projects
can bundle `.trk` files together with videos and supporting files in
`.trz` archives. This package is designed around a narrow
workflow: read Tracker-derived movement data for *Scutigera*,
standardize the coordinate output, and generate reproducible top-view
visualizations suitable for exploratory analysis and quality
control.

## Scope

The package is intended for datasets in which motion has been
digitized in Tracker using a calibrated coordinate system and exported
or saved in a Tracker-compatible format.

## Suggested workflow

1. Record or prepare a top-view video of the animal moving in a mostly planar arena.
2. In Tracker, import the video, define scale, and place the coordinate system so that the exported x-y coordinates have a biologically interpretable frame of reference.
3. Track the focal point or body landmark frame by frame, then save the analysis in a `.trk` file or a `.trz` project containing `.trk` content.
4. In R, use `scutKinDev` to ingest the file, tidy the coordinate time series, and produce top-view plots for path inspection, coverage, and movement summaries.

## Installation

If the package is hosted on GitHub, installation can be documented with `remotes`.

```r
# install.packages("remotes")
remotes::install_github("whrl/scutKinDev")
```

## TODO: Minimal example

