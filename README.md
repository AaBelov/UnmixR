### unmixR
#### An R package for unmixing of hyperspectral images

`unmixR is a WORK IN PROGRESS`.  The fundamental structures & behavior are changing frequently.

`unmixR` is supported by [Google Summer of Code](http://www.google-melange.com/gsoc/homepage/google/gsoc2013).  Thank you!

Note: hyperspectral data are also called 'imaging spectroscopy' and 'imaging spectrometer data' depending upon the discipline.  Such data consists of spectra collected over an x, y grid.  Data sets like this are found in airborne land imaging studies, biomedical studies and art history investigations.  The spectra are often visible, infrared, near-infrared, raman spectra or mass spectrometer data sets.

#### Things to do and Things to remember + Misc Notes

An informal list to keep us on track.  Naturally, edit as you see fit.

##### Top priority
* nfindr variants
    * nfindr99: In working order, could be optimized
    * nfindrLDU
    * nfindrSeqLDU
* predict.nfindr
* class definition(s)
* changes needed for build & check
* vignettes

##### Lower priority

* methods for finding the number of endmembers (recommended in Harsayni2012)
    * hfc (prelim version seems OK)
    * nwhfc
    * sml (second moment linear)
    * Hysime (hyperspectral signal subspace identification by minimum error)

##### Misc notes

* A possible name for the package is `unmixR`

