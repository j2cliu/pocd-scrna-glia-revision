# Frozen-display reproduction audit

- Status: **PASS**
- Final displays re-rendered: 7
- Submission TIFF hashes: 7/7 identical
- Convenience PNG preview hashes: 5/7 byte-identical
- Execution location: temporary directory removed after audit

The 300-dpi TIFF files are the journal-submission raster artifacts and are
the release-blocking comparison. PNG files are convenience previews; two
frozen previews were resized or re-exported after the canonical TIFF render
and are therefore reported but not used as a release criterion. PDF hashes
were not used because PDF device metadata can vary by render time.

This audit verifies deterministic reconstruction of the final submission
rasters from repository-contained derived inputs; it does not rerun the
upstream biological analyses or add independent validation.
