README:

Python3 compatible loading routines for LINCS (gctx files)

- Based on cmapPy, commit 7a2e18030, which requires Python2

Changes compared to cmapPy:
- Only contains loading routines relevant to GSE70138 (which also was used for testing), some other routines are commented out (as not tested).
- Introduces Python3 compatiblility by introducing a single line into parse_metadata_df [ temp_array = np.core.defchararray.decode(temp_array, 'utf8')  ]
- Comments out generation of GCToo objects, which would otherwise be the returned data, and instead directly provides the contained data frame (rationale: Thomas Stoeger never used GCToo objects directly; did not test or ensure Python3 compliance of GCToo objects)