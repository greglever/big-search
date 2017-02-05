# big-search
Case Study in Python searching datasets larger than RAM.

### Proposed Solution (Draft)

Generate a Burrows Wheeler Transform of the text T, `BWT(T)`, from the suffix array.

Take the `BWT(T)` and generate the FM index.

When searching through the FM index:

- If n+1 character is not in next index, increment mismatch count until the count is at maximum allowed.

In order to generate `BWT(T)`, read in reference genome sequence, output rotations into multiple files.

Sort multiple files lexicographically, using external merge sort techniques.

Then create FM index from these sorted files.
