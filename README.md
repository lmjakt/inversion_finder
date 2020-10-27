# inversion_finder

Find inverted repeats in long sequences quickly.

inversion_finder provides a function to R that performs a cancelling
k-mer count, decrementin an overall k-mer number when encountering
k-mers whose reverse complement have been observed previously in
the sequence analysed.

## Usage

To compile:

```sh
cd src
R CMD SHLIB inversion_finder.cpp
```

Then follow the example in the `example.R` file.

Note that the function provided by `inversion_finder.so` is likely to
be unsafe; it may crash your R-session if you do not provide it with
the required arguments. This should be fixed in the future.