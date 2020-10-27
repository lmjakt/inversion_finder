## get some sequence:

seq <- readLines( "example_seq.fa" )
tmp <- seq[1]
seq <- seq[2]
names(seq) <- tmp

dyn.load('src/inversion_finder.so')
sums <- .Call( "find_inversions", seq )

plot( sums[[1]], type='l' )
