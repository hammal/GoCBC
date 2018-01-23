# Matrix and Vector Extensions
This package includes convenient extensions to the methods provided in [Gonum](https://www.gonum.org)
many of which are inspired from the [Numpy](http://www.numpy.org) python package.

#### Matrix creating routines
Routines used for easily creating matrices
###### Eye
[Numpy Eye Interface](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.eye.html)

###### Ones
[Numpy Ones Interface](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.ones.html
###### Full
[Numpy Full Interface](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.full.html#numpy.full)

#### Matrix Vector Conversion
These are matrix vector conversion tools
This is redundant see the mat.VecDense receiver functions
- ColViewOf
- RowViewOf

And since the Vector interface includes the matrix one
also VecToMat is redundant as it is an matrix already.

###### ~~MatColToVec~~

###### ~~MatRowToVec~~

###### ~~VecToMat~~

#### Matrix checks norm etc
###### check if INF or NAN

## Todos
- TODO: Write tests for gonumExtensions
