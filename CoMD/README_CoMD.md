
The [CoMD](http://www.exmatex.org/comd.html) miniapp performs molecular dynamics using a cell method
for improved scaling of the force loop.

Documentation and comments are pretty sparse.  See the original CoMD code for more documentation.

So far there are three variants, ordered in increasing complexity, to make it easier to port to
new languages.

* nsquared  This variant uses the basic N^2 algorithm for particle interactions.    It should be the easiest to port to a new language.
* cell This variant uses the cell method.
* mpi Builds on the cell variant and runs in parallel using MPI.





