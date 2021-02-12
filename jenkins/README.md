# opm-flowdiagnostics-applications jenkins build scripts:

**build-opm-flowdiagnostics-applications.sh**:
This is a helper script which contains a function for building,
testing and cloning opm-flowdiagnostics-applications and its dependencies.

**build.sh**:
This script will build dependencies, then build opm-flowdiagnostics-applications
and execute its tests.
It is intended for post-merge builds of the master branch.

**build-pr.sh**:
This script will build dependencies, then build opm-flowdiagnostics-applications
and execute its tests. It inspects the $ghbPrBuildComment environmental
variable to obtain a pull request to use for ert, opm-common,
opm-flowdiagnostics, opm-parser, opm-material and opm-core (defaults to master)
and then builds $sha1 of opm-flowdiagnostics-applications.

It is intended for pre-merge builds of pull requests.

You specify a given pull request to use for ert, opm-common,
opm-flowdiagnostics, opm-parser, opm-material and opm-core through the trigger.
The trigger line needs to contain &lt;modulename&gt;=&lt;pull request number&gt;.
