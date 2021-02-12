#!/bin/bash

source `dirname $0`/build-opm-flowdiagnostics-applications.sh

# Upstream revisions
declare -a upstreams
upstreams=(opm-flowdiagnostics
           opm-parser
           opm-material
           opm-core)

declare -A upstreamRev
upstreamRev[opm-flowdiagnostics]=master
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

build_opm_flowdiagnostics_applications
test $? -eq 0 || exit 1

cp serial/build-opm-flowdiagnostics-applications/testoutput.xml .
