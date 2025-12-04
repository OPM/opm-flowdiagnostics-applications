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

if grep -q "ert=" <<< $ghprbCommentBody
then
  ERT_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ert=([0-9]+).*/\1/g'`/merge
fi

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  OPM_COMMON_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

for upstream in $upstreams
do
  if grep -q "$upstream=" <<< $ghprbCommentBody
  then
    upstreamRev[$upstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*$upstream=([0-9]+).*/\1/g"`/merge
  fi
done

echo "Building with ert=$ERT_REVISION opm-common=$OPM_COMMON_REVISION opm-flowdiagnostics=${upstreamRev[opm-flowdiagnostics]} opm-parser=${upstreamRev[opm-parser]} opm-material=${upstreamRev[opm-material]} opm-core=${upstreamRev[opm-core]} opm-flowdiagnostics-applications=$sha1"

build_opm_flowdiagnostics_applications
test $? -eq 0 || exit 1

# If no downstream builds we are done
if ! grep -q "with downstreams" <<< $ghprbCommentBody
then
  cp serial/build-opm-flowdiagnostics-applications/testoutput.xml .
  exit 0
fi

echo "No registered downstreams."

test $? -eq 0 || exit 1
