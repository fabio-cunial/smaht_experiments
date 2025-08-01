#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l MergeAndSubsample.wdl
java -jar ${WOMTOOL_PATH} validate -l DownloadBams.wdl
