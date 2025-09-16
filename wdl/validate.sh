#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l ClippedAlignments.wdl
java -jar ${WOMTOOL_PATH} validate -l SortByReadId.wdl
java -jar ${WOMTOOL_PATH} validate -l Merge.wdl
java -jar ${WOMTOOL_PATH} validate -l AbPoa.wdl
java -jar ${WOMTOOL_PATH} validate -l Gridss.wdl
java -jar ${WOMTOOL_PATH} validate -l Manta.wdl
java -jar ${WOMTOOL_PATH} validate -l Sniffles.wdl
java -jar ${WOMTOOL_PATH} validate -l MergeAndSubsample.wdl
java -jar ${WOMTOOL_PATH} validate -l DownloadBams.wdl
