#!/bin/bash
#
set -euxo pipefail

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi

TAG="$1"

# cp ../scripts/*.java .
docker build --progress=plain -t fcunial/smaht_experiments .
docker tag fcunial/smaht_experiments fcunial/smaht_experiments:${TAG}
docker push fcunial/smaht_experiments:${TAG}
# rm -f *.java *.class
