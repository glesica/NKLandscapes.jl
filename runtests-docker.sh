#!/bin/sh

docker build -f Dockerfile-release -t julia-release .
docker build -f Dockerfile-latest -t julia-latest .

echo "Running tests on release version..."
docker run -v "$PWD":/opt/src julia-release
echo "Running tests on latest version..."
docker run -v "$PWD":/opt/src julia-latest
