#!/bin/sh
# Requires a JDK (>=11). Output: build/ ; run with: java -cp build npctransit.Main <command>
mkdir -p build && javac -d build src/main/java/npctransit/*.java
