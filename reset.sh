#!/usr/bin/env bash

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "! WARNING: This will delete all output and logs !"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo ""
echo "Do you want to continue? y/[n]"

read response

if (( $response != "y" )); then
	exit
fi

git pull

rm *.log
rm -r output/
rm data/*.ENCFF*.fa

