#!/bin/bash
set -euo pipefail

input_file=$1

while IFS= read -r cmd; do
  eval "$cmd"
done < "$input_file"
