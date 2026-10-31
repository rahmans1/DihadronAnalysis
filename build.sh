#!/bin/bash
# Compile the macro once, before submitting. Jobs reuse the library.
set -euo pipefail
HERE=$(dirname "$(realpath -e "$0")")
source "$HERE/config.sh"
cd "$HERE"
singularity exec -B "$BINDS" "$IMAGE" root -l -b -q -e '.L ePIC_multikt2.C+'
