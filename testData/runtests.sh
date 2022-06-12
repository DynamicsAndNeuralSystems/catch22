#!/bin/bash

"$(dirname -- $0)/../C/runAllTS.sh" -i "$(dirname -- $0)" -o "$(dirname -- $0)" -a "_output" -s 1 -i "$(dirname -- $0)" -o "$(dirname -- $0)" -a "_output" -s 1
