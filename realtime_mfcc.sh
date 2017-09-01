#!/bin/sh
rec -t wav - 2>/dev/null | sox -t wav - -c 1 -r 8k -b 64 -e float -t raw - 2>/dev/null | ./mfcc
