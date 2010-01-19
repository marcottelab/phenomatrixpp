#!/bin/sh

aclocal -I ./aclocal
automake -a
autoconf
./configure --config-cache $*
