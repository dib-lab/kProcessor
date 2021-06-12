#!/bin/bash

strace -ff -e trace=file "$@" 2>&1 | perl -ne 's/^[^"]+"(([^\\"]|\\[\\"nt])*)".*/$1/ && print'
