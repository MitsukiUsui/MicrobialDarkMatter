#!/bin/bash
set -ue

ftp_path=${1}
local_path=${2}

if [ ! -e ${local_path/.gz/} ]; then
    mkdir -p `dirname ${local_path}`
    wget --timeout 120 -O ${local_path} ${ftp_path}
    gunzip ${local_path}
fi
