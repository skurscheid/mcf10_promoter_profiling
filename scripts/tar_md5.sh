#!/bin/bash

function tar_md5 {

DIR=${1}

if [[ -d $DIR ]]
  then 
  tar -cvpf ${DIR}.tar ${DIR}/| xargs -I '{}' sh -c "test -f '{}' &&  md5sum '{}'" | tee ${DIR}.md5
fi

}

export -f tar_md5