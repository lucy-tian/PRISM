#!/bin/bash

get_repo_d () {
    local repo_dir="/home/lucytian/data/1_Single_Cell_PRS/0_trial"
    readlink -f "${repo_dir}"
}

source_paths () {
    # read paths / constants file in the specified directory and their parent directories
    local srcdir=$(readlink -f $1)
    local repo_dir=$(get_repo_d)

    # list of basenames
    paths_files=('paths.sh')

    while
        [ "${srcdir}" != '/' ] &&
        [ "${srcdir}" != $(dirname ${repo_dir}) ] ; do

        for paths_file in "${paths_files[@]}" ; do
            if [ -s "${srcdir}/${paths_file}" ] ; then
                source "${srcdir}/${paths_file}"
            fi
        done

        srcdir=$(dirname ${srcdir})
    done
}

get_slurm_log_d () {
    local srcdir=$(readlink -f $1)
    local log_d_root=$( readlink -f ${HOME}/slurm-logs )
    local repo_root=$(  readlink -f ${HOME}/repos )
    repo_path=$(echo ${srcdir} | sed -e "s%${repo_root}/%%")
    echo "${log_d_root}/repos/${repo_path}"
}

get_family () {
    local pheno=$1
    pheno_prefix=$(echo ${pheno} | sed -e 's/[0-9]//g')
    if [ ${pheno_prefix} == "INI" ] ; then
        echo "gaussian"
    else
        echo "binomial"
    fi
}
