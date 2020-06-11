
die (){
    errMsg=$1
    errCode=$2

    if [ -n "$errCode" ]; then
        errMsg="$errMsg - exiting with error code $errCode"
    else
        errCode=1
    fi

    echo -e "$errMsg" 1>&2
    exit $errCode
}

# Check the difference between config files or derived SDRFs, ignoring headers

diff_configs (){
    file1=$1
    file2=$2

    diff -I '^//' $file1 $file2 > /dev/null 2>&1
    return $? 
}
