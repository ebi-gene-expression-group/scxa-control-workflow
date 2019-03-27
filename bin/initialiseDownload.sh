confFile=$1
downloadConfig=$2
probeFile=$3
sshUser=${4:-''}

# Fetch params from global config

retries=$(parseNfConfig.py --paramFile $confFile --paramKeys params,downloadRetries)
fetchFreq=$(parseNfConfig.py --paramFile $confFile --paramKeys params,fetchFreqMillis)
allowedDownloadMethods=$(parseNfConfig.py --paramFile $confFile --paramKeys params,allowedDownloadMethods)
tempDir=$NXF_TEMP/atlas-fastq-provider

# Create downloader config

echo "ENA_RETRIES='$retries'" > $downloadConfig
echo "FETCH_FREQ_MILLIS='$fetchFreq'" >> $downloadConfig
echo "FASTQ_PROVIDER_TEMPDIR='$NXF_TEMP/atlas-fastq-provider'" >> $downloadConfig
echo "ALLOWED_DOWNLOAD_METHODS='$allowedDownloadMethods'" >> $downloadConfig

if [ -n "$sshUser" ]; then
    echo "ENA_SSH_USER='$sshUser'" >> $downloadConfig
fi

# Initialise the probe testing which download methods are viable

initialiseEnaProbe.sh -c $downloadConfig -t $probeFile
