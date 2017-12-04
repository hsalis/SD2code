CONTAINER_IMAGE="index.docker.io/hsalis/sd2docker:0.11"

. _util/container_exec.sh

COMMAND="python searchReads.py"
PARAMS="${pathToDirectory}"

container_exec ${CONTAINER_IMAGE} ${COMMAND} ${PARAMS}
