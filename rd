#/bin/sh
image=rootproject/root-ubuntu16

if [ ! -z "$DISPLAY" ]
then
    sudo docker run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $PWD:/excercise --rm -it --user $(id -u) $image /bin/bash -c "cd /excercise/; $*"
else
    sudo docker run  -v $PWD:/home/builder  --rm -it --user $(id -u) $image /bin/bash -c "cd /home/builder; $*"
fi
