# MFT-MUON-Hackathon

## Run AliRoot inside a docker container

1) Pull the docker image from the docker hub

```bash
docker pull  bovulpes/mft-muon-hackathon:v01
```

2) Define the working directory 

```bash
export ALIROOT_HOME=/<where-you-want>/
```

3) start a container named "aliroot" from the downloaded image

```bash
docker run -itd --name aliroot -v $ALIROOT_HOME:$ALIROOT_HOME --net=host bovulpes/mft-muon-hackathon:v01 /bin/bash
```

4) add your local user credentials to the running container, in order to be
able to share the files created inside the container with your local user
outside the container

```bash
docker exec aliroot useradd -u $UID $USER
```

5) open an interactive  shell session inside the container

```bash
docker exec -it --user $USER -w $PWD -e DISPLAY -e XAUTHORITY aliroot /bin/bash
```

6) define also inside the container the working directory from the point 2)

```bash
export ALIROOT_HOME=/<same-as-at-2>/
```

7) create some links to prepare to use the already installed AliRoot
environment

```bash
cd $ALIROOT_HOME

mkdir -p root/v5-34-30

mkdir -p geant3/v2-0

mkdir -p aliroot/feature-itsmft

ln -s /opt/alice/ROOT/v5-34-30 root/v5-34-30/inst

ln -s /opt/alice/GEANT3/v2-0_v5-34-30 geant3/v2-0/inst

ln -s /opt/alice/AliRoot-feature-itsmft aliroot/feature-itsmft/inst

cp /opt/alice/alice-env.sh .

cp /opt/alice/alice-env.conf .

source alice-env.sh -n 9
```

if everything goes fine, the output should finish with these lines:

```bash
  AliEn          <not found>

  ROOT           /home/vulpescu/alice/aliroot/opt/root/v5-34-30/inst

  Geant3         /home/vulpescu/alice/aliroot/opt/geant3/v2-0/inst

  AliRoot Core   /home/vulpescu/alice/aliroot/opt/aliroot/feature-itsmft/inst
```

8) repeat 5) if you want to open other terminals inside the same container; at
the end, when you don't need it anymore, stop end remove it

```bash
docker stop aliroot

docker rm aliroot

```

## Run a simulation + reconstruction

```bash
. runMFT.sh <n_mult> <n_events> 1
```