### 2nd MFT-MUON hackathon, 1 & 4 September 2020

## MC simulation

```bash
o2-sim -m MFT -e TGeant3 -n 5 -g boxgen --field=-5 --configKeyValues 'BoxGun.pdg=13; BoxGun.eta[0]=-3.8; BoxGun.eta[1]=-2.1; BoxGun.prange[0]=0.1; BoxGun.prange[1]=5.0; BoxGun.number=10'

o2-sim-digitizer-workflow -b 

o2-mft-reco-workflow -b
```

## The macros

```bash
root.exe -b -q Read_Kine.C+
```

Read_Kine = read the MC generated particles

Read_Hits = read the hits produced in the MFT detector

Read_digROF_1 = read the vector with digits and associate labels

Read_digROF_2 = read the digits inside each Read-Out-Frame (ROF)

Read_digROF_3 = find how the generated MC events give signals over the different Read-Out-Frames

(and similar with the macros Read_clsROF for the clusters)

Read_clsROF_1 = calculates the (x,y,z) coordinates of the cluster using the pattern ID and the pattern dictionary

Read_clsHits = associates MC hits to clusters using the MC labels

Read_Tracks = read reconstructed tracks and the clusters

Read_CompClsDict = read the cluster patterns stored in the dictionary file and prints the entries