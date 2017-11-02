# lblb-resonance

First need to put links to a cms3 version of CORE and a recent MT2Analysis in base directory:
```
ln -s <path to cms3 branch of CORE> CORE
ln -s <path to MT2Analysis directory> MT2Analysis
```

Everything else should work out of the box. Looper is based on MT2 code and works similarly. There are skimmed babies at

`/nfs-6/userdata/bemarsh/leptoquark/babies/v2`

with the skim

`nJet30>=2 && jet_pt[1]>=40 && nlep>=2 && lep_charge[0]*lep_charge[1]<0`

The script `plotting/make_plots.py` will make all of the plots in various regions, just like the CRPlotMaker in the MT2 repository.
