# SIDM
Repo to start playing with SIDM analysis at coffea-casa. This may become the full-blown analysis repo, or I might wind up forking github.com/phylsix/Firefighter and/or github.com/phylsix/FireROOT

## Getting started
```
# log in to coffea.casa
# File >> New >> Terminal
git clone https://github.com/btcardwell/SIDM.git
cd SIDM
source setup.sh
```

## How to update requirements.txt
```
cd SIDM/
pip install pipreqs
pipreqs . --force # overwrites current requirements.txt 
```
