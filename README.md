# SIDM
Repo to start playing with SIDM analysis at coffea-casa. This may become the full-blown analysis repo, or I might wind up forking github.com/phylsix/Firefighter and/or github.com/phylsix/FireROOT

## Getting started
- Fork this repository ([here's a nice guide to follow](https://gist.github.com/Chaser324/ce0505fbed06b947d962))
- Log in to coffea.casa as described [here](https://coffea-casa.readthedocs.io/en/latest/cc_user.html#cms-authz-authentication-instance)
- Clone your fork of this repository using the [coffee-casa git interface](https://coffea-casa.readthedocs.io/en/latest/cc_user.html#using-git). Note the following:
  - Use the https link, not the ssh one (e.g. `https://github.com/btcardwell/SIDM.git`, not `git@github.com:btcardwell/SIDM.git`)
  - You will likely need to generate a [github personal access token](https://docs.github.com/en/enterprise-server@3.4/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) to log in to your github account on coffea-casa
- Navigate to the newly created SIDM directory using the coffea-casa file browser
- You should be good to go! If you want to test that your environment is set up correctly, try running any of the existing notebooks in SIDM/analysis/studies and comparing your output with the output shown [here](https://github.com/btcardwell/SIDM/tree/main/analysis/studies)

## How to update requirements.txt
```
cd SIDM/
pip install pipreqs
pipreqs . --force # overwrites current requirements.txt 
```
