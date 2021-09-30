# Release history

## 0.2.1
- fix: enable to de-zoom of 10%
- fix: better display score for uniprot IDs
- fix(viewer.js): reduce length of reaction labels
- fix(viewer.js): reaction label is chosen only on python side
- feat: use template reaction ID as label
- feat: show reaction template ID in side panel

## 0.2.0
- chore: put cli instructions into methods
- fix: handle -1 global score as a specific case
- feat: generate sphinx documentation
- fix: also detect cofactor based on their IDs
- feat: updating js librairies
- fix: handle list of rule IDs
- chore(test.mk): remove error when no test collected
- fix: javascript variable scope
- feat: add uniprot IDs
- fix: rule score for pathways
- fix(utils): get dict values instead of the dict itself
- fix(utils): do not try to process empty inchi

## 0.1.3
- fix(utils.py): update according to new rpSBML interface
- fix(Viewer.py): update path to templates
- fix(viewer.js): update value type
- docs(utils): add type hints in methods
- build: rename conda dev env

## 0.1.2
- fix(cli): remove unused argument template folder
- fix(meta.yaml): prevent use of most recent rptools
- chore(setup.py): make use centralized metadata

## 0.1.1
- feat: accept both folder and tarfile as input
- tests: input tarfile example
- docs(README): update

## 0.1.0
- feat: update code according to the new rpSBML interface
- feat: update json specifications

## 0.0.1
- first working version