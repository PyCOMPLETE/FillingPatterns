# Tool for the analysis of LHC filling schemes

**Authors:** G. Iadarola, A. Poyet, G. Sterbini

## To install:
```bash
git clone https://github.com/giadarol/fillingpatterns.git
pip install ./fillingpatterns
```

## To run an example:
```bash
cd fillingpatterns/examples
python 001_from_csv_analyze_bb.py
```

## Usage
The filling scheme can be loaded in different ways:
 * From a json file (as provided by the LPC filling scheme tool):
```python
import fillingschemes as fp
fp.FillingPattern.from_json('fname.json')
```
 * From a csv file (which can be benerated by this tool):
```python
import fillingschemes as fp
fp.FillingPattern.from_csv('fname.csv')
```

 * By providing two boolian arrays with the scheme:
```python
import fillingschemes as fp
fp.FillingPattern(pattern_b1, pattern_b2)
```
