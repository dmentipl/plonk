---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

*Please keep the section headings but replace all the text with your own responses. And delete any optional sections that you don't need.*

### Describe the bug

A clear and concise description of what the bug is.

### To Reproduce

Ideally, please provide a [minimal, reproducible example](https://stackoverflow.com/help/minimal-reproducible-example).

Here is a contrived example.

```python
import plonk
snap = plonk.load_snap('path_to_file.h5')
snap['abc']
```

This produces the following error.

```python
---------------------------------------------------------------------------
ValueError                                Traceback (most recent call last)
<ipython-input-19-40e938fb10d8> in <module>
----> 1 snap['abc']

~/repos/plonk/plonk/snap/snap.py in __getitem__(self, inp)
   1060     ) -> Union[ndarray, SubSnap]:
   1061         """Return an array, or family, or subset."""
-> 1062         return self._getitem(inp, sinks=False)
   1063
   1064     def __setitem__(self, name: str, item: ndarray):

~/repos/plonk/plonk/snap/snap.py in _getitem(self, inp, sinks)
   1035         """Return an array, or family, or subset."""
   1036         if isinstance(inp, str):
-> 1037             return self._getitem_from_str(inp, sinks)
   1038         if sinks:
   1039             raise ValueError('Cannot return sinks as SubSnap')

~/repos/plonk/plonk/snap/snap.py in _getitem_from_str(self, inp, sinks)
   1028                 return self._get_array(inp_root).sum(axis=1)
   1029
-> 1030         raise ValueError('Cannot determine item to return.')
   1031
   1032     def _getitem(

ValueError: Cannot determine item to return.
```

The basic idea is to provide an easy way to reproduce the behaviour and see the output. E.g.

1. Load data '...'
2. Visualise using '...'
3. '....'
4. See error

### Expected behaviour

A clear and concise description of what you expected to happen.

### Screenshots (optional)

If applicable, add screenshots to help explain your problem.

### Environment (please complete the following information)

+ OS: [e.g. macOS 10.14.6, Ubuntu 18.04]
+ Python: [e.g. Python 3.7, Anaconda with conda 4.7.5, using Jupyter notebook]
+ Plonk version: [e.g. 0.3.0]

### Data source

If the dataset is small you can attach it to this issue. Alternatively, it can be helpful to provide a quick way to generate a similar dataset.

#### Example

Clone Phantom, checkout version `6666c55f`, and build with `SETUP=disc`. Then use the `mydisc.setup` file (attached) to generate the initial snapshot.

### Additional context (optional)

Add any other context about the problem here.
