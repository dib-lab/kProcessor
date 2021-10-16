# How to build Python docs?

```bash
pip install recommonmark
pip install nbsphinx
pip install sphinx_rtd_theme
mamba install pandoc
python -m sphinx -T -b html -d _build/doctrees -D language=en . _build/html
```
