# Building cdx-reader

## Prerequisites

```bash
pip install build twine
```

## Build wheel

```bash
python -m build
```

This creates two files in `dist/`:
- `cdx_reader-0.1.0-py3-none-any.whl` (wheel)
- `cdx_reader-0.1.0.tar.gz` (source distribution)

## Install locally from wheel

```bash
pip install dist/cdx_reader-0.1.0-py3-none-any.whl
```

## Install in development mode

```bash
pip install -e .
```

## Run tests before release

```bash
pip install -e ".[dev]"
pytest
```

## Upload to PyPI

### Test PyPI (recommended first)

```bash
twine upload --repository testpypi dist/*
```

### Production PyPI

```bash
twine upload dist/*
```

## Clean build artifacts

```bash
rm -rf build/ dist/ *.egg-info cdx_reader/*.egg-info
```
