# Distribution Guide for pyIrena

This guide explains how to distribute pyIrena on PyPI and Conda.

## Prerequisites

### For PyPI Distribution

```bash
pip install build twine
```

### For Conda Distribution

```bash
conda install conda-build anaconda-client
```

## Version Management

Before each release, update the version number in:

1. `pyproject.toml`:
   ```toml
   [project]
   version = "0.1.0"  # Update this
   ```

2. `pyirena/__init__.py`:
   ```python
   __version__ = "0.1.0"  # Update this
   ```

3. `conda/meta.yaml`:
   ```yaml
   {% set version = "0.1.0" %}  # Update this
   ```

4. `CHANGELOG.md`:
   ```markdown
   ## [0.1.0] - 2024-02-13
   ### Added
   - New feature description
   ```

## PyPI Distribution

### 1. Prepare for Release

```bash
# Ensure all tests pass
pytest

# Ensure code is formatted
black pyirena/

# Check for any issues
flake8 pyirena/
```

### 2. Build Distribution

```bash
# Clean previous builds
rm -rf dist/ build/ *.egg-info

# Build source and wheel distributions
python -m build
```

This creates:
- `dist/pyirena-0.1.0.tar.gz` (source distribution)
- `dist/pyirena-0.1.0-py3-none-any.whl` (wheel)

### 3. Check Distribution

```bash
# Check package metadata
twine check dist/*

# Test installation locally
pip install dist/pyirena-0.1.0-py3-none-any.whl
```

### 4. Upload to Test PyPI (Optional)

```bash
# Create account at https://test.pypi.org/

# Upload to test PyPI
twine upload --repository testpypi dist/*

# Test installation from test PyPI
pip install --index-url https://test.pypi.org/simple/ pyirena
```

### 5. Upload to PyPI

```bash
# Create account at https://pypi.org/
# Set up API token at https://pypi.org/manage/account/token/

# Upload to PyPI
twine upload dist/*

# Or with API token:
# TWINE_USERNAME=__token__ TWINE_PASSWORD=<your-token> twine upload dist/*
```

### 6. Verify Installation

```bash
pip install pyirena
python -c "import pyirena; print(pyirena.__version__)"
```

## Conda Distribution (conda-forge)

### Initial Setup (One-time)

1. Fork the conda-forge feedstock template
2. Create a new feedstock repository
3. Submit to conda-forge

### For Each Release

1. **Update recipe** in `conda/meta.yaml`:
   ```yaml
   {% set version = "0.1.0" %}  # New version
   {% set sha256 = "..." %}     # SHA256 of source tarball
   ```

2. **Get SHA256 hash** of PyPI source:
   ```bash
   # After uploading to PyPI
   wget https://pypi.io/packages/source/p/pyirena/pyirena-0.1.0.tar.gz
   sha256sum pyirena-0.1.0.tar.gz
   ```

3. **Build locally** (test):
   ```bash
   conda build conda/
   ```

4. **Submit to conda-forge**:
   - Update your feedstock fork
   - Create pull request to conda-forge

### Manual Conda Build (Alternative)

```bash
# Build package
conda build conda/

# Upload to anaconda.org
anaconda upload /path/to/pyirena-0.1.0.tar.bz2

# Or upload to your channel
anaconda upload -u YOUR_USERNAME /path/to/pyirena-0.1.0.tar.bz2
```

## GitHub Release

### Create Release Tag

```bash
# Tag the release
git tag -a v0.1.0 -m "Release version 0.1.0"
git push origin v0.1.0
```

### Create GitHub Release

1. Go to GitHub repository
2. Click "Releases" → "Draft a new release"
3. Choose tag `v0.1.0`
4. Fill in release notes from `CHANGELOG.md`
5. Attach distribution files (optional)
6. Publish release

## Pre-release Checklist

- [ ] All tests passing
- [ ] Code formatted with black
- [ ] Version numbers updated in all files
- [ ] CHANGELOG.md updated
- [ ] Documentation reviewed and updated
- [ ] Examples tested
- [ ] README.md accurate
- [ ] Dependencies list current

## Release Checklist

- [ ] Git tag created
- [ ] PyPI upload successful
- [ ] conda-forge updated (if applicable)
- [ ] GitHub release created
- [ ] Documentation deployed
- [ ] Announcement sent (if applicable)

## Common Issues

### Import errors after installation

Check that package name matches in all configuration files.

### Missing files in distribution

Update `MANIFEST.in` to include necessary files.

### Version mismatch

Ensure version is consistent in:
- `pyproject.toml`
- `pyirena/__init__.py`
- `conda/meta.yaml`
- Git tag

### Build fails on certain platforms

Test builds on:
- Linux (Ubuntu)
- macOS
- Windows

Use GitHub Actions for automated testing.

## Updating Package After Release

### For bug fixes (0.1.0 → 0.1.1)

```bash
# Update version numbers
# Update CHANGELOG.md
git commit -am "Bump version to 0.1.1"
git tag v0.1.1
git push origin main --tags
# Rebuild and upload
```

### For new features (0.1.0 → 0.2.0)

```bash
# Update version numbers
# Update CHANGELOG.md
git commit -am "Bump version to 0.2.0"
git tag v0.2.0
git push origin main --tags
# Rebuild and upload
```

### For breaking changes (0.1.0 → 1.0.0)

```bash
# Update version numbers
# Update CHANGELOG.md with migration guide
git commit -am "Bump version to 1.0.0"
git tag v1.0.0
git push origin main --tags
# Rebuild and upload
```

## Documentation Hosting

### GitHub Pages (Option 1)

Enable GitHub Pages in repository settings.

### Read the Docs (Option 2)

1. Sign up at https://readthedocs.org/
2. Import repository
3. Configure build settings

## Automation

Consider setting up:

- **GitHub Actions** for automatic PyPI upload on tag push
- **Dependabot** for dependency updates
- **Pre-commit hooks** for code quality

Example GitHub Action for PyPI upload:

```yaml
name: Publish to PyPI

on:
  release:
    types: [published]

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    - name: Build
      run: |
        pip install build
        python -m build
    - name: Publish
      uses: pypa/gh-action-pypi-publish@release/v1
      with:
        user: __token__
        password: ${{ secrets.PYPI_API_TOKEN }}
```

## Support

For issues with distribution:
- PyPI: https://pypi.org/help/
- conda-forge: https://conda-forge.org/docs/
- Packaging guide: https://packaging.python.org/
