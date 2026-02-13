# GitHub Repository Setup Checklist

Your repository is now live at: **https://github.com/jilavsky/pyirena**

## ✅ Completed

- [x] Repository created on GitHub
- [x] Code pushed to main branch
- [x] All documentation files uploaded
- [x] GitHub Actions workflow configured

## 📋 Recommended GitHub Configuration

### 1. Update Repository Description and Topics

Go to: https://github.com/jilavsky/pyirena

Click **⚙️ Settings** (or edit repository details):

**Description:**
```
Python tools for small-angle scattering (SAXS/SANS/USAXS) data analysis, including the Unified Fit model
```

**Website:**
```
https://github.com/jilavsky/pyirena#readme
```

**Topics (tags for discoverability):**
```
saxs, sans, usaxs, small-angle-scattering, python, scientific-computing,
unified-fit, beaucage, x-ray-scattering, data-analysis, physics, materials-science
```

### 2. Enable GitHub Features

#### Enable Discussions (Optional)
- Go to Settings → General → Features
- Check ✅ "Discussions"
- This allows users to ask questions without creating issues

#### Enable GitHub Pages (for documentation)
- Go to Settings → Pages
- Source: Deploy from a branch
- Branch: `main` / `/docs`
- This will make your documentation available at: https://jilavsky.github.io/pyirena/

### 3. Configure Branch Protection (Recommended)

Go to: Settings → Branches → Add branch protection rule

**Branch name pattern:** `main`

Recommended settings:
- ✅ Require a pull request before merging
- ✅ Require status checks to pass before merging
  - Search for: `test` (from GitHub Actions)
- ✅ Require conversation resolution before merging

This ensures code quality by requiring reviews and passing tests.

### 4. Add Repository Secrets (for PyPI auto-upload)

When you're ready to publish to PyPI:

1. Create PyPI API token at: https://pypi.org/manage/account/token/
2. Go to: Settings → Secrets and variables → Actions
3. Click "New repository secret"
   - Name: `PYPI_API_TOKEN`
   - Value: `pypi-...` (your token)

This allows GitHub Actions to automatically publish releases to PyPI.

### 5. Create First Release (When Ready)

When you want to create v0.1.0:

1. Go to: Releases → Create a new release
2. Choose a tag: `v0.1.0` (create new tag)
3. Release title: `v0.1.0 - Initial Release`
4. Description: Copy from CHANGELOG.md
5. Check "Set as a pre-release" (since it's 0.1.0)
6. Publish release

### 6. Add Badges to README

Add these badges to the top of README.md for a professional look:

```markdown
[![GitHub release](https://img.shields.io/github/v/release/jilavsky/pyirena)](https://github.com/jilavsky/pyirena/releases)
[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/jilavsky/pyirena/workflows/Tests/badge.svg)](https://github.com/jilavsky/pyirena/actions)
[![PyPI](https://img.shields.io/pypi/v/pyirena)](https://pypi.org/project/pyirena/)
[![Downloads](https://pepy.tech/badge/pyirena)](https://pepy.tech/project/pyirena)
```

### 7. Configure Issue Labels

Go to: Issues → Labels

Add custom labels:
- `help-wanted` - Good for contributors
- `good-first-issue` - For newcomers
- `question` - User questions
- `priority-high` - Important issues
- `needs-testing` - Requires testing

### 8. Create a CODEOWNERS file (Optional)

Create `.github/CODEOWNERS`:

```
# These owners will be requested for review on PRs
*       @jilavsky
```

### 9. Add a Citation File

Create `CITATION.cff` for academic citations:

```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
  - family-names: "Ilavsky"
    given-names: "Jan"
    email: "ilavsky@aps.anl.gov"
title: "pyIrena: Python tools for small-angle scattering data analysis"
version: 0.1.0
date-released: 2024-02-13
url: "https://github.com/jilavsky/pyirena"
repository-code: "https://github.com/jilavsky/pyirena"
license: MIT
```

### 10. Monitor Your Repository

Set up notifications:
- Watch → Custom → Issues, Pull Requests, Releases

Regular tasks:
- Respond to issues within 48 hours
- Review pull requests promptly
- Keep dependencies updated
- Monitor GitHub Actions for failed builds

## 🔍 Verify Everything Works

### Check GitHub Actions
1. Go to: Actions tab
2. Verify workflows are enabled
3. Check if any workflows have run (they trigger on push)

### Test Clone and Install
```bash
# In a new directory
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
pip install -e .
python -c "from pyirena import UnifiedFitModel; print('Success!')"
```

### Verify Documentation
- Check README renders correctly on GitHub
- Verify all links work
- Check code blocks display properly

## 📢 Share Your Project

Once you're ready:

1. **Social Media:**
   - Twitter/X: Announce the release
   - LinkedIn: Share with professional network
   - ResearchGate: Post to relevant groups

2. **Scientific Community:**
   - APS mailing lists
   - Small-angle scattering forums
   - Relevant Slack/Discord communities

3. **Python Community:**
   - Reddit: r/Python, r/datascience, r/Physics
   - Python Weekly newsletter
   - Science Python newsletter

4. **Academic:**
   - Present at conferences (APS, ACS, etc.)
   - Write a paper about the software
   - Create tutorial videos

## 🎯 Next Development Steps

1. **Test installation** from GitHub
2. **Add more examples** with real data
3. **Create tutorial notebooks**
4. **Build documentation** with Sphinx
5. **Gather user feedback**
6. **Fix bugs and iterate**
7. **Publish to PyPI** when stable

## 📚 Useful GitHub Resources

- [GitHub Docs](https://docs.github.com/)
- [GitHub Actions](https://docs.github.com/en/actions)
- [GitHub Pages](https://pages.github.com/)
- [Open Source Guides](https://opensource.guide/)

## ✨ Your Repository is Live!

Visit: **https://github.com/jilavsky/pyirena**

Congratulations on publishing your open-source project! 🚀
