# Contributing to pyIrena

Thank you for your interest in contributing to pyIrena! This document provides guidelines and instructions for contributing.

## Code of Conduct

Please be respectful and constructive in all interactions. We aim to create a welcoming environment for all contributors.

## How to Contribute

### Reporting Bugs

If you find a bug:

1. Check if the issue already exists in [GitHub Issues](https://github.com/jilavsky/pyirena/issues)
2. If not, create a new issue with:
   - Clear title and description
   - Steps to reproduce
   - Expected vs actual behavior
   - Python version and OS
   - Minimal code example

### Suggesting Features

Feature requests are welcome! Please:

1. Check existing issues and discussions
2. Create a new issue describing:
   - The feature and its benefits
   - Possible implementation approach
   - Any relevant examples or references

### Contributing Code

#### Setup Development Environment

1. Fork the repository on GitHub
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/pyirena.git
   cd pyirena
   ```

3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

4. Install in development mode:
   ```bash
   pip install -e ".[dev]"
   ```

#### Making Changes

1. Create a new branch:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/your-bugfix-name
   ```

2. Make your changes following our coding standards (see below)

3. Add or update tests as needed

4. Run tests:
   ```bash
   pytest
   ```

5. Format code:
   ```bash
   black pyirena/
   ```

6. Check linting:
   ```bash
   flake8 pyirena/
   ```

#### Submitting Changes

1. Commit your changes:
   ```bash
   git add .
   git commit -m "Brief description of changes"
   ```

2. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

3. Open a Pull Request on GitHub with:
   - Clear description of changes
   - Reference to related issues
   - Screenshots if applicable

## Coding Standards

### Python Style

- Follow [PEP 8](https://pep8.org/) style guide
- Use [Black](https://black.readthedocs.io/) for code formatting (line length: 100)
- Use meaningful variable and function names
- Add docstrings to all public functions and classes

### Documentation

- Use Google-style docstrings:
  ```python
  def my_function(param1: int, param2: str) -> bool:
      """
      Brief description.

      Longer description if needed.

      Args:
          param1: Description of param1
          param2: Description of param2

      Returns:
          Description of return value

      Raises:
          ValueError: When input is invalid
      """
      pass
  ```

### Testing

- Write tests for new features
- Maintain or improve code coverage
- Use pytest framework
- Place tests in `pyirena/tests/`
- Name test files `test_*.py`

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb (Add, Fix, Update, Remove, etc.)
- Reference issue numbers when applicable

Example:
```
Add support for custom correlation functions (#42)

- Implement generic correlation function interface
- Add examples for common correlation types
- Update documentation
```

## Development Workflow

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pyirena

# Run specific test file
pytest pyirena/tests/test_unified.py

# Run with verbose output
pytest -v
```

### Building Documentation

(When documentation is set up)

```bash
cd docs
make html
```

### Type Checking

```bash
mypy pyirena/
```

## Areas for Contribution

We especially welcome contributions in these areas:

1. **Additional Analysis Tools**
   - Size distribution calculations
   - Form factor models
   - Structure factor models

2. **Improved Slit Smearing**
   - Full Lake integration implementation
   - Resolution smearing

3. **GUI Development**
   - Interactive fitting interface
   - Parameter adjustment tools

4. **Performance Optimization**
   - Parallel fitting for batch processing
   - Optimization of calculation routines

5. **Documentation**
   - Tutorial notebooks
   - Theory documentation
   - API documentation improvements

6. **Testing**
   - Additional test cases
   - Integration tests
   - Benchmark tests

## Questions?

Feel free to:
- Open an issue for questions
- Start a discussion on GitHub
- Contact the maintainer: ilavsky@aps.anl.gov

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## Acknowledgments

All contributors will be acknowledged in the project documentation and release notes.

Thank you for contributing to pyIrena!
