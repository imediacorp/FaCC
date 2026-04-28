# Contributing to FaCC

Thank you for your interest in contributing to **Fibonacci as Cosmological Constant (FaCC)**.

## Getting started

1. **Fork** the repository and clone your fork locally.
2. Create a branch from `main`:
   ```bash
   git checkout -b feat/your-feature-name
   # or
   git checkout -b fix/short-description
   ```
3. Install the full development environment:
   ```bash
   python -m venv .venv && source .venv/bin/activate
   pip install -r requirements.txt -r requirements-dev.txt
   ```

## Code style

- Follow **PEP 8**. Run `flake8 src/ tests/` before committing — fatal errors (`E9`, `F63`, `F7`, `F82`) will block CI.
- Use **type hints** (PEP 484) on all public function signatures. The existing modules (`src/phi_modulation.py` etc.) show the expected style.
- Write **no comments unless the why is non-obvious**. Well-named functions and variables are self-documenting.
- Keep **docstrings short** — one-line summary plus Parameters / Returns blocks in NumPy style where needed.

## Testing

All PRs must include tests. Run the suite locally before pushing:

```bash
pytest tests/ -v --tb=short --cov=src --cov-report=term-missing
```

- Place tests in `tests/test_<module_name>.py`.
- Mock external calls (CAMB, network requests) using `unittest.mock.patch` — see `tests/test_phi_modulation.py` for CAMB mocking and `tests/test_desi_bao.py` for HTTP mocking.
- Aim to keep existing coverage; a PR that drops coverage by more than 2% will be asked to add tests.

## Pull requests

1. Keep PRs **focused** — one logical change per PR. Refactors and feature additions in separate PRs.
2. Update `CHANGELOG.md` under `[Unreleased]` with a one-line summary of your change.
3. Fill out the PR template — the checklist is there to help, not slow you down.
4. CI must be green before merge. All matrix builds (Python 3.9–3.12) must pass.

## Reporting bugs

Use the **Bug Report** issue template. Include:
- Python version and OS
- Full traceback
- Minimal reproducible example

## Suggesting features

Use the **Feature Request** issue template. Science-related feature requests should include a brief description of the physical motivation.

## Scientific contributions

If you are proposing a new analysis module (e.g., a new observable or dataset), please open a **Discussion** first to agree on the interface before writing code. New modules should:

- Live in `src/` with a `run_<name>_analysis()` entry point returning a plain dict.
- Export their entry point from `src/__init__.py`.
- Have a corresponding example script in `examples/`.
- Include unit tests that do not require network access or proprietary data files.

## Licence

By contributing you agree that your contributions will be licensed under the project's MIT licence.
