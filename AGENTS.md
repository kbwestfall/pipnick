# pipnick — working notes for agents

## Conventions

- **Test data goes in `pipnick/data/tests/`, never `pipnick/tests/data/`.**
  `pipnick/tests/` holds test *code* only.
- **Never use NaN intentionally.** A value that is invalid, missing or to be
  ignored is flagged with a mask, never by setting it to NaN. A NaN always
  signifies a computational error. Do not write NaN as a placeholder, a fill
  value, an "undefined" result, or a gap in a plotted line.
- **Always use explicit condition checks.** Never rely on the truthiness or
  falsiness of a non-boolean object: write `if x is None`, `if len(x) == 0`,
  `if x != 0`, not `if not x` or `if x`. Truthiness on an actual `bool` is fine.
- **Every `assert` in a test carries a message describing the expected
  behaviour.** A bare assertion reports only that a comparison failed; the
  message is what tells the next reader what was supposed to happen.
- **Anything that wraps onto multiple lines** — a signature, a call, an `if`
  condition, a collection literal — puts its contents on the line *below* the
  opening construct, indented, with the closing bracket (and any trailing colon)
  alone on a line aligned with the line that opened it:

  ```python
  def build_camera_datatable(
      root, ext='.fits', path_key='path', camera=None, classify='header', skip_files=None
  ):

  tbl, rejected = build_camera_datatable(
      rawdir, ext=ext, path_key='rawdir', camera=camera, classify=classify
  )

  if (
      'max_temperature' in quality and temperature > quality['max_temperature']
  ):
  ```
- **Lines are at most 99 characters.** Enforced by `tox -e codestyle`, which
  configures `pycodestyle` from `tox.ini` to check `E9` (syntax/runtime errors)
  and `E501` (line length) and nothing else.
- **No nested function definitions.** A helper belongs at module level, private
  by a leading underscore if it is not part of the interface. The exception is
  test code, where a closure over the fixture under test is often the clearest
  way to express a check.
- **Keep module-level globals to a minimum.** Do not define a global constant
  unless (1) it is used in many different functions, or (2) it is in test
  code. Otherwise compute the value where it is used and explain it in that
  function's docstring.
- **Imports go at the top of the file.** The one exception is a script: the
  imports a `ScriptBase` subclass needs belong inside its `main`, because every
  registered entry point is imported when the console scripts are resolved, so
  a module-level import of the heavy machinery makes `--help` pay for it.
- **Avoid `import ... as`.** Industry-standard aliases (`import numpy as np`)
  are fine; otherwise import the module or name directly.
- **Any docstring containing a backslash is an r-string** (`r"""`), so that
  LaTeX in `:math:` roles reads as written rather than being double-escaped.
- Docstrings are **NumPy style**, rendered by `numpydoc`. References to objects
  in other packages take an **explicit role** so that intersphinx resolves them
  — `` :class:`numpy.ndarray` ``, `` :func:`astropy.stats.sigma_clip` `` — not
  the default role. Bare package names in prose are literals: ``` ``photutils`` ```.
- **The user makes every commit.** Leave changes in the working tree for
  review; never run `git commit` (or `git commit --amend`). The one exception
  is a rebase, when one is needed and agreed.
- **Two places record a change, and they do not overlap:**

  | Where | What goes there |
  |---|---|
  | Commit message | Very short — a single subject line. **No attribution or co-author trailers.** |
  | `docs/releases/` | Brief, simple statements of the work completed. No rationale. There is no root `CHANGES.rst`. |

- Configuration files are **TOML**, read with the standard-library `tomllib`.
- Adapted code carries an attribution note in its module docstring naming the
  source project and file.

## Environment

**Development is on macOS.** Shell commands must work with the BSD userland, not
the GNU one. The differences that actually bite:

| Tool | Watch for |
|---|---|
| `sed` | No `\b` word boundary, and `-i` requires an explicit backup suffix (`sed -i ''`). For anything beyond a literal substitution, use Python instead. |
| `grep` | Use `-E` for extended regular expressions; `\b` is unreliable. |
| `date`, `stat`, `readlink` | Different flags from the GNU versions. |

Python itself is installed via Homebrew.

Development install (Python ≥ 3.12), into the virtual environment `$ENVS_HOME/pipnick`:

```
python3 -m venv $ENVS_HOME/pipnick
source $ENVS_HOME/pipnick/bin/activate
pip install -e ".[dev]"
```

Virtual environments live in `$ENVS_HOME` (`~/.virtualenvs`), never in the
repository — an in-tree environment named `pipnick` would collide with the `pipnick/`
package directory.
