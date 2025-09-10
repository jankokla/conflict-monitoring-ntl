# Conflict Monitoring - Night Time Light

Semester project at the EPFL ECEO Lab on conflict monitoring using nighttime light data.

## Installation

Make sure you have **uv** installed on your system. If not, install it via:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then install project dependencies:

```bash
uv sync
```

## Documentation

Documentation dependencies are managed in the **docs** dependency group. You can serve the MkDocs site locally:

```bash
uv run --only-group docs mkdocs serve
```

This will start a local development server, usually available at:

```
http://127.0.0.1:8000/
```
