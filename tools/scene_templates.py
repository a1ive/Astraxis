"""Expand body templates for offline tools. Property validation lives in scene_loader.cpp."""

import copy
from pathlib import Path, PurePosixPath
import tomllib


def _overlay(target, source):
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _overlay(target[key], value)
        else:
            target[key] = copy.deepcopy(value)


def load_scene(path):
    """Read an assets/scenes/*.toml file, applying single-level body templates."""
    path = Path(path)
    with path.open('rb') as file:
        scene = tomllib.load(file)
    for index, body in enumerate(scene.get('bodies', [])):
        if 'template' not in body:
            if 'unset' in body:
                raise ValueError(f'{path}: unset needs a body template')
            continue
        name = body.get('name', f'bodies[{index}]')
        template = body['template']
        context = f'{path}: body {name!r} template {template!r}'
        if (not isinstance(template, str) or not template or '\\' in template or ':' in template
                or PurePosixPath(template).is_absolute() or '..' in PurePosixPath(template).parts):
            raise ValueError(f"{context}: path must be relative to assets, using '/' and no '..'")
        try:
            with (path.parent.parent / template).open('rb') as file:
                merged = tomllib.load(file)
            if 'template' in merged or 'unset' in merged:
                raise ValueError('templates cannot reference other templates or use unset')
            unset = body.get('unset', [])
            if not isinstance(unset, list):
                raise ValueError('unset must be an array of top-level template field names')
            for key in unset:
                if not isinstance(key, str) or key not in merged:
                    raise ValueError(f'unset field {key!r} is absent from the template or repeated')
                del merged[key]
            _overlay(merged, {k: v for k, v in body.items() if k not in ('template', 'unset')})
        except (OSError, ValueError) as error:
            raise ValueError(f'{context}: {error}') from error
        scene['bodies'][index] = merged
    return scene
