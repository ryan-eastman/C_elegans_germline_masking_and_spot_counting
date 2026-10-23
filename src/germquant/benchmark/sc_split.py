"""The registered SC validation split (tests/golden/sc_split.json; docs/SC_VALIDATION_PROTOCOL.md section 3).

``load_split()`` reads the byte copy shipped with the tests and checks its sha256 against the registration
(protocol changelog entry 1), ``role_of(stem_or_gonad)`` answers dev / test / reserved / excluded, and
``guard(...)`` raises ``SealedError`` for a test or reserved gonad unless ``final=True`` (only the scoring
script passes it, under ``--final`` with a reason). The sc-straighten command uses ``sealed_stems()`` to refuse
sealed images without ``--allow-sealed``.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path

REGISTERED_SHA256 = "db328c1cf4bc2b08dda32027a0e2567a128700a26facac4f7dc47738f67ae9eb"
SEALED_ROLES = ("test", "reserved")


class SealedError(PermissionError):
    """A test or reserved gonad (or a sealed seed) was asked for without ``final=True``."""

    __test__ = False


def default_path() -> Path | None:
    from germquant import provenance

    root = provenance._repo_root()
    if root is None:
        return None
    p = root / "tests" / "golden" / "sc_split.json"
    return p if p.is_file() else None


@lru_cache(maxsize=4)
def _load(path: str, verify: bool) -> dict:
    # the registered file has LF line endings; a CRLF checkout (core.autocrlf) is hashed after normalising
    data = Path(path).read_bytes().replace(b"\r\n", b"\n")
    if verify:
        h = hashlib.sha256(data).hexdigest()
        if h != REGISTERED_SHA256:
            raise ValueError(f"{path}: sha256 {h} is not the registered split {REGISTERED_SHA256} "
                             "(docs/SC_VALIDATION_PROTOCOL.md, changelog entry 1)")
    return json.loads(data.decode("utf-8"))


def load_split(path: str | Path | None = None, verify: bool = True) -> dict:
    p = Path(path) if path is not None else default_path()
    if p is None:
        raise FileNotFoundError("tests/golden/sc_split.json not found next to the package")
    return _load(str(p), bool(verify))


def _norm(s: str) -> str:
    return str(s).replace("\\", "/").rsplit("/", 1)[-1]


def gonad_record(key: str, split: dict | None = None) -> dict | None:
    """The split record of a gonad, by Madeleine's gonad id, by image stem or by RAD-51 tag."""
    split = split or load_split()
    k = _norm(key)
    for g in split["gonads"]:
        if k in (g.get("gonad"), g.get("image_stem"), g.get("rad51_tag")):
            return g
    return None


def role_of(key: str, split: dict | None = None) -> str | None:
    g = gonad_record(key, split)
    return None if g is None else g.get("role")


def sealed_stems(split: dict | None = None) -> set[str]:
    """Image stems of the test and reserved gonads while the split says ``sealed: true`` (else empty)."""
    split = split or load_split()
    if not split.get("sealed", False):
        return set()
    return {g["image_stem"] for g in split["gonads"] if g.get("role") in SEALED_ROLES and g.get("image_stem")}


def guard(keys, *, final: bool = False, split: dict | None = None) -> None:
    """Raise SealedError when any key is a test or reserved gonad and ``final`` is false."""
    split = split or load_split()
    bad = []
    for k in ([keys] if isinstance(keys, str) else list(keys)):
        r = role_of(k, split)
        if r in SEALED_ROLES:
            bad.append(f"{k} ({r})")
    if bad and not final:
        raise SealedError(
            f"refused: {', '.join(bad)} {'is' if len(bad) == 1 else 'are'} sealed by the SC validation protocol "
            "(docs/SC_VALIDATION_PROTOCOL.md section 3.4); test gonads are scored once, by the orchestrator, with "
            "scripts/benchmark_sc.py --final --reason TEXT")


def dev_gonads(split: dict | None = None, session: str | None = None) -> list[dict]:
    split = split or load_split()
    return [g for g in split["gonads"] if g.get("role") == "dev" and (session is None or g.get("session") == session)]
