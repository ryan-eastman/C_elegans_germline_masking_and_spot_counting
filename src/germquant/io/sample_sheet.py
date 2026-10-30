"""The sample sheet: a study's design columns, image by image (docs/PROJECTS.md, "Sample sheet";
docs/DESIGN_system_packs.md section 7.4; PLAN.md decision D21).

A sample sheet is a CSV file. Its first line names the columns:

* the key, ``image_id`` (the file name without its format suffix, io.formats.image_stem, which names the results
  folder) or ``file`` (a file name, whose image id is its image_stem, or a path, relative to the sheet's folder, with
  / or \\ between folders on every system, which then names that one file);
* any of the design columns ``genotype``, ``sex``, ``germ_cell``, ``treatment``, ``replicate`` and
  ``acquisition_date``: process_image(sample=row) writes them into every per-image table instead of what the file
  name gives (parse_sample with metadata.filename_regex); a cell left empty keeps the file name's value, and when a
  row gives sex but no germ_cell, germ_cell follows the row's sex by parse_sample's rule (herm: oocyte, male:
  spermatocyte);
* extra columns of the lab's own, which no run reads: batch.collect(sample_sheet=) and `germquant collect --samples`
  join them into every batch_<table>.csv.

Values are taken as given, with spaces around them stripped. A non-ISO acquisition_date (anything but YYYY-MM-DD or
YYYYMMDD) is kept, with a SampleSheetWarning: Excel rewrites dates in the computer's own format when it saves. So
are a germ_cell other than oocyte or spermatocyte and, on a row without germ_cell, a sex other than herm or male (H,
HERM, M or MALE in any letter case), each with a SampleSheetWarning: germ_cell then gives the image no expected
synaptonemal complex or crossover count (io.sample_metadata.expected_sc_count), which the C. elegans recipes take
from it.

The file is read as UTF-8 (a byte-order mark is allowed), else as cp1252 with a warning (Excel's plain "CSV" on
Windows), so either save reads the same letters. Line ends, the byte-order mark, the encoding, the column order,
quoting and the padding of cells never change what is read. write_sample_sheet writes ASCII with LF line ends,
atomically through fsutil.long_path, and refuses a value outside ASCII, naming its cell.

The sheet is not hashed, like the file name it replaces. A row's identity is its row_sha256: the sha256 of the
canonical JSON of the row's design cells (keys sorted, ASCII escapes, no spaces; PLAN.md 9.3), so a corrected
genotype reruns that image only, and an edited extra column, which no run reads, reruns nothing.
process_image records it in run_geometry (sample_sheet {source, row_sha256}) and in the done marker
(sample_row_sha256), in sheet runs only, and batch.is_done compares it when given. An image the sheet has no row
for keeps the file name's design columns and raises the flag metadata:not_in_sample_sheet, in sheet runs only.

A time series takes no sheet row in this release (SERIES_REFUSAL): the time-lapse loop passes none to its
timepoints. The module imports only the standard library and germquant's light helpers.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import time
import uuid
import warnings
from dataclasses import dataclass
from datetime import date
from pathlib import Path, PurePosixPath, PureWindowsPath

from ..fsutil import long_path
from ..schema import SHARED_META
from .formats import image_stem

DESIGN_COLUMNS = ("genotype", "sex", "germ_cell", "treatment", "replicate", "acquisition_date")
KEY_COLUMNS = ("image_id", "file")
NOT_IN_SHEET_FLAG = "metadata:not_in_sample_sheet"
SERIES_REFUSAL = ("a sample sheet does not apply to a time series in this release: the time-lapse loop gives its "
                  "timepoints no sample row. Run the series without --samples (its design columns then come from "
                  "the file name, through metadata.filename_regex and metadata.defaults), or run its timepoints as "
                  "separate images (docs/PROJECTS.md, \"Sample sheet\")")
ENCODINGS = ("utf-8-sig", "cp1252")
RETRIES = 5                      # after the first attempt: five more over RETRIES * RETRY_WAIT_S = 2 s
RETRY_WAIT_S = 0.4
_ISO_DATE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})$|^(\d{4})(\d{2})(\d{2})$")
_GERM_CELLS = ("oocyte", "spermatocyte")       # the words expected_sc_count counts (io/sample_metadata.py)
_NO_VALUE = ("unknown", "na")                  # a design cell that says the value is not known (any letter case)

# seams the tests replace (a failing or a briefly locked file system)
_replace = os.replace
_sleep = time.sleep


class SampleSheetError(ValueError):
    """A sample sheet that cannot be read or written; the message names the file and the cell or row."""


class SampleSheetWarning(UserWarning):
    """Something about a sample sheet that is read anyway (a cp1252 file, a non-ISO date, a sex or germ_cell no
    expected count can be made from)."""


@dataclass(frozen=True)
class SampleRow:
    """One image's row of a sample sheet, as process_image(sample=) takes it. ``fields``: the design cells the row
    sets (stripped, empty cells left out); ``source``: the sheet's path; ``row_sha256``: row_sha256(fields), or None
    when the sheet has no row for the image (then ``fields`` is empty and the run raises metadata:not_in_sample_sheet).
    """
    fields: dict
    source: str
    row_sha256: str | None


def row_sha256(fields: dict) -> str:
    """The sha256 of the canonical JSON of a row's design cells: json.dumps(fields, sort_keys=True,
    ensure_ascii=True, separators=(",", ":")) encoded as ASCII (PLAN.md 9.3, D21)."""
    text = json.dumps({str(k): str(v).strip() for k, v in fields.items()}, sort_keys=True, ensure_ascii=True,
                      separators=(",", ":"))
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def sample_row(fields: dict, source: str = "api") -> SampleRow:
    """A SampleRow made in code (a notebook, a test): ``fields`` are design cells, stripped, empty ones left out."""
    clean = {k: str(v).strip() for k, v in dict(fields).items() if v is not None and str(v).strip()}
    row = SampleRow(clean, str(source), row_sha256(clean))
    check_row(row)
    return row


def check_row(row) -> None:
    """Refuse anything process_image(sample=) cannot take: TypeError for an object that is not a SampleRow, and
    SampleSheetError for a field that is not a design column, a value that is not text, a row_sha256 that does not
    match the fields, or fields on a row that is not in the sheet."""
    if not isinstance(row, SampleRow):
        raise TypeError(f"sample must be a germquant.io.sample_sheet.SampleRow (read_sample_sheet(path).row(image_id, "
                        f"path), or sample_row(fields)), not {type(row).__name__}")
    if not isinstance(row.fields, dict):
        raise SampleSheetError(f"a sample row's fields must be a dict of design columns, not "
                               f"{type(row.fields).__name__}")
    unknown = [k for k in row.fields if k not in DESIGN_COLUMNS]
    if unknown:
        raise SampleSheetError(f"a sample row sets only the design columns ({', '.join(DESIGN_COLUMNS)}); "
                               f"{', '.join(map(repr, unknown))} is not one of them")
    if any(not isinstance(v, str) for v in row.fields.values()):
        raise SampleSheetError("a sample row's values are text (str), as a sheet's cells are")
    if row.row_sha256 is None:
        if row.fields:
            raise SampleSheetError("a sample row that is not in the sheet (row_sha256 None) has no fields")
        return
    if row.row_sha256 != row_sha256(row.fields):
        raise SampleSheetError(f"the sample row's row_sha256 {row.row_sha256} does not match its fields (made from "
                               f"them: {row_sha256(row.fields)}); make rows with read_sample_sheet or sample_row")


def germ_cell_of(sex: str) -> str:
    """parse_sample's rule (io/sample_metadata.py): herm (H, HERM in any case) is oocyte, male (M, MALE) is
    spermatocyte, anything else unknown."""
    from .sample_metadata import _SEX_CANON, _SEX_TO_GERMCELL

    s = str(sex).strip().lower()
    return _SEX_TO_GERMCELL.get(_SEX_CANON.get(s, s), "unknown")


def apply_row(sample: dict, row: SampleRow) -> list[str]:
    """Write a row's design cells over parse_sample's fields (``sample``, changed in place), germ_cell following the
    row's sex when the row leaves it out. The flags to raise: metadata:not_in_sample_sheet for an image the sheet has
    no row for, whose fields stay the file name's; none otherwise."""
    if row.row_sha256 is None:
        return [NOT_IN_SHEET_FLAG]
    for key in DESIGN_COLUMNS:
        if key in row.fields:
            sample[key] = row.fields[key]
    if "sex" in row.fields and "germ_cell" not in row.fields:
        sample["germ_cell"] = germ_cell_of(row.fields["sex"])
    return []


def iso_date(value: str) -> bool:
    """True for an ISO 8601 calendar date, YYYY-MM-DD or YYYYMMDD, that exists."""
    m = _ISO_DATE.match(str(value).strip())
    if not m:
        return False
    y, mo, d = (int(g) for g in (m.groups()[:3] if m.group(1) else m.groups()[3:]))
    try:
        date(y, mo, d)
    except ValueError:
        return False
    return True


def _germ_cell_warnings(where: str, n: int, image_id: str, fields: dict) -> list[str]:
    """The warning for a row whose germ cell gives no expected count (expected_sc_count): a germ_cell other than
    oocyte or spermatocyte, or, on a row without germ_cell, a sex that parse_sample's map does not know (germ_cell
    then reads unknown). unknown and NA say the value is not known, and are not warned about. The value is kept as
    given (D21)."""
    if "germ_cell" in fields:
        value = fields["germ_cell"]
        if value in _GERM_CELLS or value.lower() in _NO_VALUE:
            return []
        return [f"{where}: row {n} (image {image_id}): germ_cell {value!r} is not oocyte or spermatocyte (in lower "
                f"case), so a recipe that counts synaptonemal complexes or crossovers gives the image no expected "
                f"count; it is kept as given"]
    sex = fields.get("sex")
    if sex is None or sex.lower() in _NO_VALUE or germ_cell_of(sex) != "unknown":
        return []
    return [f"{where}: row {n} (image {image_id}): sex {sex!r} is not herm or male (H, HERM, M or MALE in any letter "
            f"case), so germ_cell reads unknown and a recipe that counts synaptonemal complexes or crossovers gives the "
            f"image no expected count; for C. elegans write herm or male, or give a germ_cell column (oocyte or "
            f"spermatocyte)"]


@dataclass(frozen=True)
class _Row:
    line: int                    # the row's number in the file, the header being row 1 (a spreadsheet's numbering)
    image_id: str
    file_path: str | None        # with a file key that holds a folder part: that file, absolute; else None
    values: dict                 # every column, stripped
    fields: dict                 # the design cells the row sets
    sha: str


class SampleSheet:
    """A read sample sheet (read_sample_sheet). ``row(image_id, path)`` gives one image's SampleRow."""

    def __init__(self, path, *, columns: list[str], key: str, rows: list[_Row], encoding: str,
                 warnings_: list[str]):
        self.path = Path(path)
        self.source = os.path.abspath(path)
        self.columns = list(columns)
        self.key = key
        self.encoding = encoding
        self.warnings = list(warnings_)
        self._rows = list(rows)
        self._by_id = {r.image_id: r for r in rows}

    def __repr__(self) -> str:
        return f"SampleSheet({self.source!r}, {len(self._rows)} rows, key {self.key})"

    def __len__(self) -> int:
        return len(self._rows)

    @property
    def rows(self) -> list[dict]:
        """Every row with every column, as read (stripped), in the file's order."""
        return [dict(r.values) for r in self._rows]

    @property
    def image_ids(self) -> list[str]:
        return [r.image_id for r in self._rows]

    @property
    def extra_columns(self) -> list[str]:
        """The columns no run reads (neither the key nor a design column), in the file's order: joined at collect."""
        return [c for c in self.columns if c != self.key and c not in DESIGN_COLUMNS]

    def extra_values(self) -> dict[str, dict[str, str]]:
        """``{image_id: {extra column: value}}`` for every row (an empty cell is "")."""
        cols = self.extra_columns
        return {r.image_id: {c: r.values.get(c, "") for c in cols} for r in self._rows}

    def row(self, image_id: str, path=None) -> SampleRow:
        """The SampleRow of the image ``image_id`` (at ``path``, when known): its row, or a row with no fields and
        row_sha256 None when the sheet has none for it. A row keyed by a file path (a ``file`` cell with a folder
        part) is the row of that file only, so with ``path`` given another file with the same id is not in the sheet."""
        r = self._by_id.get(str(image_id))
        if r is not None and r.file_path is not None and path is not None and not _same_file(r.file_path, path):
            r = None
        if r is None:
            return SampleRow({}, self.source, None)
        return SampleRow(dict(r.fields), self.source, r.sha)


def _same_file(a, b) -> bool:
    return os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))


def _has_folder(value: str) -> bool:
    return len(PureWindowsPath(value).parts) > 1 or len(PurePosixPath(value).parts) > 1


def _decode(data: bytes, where: str) -> tuple[str, str]:
    """(text, encoding): UTF-8 (with or without a byte-order mark), else cp1252."""
    if data[:2] in (b"\xff\xfe", b"\xfe\xff") or b"\x00" in data:
        raise SampleSheetError(f"{where}: the file looks like UTF-16 (Excel's \"Unicode Text\"), which a sample sheet "
                               f"is not: save it from Excel as \"CSV UTF-8 (Comma delimited)\"")
    try:
        return data.decode("utf-8-sig"), "utf-8-sig"
    except UnicodeDecodeError:
        if data.startswith(b"\xef\xbb\xbf"):          # a UTF-8 byte-order mark, then bytes that are not UTF-8
            raise SampleSheetError(f"{where}: the file starts with a UTF-8 byte-order mark but is not UTF-8 "
                                   f"throughout: save the sheet again as \"CSV UTF-8 (Comma delimited)\"") from None
    try:
        return data.decode("cp1252"), "cp1252"
    except UnicodeDecodeError as e:
        raise SampleSheetError(f"{where}: byte {e.start} (0x{data[e.start]:02x}) cannot be read as UTF-8 or as "
                               f"cp1252: save the sheet from Excel as \"CSV UTF-8 (Comma delimited)\"") from None


def _parse(text: str, where: str, base_dir: Path, encoding: str) -> SampleSheet:
    """The SampleSheet of a sheet's text; SampleSheetError naming ``where`` and the row or column for what cannot be
    read."""
    found: list[str] = []
    if encoding == "cp1252":
        found.append(f"{where}: not UTF-8, so it was read as cp1252 (Excel's plain \"CSV\" format on Windows); save it "
                     f"from Excel as \"CSV UTF-8 (Comma delimited)\" so every letter reads the same on every computer")
    lines = list(csv.reader(io.StringIO(text, newline="")))
    if not lines or not any(cell.strip() for cell in lines[0]):
        raise SampleSheetError(f"{where}: the sheet is empty, or its first line is: the first line must name the "
                               f"columns (image_id or file, then genotype, sex, ... and any columns of your own)")
    header = [h.strip() for h in lines[0]]
    if len(header) == 1 and (";" in header[0] or "\t" in header[0]):
        sep = "semicolons" if ";" in header[0] else "tabs"
        raise SampleSheetError(f"{where}: the columns are separated by {sep}, not commas (a semicolon list is what Excel "
                               f"saves in some languages): save it as \"CSV UTF-8 (Comma delimited)\"")
    body = [(n, cells) for n, cells in enumerate(lines[1:], 2) if any(c.strip() for c in cells)]
    keep: list[int] = []
    for i, name in enumerate(header):
        if name:
            keep.append(i)
            continue
        used = [n for n, cells in body if i < len(cells) and cells[i].strip()]
        if used:
            raise SampleSheetError(f"{where}: column {i + 1} has no name in the first line but holds a value in row "
                                   f"{used[0]}: name it, or delete the column")
    columns = [header[i] for i in keep]
    reserved = set(KEY_COLUMNS) | set(DESIGN_COLUMNS)
    for name in columns:
        if columns.count(name) > 1:
            raise SampleSheetError(f"{where}: the column {name} appears twice in the first line")
        if name != name.lower() and name.lower() in reserved:
            raise SampleSheetError(f"{where}: the column {name!r} must be written in lower case, {name.lower()!r}, to "
                                   f"be read as that column (the key and design columns are lower case)")
    key = "image_id" if "image_id" in columns else "file" if "file" in columns else None
    if key is None:
        raise SampleSheetError(f"{where}: no image_id or file column: the first line names {', '.join(columns)}; a "
                               f"sample sheet names each image by image_id (the file name without its suffix) or by "
                               f"file")
    for name in columns:
        if name in SHARED_META and name not in reserved:
            raise SampleSheetError(f"{where}: the column {name} is one every results table already has (the shared "
                                   f"metadata block of germquant.schema): rename it")
    rows: list[_Row] = []
    seen: dict[str, int] = {}
    for n, cells in body:
        extra = [c for i, c in enumerate(cells) if i >= len(header) and c.strip()]
        if extra:
            raise SampleSheetError(f"{where}: row {n} has more cells than the first line has columns ({extra[0]!r} is "
                                   f"beyond the last column)")
        values = {header[i]: (cells[i].strip() if i < len(cells) else "") for i in keep}
        ident = values[key]
        if not ident:
            raise SampleSheetError(f"{where}: row {n} has no {key}: every row names its image")
        if key == "file":
            # a Windows separator names the same file on Linux and macOS, where pathlib keeps a backslash in a name
            ident = ident.replace("\\", "/")
        image_id = image_stem(ident) if key == "file" else ident
        file_path = None
        if key == "file" and _has_folder(ident):
            p = Path(ident)
            file_path = str(p if p.is_absolute() else base_dir / p)
        if image_id in seen:
            raise SampleSheetError(f"{where}: rows {seen[image_id]} and {n} name the image {image_id} twice (an image "
                                   f"has one row)")
        seen[image_id] = n
        fields = {c: values[c] for c in DESIGN_COLUMNS if c in values and values[c]}
        if "acquisition_date" in fields and not iso_date(fields["acquisition_date"]):
            found.append(f"{where}: row {n} (image {image_id}): acquisition_date {fields['acquisition_date']!r} is not "
                         f"an ISO date (YYYY-MM-DD or YYYYMMDD); it is kept as given. Excel rewrites dates in the "
                         f"computer's own format when it saves: format the column as text, or write ISO dates")
        found += _germ_cell_warnings(where, n, image_id, fields)
        rows.append(_Row(n, image_id, file_path, values, fields, row_sha256(fields)))
    return SampleSheet(where, columns=columns, key=key, rows=rows, encoding=encoding, warnings_=found)


def read_sample_sheet(path) -> SampleSheet:
    """Read a sample sheet (module docstring). SampleSheetError for a sheet that cannot be read, OSError for a file
    that cannot be opened; each warning is in ``.warnings`` and is also issued as a SampleSheetWarning."""
    where = os.path.abspath(path)
    with open(long_path(path), "rb") as fh:
        data = fh.read()
    text, encoding = _decode(data, where)
    sheet = _parse(text, where, Path(where).parent, encoding)
    for w in sheet.warnings:
        warnings.warn(w, SampleSheetWarning, stacklevel=2)
    return sheet


def _cell_text(value, where: str) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        raise SampleSheetError(f"{where}: {value!r} is not text or a number")
    if isinstance(value, str):
        return value
    if isinstance(value, (int, float)) or type(value).__module__ == "numpy":
        try:
            f = float(value)
        except (TypeError, ValueError):
            raise SampleSheetError(f"{where}: {value!r} is not text or a number") from None
        if f != f:                                   # NaN: an empty cell of a table
            return ""
        return str(value)
    raise SampleSheetError(f"{where}: {value!r} is not text or a number")


def write_sample_sheet(path, rows: list[dict]) -> Path:
    """Write ``rows`` (one dict per image: the key, image_id or file, design columns and columns of your own) as a
    sample sheet: the columns in the order they first appear, ASCII, LF line ends, written atomically (a temporary
    file next to it, fsync, then os.replace, all through fsutil.long_path, retried for 2 s when Windows reports the
    file busy). SampleSheetError, before anything is written, for a value or column name outside ASCII (naming its
    cell), a value that is not text or a number, and a sheet read_sample_sheet would refuse."""
    path = Path(path)
    where = os.path.abspath(path)
    rows = list(rows or [])
    if not rows:
        raise SampleSheetError(f"{where}: no rows to write")
    columns: list[str] = []
    for r in rows:
        if not isinstance(r, dict):
            raise SampleSheetError(f"{where}: each row is a dict of column: value, not {type(r).__name__}")
        columns += [str(c) for c in r if str(c) not in columns]
    for c in columns:
        if not c.isascii():
            raise SampleSheetError(f"{where}: the column {c!r} is not ASCII; a sample sheet is written in ASCII")
    key = "image_id" if "image_id" in columns else "file" if "file" in columns else None
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(columns)
    for n, r in enumerate(rows, 2):
        ident = str(r.get(key, "") or "").strip() if key else ""
        if key and not ident:
            raise SampleSheetError(f"{where}: row {n} has no {key}: every row names its image")
        label = f"row {n} ({ident})" if ident else f"row {n}"
        cells = []
        for c in columns:
            text = _cell_text(r.get(c), f"{where}: {label}, column {c}")
            if not text.isascii():
                raise SampleSheetError(f"{where}: {label}, column {c}: {text!r} is not ASCII; write_sample_sheet writes "
                                       f"ASCII only (spell the letter out, or use an ASCII stand-in)")
            cells.append(text)
        writer.writerow(cells)
    text = buf.getvalue()
    sheet = _parse(text, where, Path(where).parent, "utf-8-sig")       # what read_sample_sheet would make of it
    del sheet
    return _atomic_write(path, text.encode("ascii"))


def _atomic_write(path: Path, data: bytes) -> Path:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp")
    try:
        with open(long_path(tmp), "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        for attempt in range(RETRIES + 1):
            try:
                _replace(long_path(tmp), long_path(path))
                return path
            except PermissionError:
                if attempt == RETRIES:
                    raise
                _sleep(RETRY_WAIT_S)
    except BaseException:
        try:
            os.remove(long_path(tmp))
        except OSError:
            pass
        raise
    return path                                   # not reached: the loop returns or raises
