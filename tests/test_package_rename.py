"""The public namespace can change without invalidating archived chemistry."""
import json
from importlib.resources import files
from pathlib import Path

import graft
from graft.artifacts import (
    aam_from_record, aam_record, read_aam_checkpoint, read_graph_checkpoint,
    read_raw_cut, write_aam_checkpoint,
)
from graft.growth import native
from graft.postprocessing import decode_events

FIXTURES = Path(__file__).parent / "fixtures" / "package_rename"


def test_old_checkpoint_globals_resolve_to_graft_and_preserve_events(tmp_path):
    archived = json.loads((FIXTURES / "legacy-aam.json").read_text())
    result = read_aam_checkpoint(FIXTURES / "legacy-aam.pkl.gz")
    assert isinstance(result, graft.AAMResult)
    assert isinstance(result.graph, graft.AAMSearchGraph)
    assert type(result.graph).__module__.startswith("graft.")
    # JSON-normalize tuple/list differences, without changing the payload.
    assert json.loads(json.dumps(aam_record(result))) == archived
    assert json.loads(json.dumps(aam_record(aam_from_record(archived)))) == archived
    expected = json.loads((FIXTURES / "legacy-events.json").read_text())
    assert sorted(candidate.id for candidate in decode_events(result).candidates) == expected
    current = tmp_path / "graft.pkl.gz"
    write_aam_checkpoint(result, current)
    assert read_aam_checkpoint(current).graph == result.graph


def test_both_old_binary_cut_formats_keep_their_graph():
    expected = read_aam_checkpoint(FIXTURES / "legacy-aam.pkl.gz").graph
    assert read_graph_checkpoint(FIXTURES / "legacy-finalized.pkl.gz") == expected
    assert read_raw_cut(FIXTURES / "legacy.raw.pkl.gz") == expected


def test_installed_viewer_assets_are_available_under_graft():
    assert files("graft").joinpath("static/3Dmol-min.js").is_file()
    assert files("graft").joinpath("static/aam_search.html").is_file()


def test_native_setting_accepts_legacy_fallback_and_prefers_new_name(monkeypatch):
    monkeypatch.setattr(native, "_engine", object())
    monkeypatch.delenv("GRAFT_NATIVE", raising=False)
    monkeypatch.delenv("RXN_CORE_NATIVE", raising=False)
    assert native.available()
    monkeypatch.setenv("RXN_CORE_NATIVE", "0")
    assert not native.available()
    monkeypatch.setenv("GRAFT_NATIVE", "1")
    assert native.available()
    monkeypatch.setenv("GRAFT_NATIVE", "0")
    monkeypatch.setenv("RXN_CORE_NATIVE", "1")
    assert not native.available()
