r"""The committed core report states which forcing a run opened.

The submission README has to cite the forcing versions a run used (discussion
#37), and the freeze copy of the forcing is still re-synced by hand every
week or two, so the record has to be the run's own statement rather than the
date somebody last ran the audit. The run logs a marker line per forcing
variable and ``core_report.py`` lifts them, the way it lifts the climatology
pool. Results and logs are gitignored, so whatever the report does not carry
is lost.
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

import core_report  # noqa: E402
from icepack2_tools.climatology import CLIM_POOL_MARKER  # noqa: E402
from icepack2_tools.forcing import (FORCING_PROVENANCE_MARKER,  # noqa: E402
                                    describe_observational_forcing)


def test_the_report_lifts_each_forcing_line_once(tmp_path):
    log = tmp_path / "run.log"
    line = f"{FORCING_PROVENANCE_MARKER} ocean tf CESM2-WACCM ssp585 ocean v3"
    log.write_text(f"  {line}\nstep 1\n  {line}\n  {CLIM_POOL_MARKER} full 30/30 yr\n")
    assert core_report.forcing_provenance(str(log)) == [line]
    assert core_report.climatology_pool(str(log)) == [f"{CLIM_POOL_MARKER} full 30/30 yr"]


def test_a_missing_or_silent_log_still_writes_a_line(tmp_path):
    (absent,) = core_report.forcing_provenance(str(tmp_path / "nope.log"))
    assert "NOT RECORDED" in absent
    quiet = tmp_path / "old.log"
    quiet.write_text("step 1\n")
    (none,) = core_report.forcing_provenance(str(quiet))
    assert "predates the provenance banner" in none


def test_a_control_says_which_climatologies_it_runs_on(monkeypatch):
    monkeypatch.setenv("ISMIP7_OI_VERSION", "06_nov")
    smb, ocean = describe_observational_forcing(smb="RACMO2.4p1 SMB climatology 2000-2029", ocean=True)
    assert smb == f"{FORCING_PROVENANCE_MARKER} atmosphere RACMO2.4p1 SMB climatology 2000-2029"
    assert "release 06_nov" in ocean


def test_the_report_resolves_the_melt_knobs_left_at_their_defaults(monkeypatch):
    for k in ("ISMIP7_MELT_SLOPE", "ISMIP7_SIN_ALPHA_ANT", "ISMIP7_K_MELT"):
        monkeypatch.delenv(k, raising=False)
    env = core_report.effective_env()
    assert env["ISMIP7_MELT_SLOPE"] == "ant    # default (not exported)"
    assert env["ISMIP7_SIN_ALPHA_ANT"] == "0.005115    # default (not exported)"
    assert env["ISMIP7_K_MELT"] == "8.5e-05    # default (not exported)"
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", "local")
    assert core_report.effective_env()["ISMIP7_MELT_SLOPE"] == "local"


def test_the_report_lifts_the_front_owner_an_external_law_names(tmp_path):
    from icepack2_tools.front import FRONT_OWNER_MARKER
    log = tmp_path / "run.log"
    line = f"{FRONT_OWNER_MARKER} level-set prescribed law (external: hfb sigma_max=0.15)"
    log.write_text(f"  {line}\nstep 1\n")
    assert core_report.front_owner(str(log)) == [line]


def test_the_report_states_the_smb_feedback_left_at_its_default(monkeypatch):
    r"""The feedback is on by default and no runner exports the knob, so a
    report built from the exported environment alone would not say whether
    the run carried it."""
    monkeypatch.delenv("ISMIP7_SMB_ELEVATION_FEEDBACK", raising=False)
    env = core_report.effective_env()
    assert env["ISMIP7_SMB_ELEVATION_FEEDBACK"] == "1    # default (not exported)"
    monkeypatch.setenv("ISMIP7_SMB_ELEVATION_FEEDBACK", "0")
    assert core_report.effective_env()["ISMIP7_SMB_ELEVATION_FEEDBACK"] == "0"


def test_the_report_lifts_the_smb_feedback_banner_and_not_its_yearly_lines(tmp_path):
    from icepack2_tools.forcing import smb_feedback_banner
    log = tmp_path / "run.log"
    banner = smb_feedback_banner(None)
    log.write_text(f"  {banner}\n  dacabfdz feedback 2015: net +0.00 Gt/yr\n"
                   f"  {banner}\n")
    assert core_report.smb_feedback_record(str(log)) == [banner]
    old = tmp_path / "old.log"
    old.write_text("step 1\n")
    (none,) = core_report.smb_feedback_record(str(old))
    assert "ran without the feedback" in none
