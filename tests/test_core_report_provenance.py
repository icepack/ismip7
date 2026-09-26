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
    for k in ("ISMIP7_MELT_SLOPE", "ISMIP7_SIN_ALPHA_ANT", "ISMIP7_K_MELT",
              "ISMIP7_DELTAT_PER_BASIN_NPZ", "ISMIP7_K_PER_BASIN_NPZ"):
        monkeypatch.delenv(k, raising=False)
    env = core_report.effective_env()
    assert env["ISMIP7_MELT_SLOPE"] == "ant    # default (not exported)"
    assert env["ISMIP7_SIN_ALPHA_ANT"] == "0.005115    # default (not exported)"
    # the tracked calibration, by its path in the repository; the run's own
    # provenance line carries its sha256 and K
    assert env["ISMIP7_DELTAT_PER_BASIN_NPZ"] == (
        "antarctica/calibration/deltaT_per_basin_1000_K6.500e-05.npz"
        "    # default (not exported)")
    assert "ISMIP7_K_MELT" not in env
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", "local")
    assert core_report.effective_env()["ISMIP7_MELT_SLOPE"] == "local"
    monkeypatch.setenv("ISMIP7_K_PER_BASIN_NPZ", "/k/K_2500.npz")
    assert core_report.effective_env()["ISMIP7_DELTAT_PER_BASIN_NPZ"].startswith(
        "none, the legacy per-basin K")


def test_the_report_lifts_the_front_owner_an_external_law_names(tmp_path):
    from icepack2_tools.front import FRONT_OWNER_MARKER
    log = tmp_path / "run.log"
    line = f"{FRONT_OWNER_MARKER} level-set prescribed law (external: hfb sigma_max=0.15)"
    log.write_text(f"  {line}\nstep 1\n")
    assert core_report.front_owner(str(log)) == [line]


def test_the_report_resolves_the_step_left_at_its_default(monkeypatch):
    r"""Issue 20 moved the production step to 0.025; a report written from a
    shell that never exported it still states the step the drivers used."""
    monkeypatch.delenv("ISMIP7_DT", raising=False)
    assert core_report.effective_env()["ISMIP7_DT"] == "0.025    # default (not exported)"
    monkeypatch.setenv("ISMIP7_DT", "0.1")
    assert core_report.effective_env()["ISMIP7_DT"] == "0.1"


def test_the_report_is_named_for_the_run_resolution(monkeypatch):
    r"""The 32 km demonstration matrix keeps its names, and a production core
    is named and titled for the 1 km mesh it ran on."""
    monkeypatch.delenv("ISMIP7_LC", raising=False)
    assert core_report.resolution_km("results/hist_cesm2_waccm_32000_timeseries.csv") == "32"
    assert core_report.resolution_km("/r/ssp585_cesm2_waccm_1000_timeseries.csv") == "1"
    assert core_report.resolution_km("/r/ctrl_t1k_dthalf_2500_timeseries.csv") == "2.5"
    assert core_report.resolution_km("/r/renamed.csv") == "1"
    monkeypatch.setenv("ISMIP7_LC", "32000")
    assert core_report.resolution_km("/r/renamed.csv") == "32"
