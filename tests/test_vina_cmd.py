import re
from caverdock_lowerbound.vina import VinaBox, run_vina


class DummyProc:
    def __init__(self, stdout: str):
        self.stdout = stdout


def test_run_vina_builds_command(monkeypatch, tmp_path):
    calls = {}

    def fake_run(cmd, capture_output=False, text=False, check=False):
        # capture command arguments for assertions
        calls["cmd"] = cmd
        # minimal Vina-like output
        return DummyProc("""
-----+------------+----------+----------
    1       -7.8      0.000      0.000
""")

    monkeypatch.setattr("subprocess.run", fake_run)

    box = VinaBox(1.0, 2.0, 3.0, 10.0, 10.0, 10.0)
    score, out = run_vina(
        vina_bin="vina",
        receptor_pdbqt="r.pdbqt",
        ligand_pdbqt="l.pdbqt",
        box=box,
        out_pdbqt=str(tmp_path / "out.pdbqt"),
        log_path=str(tmp_path / "log.txt"),
        exhaustiveness=16,
        num_modes=1,
        seed=42,
        cpu=2,
    )

    assert score == -7.8
    cmd = calls["cmd"]
    assert cmd[0] == "vina"
    assert "--receptor" in cmd and "--ligand" in cmd
    assert "--center_x" in cmd and str(1.0) in cmd
    assert "--size_x" in cmd and str(10.0) in cmd
    assert "--seed" in cmd and "42" in cmd
    assert "--cpu" in cmd and "2" in cmd
