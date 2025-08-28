from caverdock_lowerbound.pipeline import run_pipeline, PipelineConfig
from caverdock_lowerbound.charge_providers import ChargeProvider, DiscContext


class DummyProvider(ChargeProvider):
    def get_receptor_pdbqt(self, context: DiscContext) -> str:
        return "rec.pdbqt"

    def get_ligand_pdbqt(self, context: DiscContext) -> str:
        return "lig.pdbqt"


class DummyProc:
    def __init__(self, stdout: str):
        self.stdout = stdout


def test_pipeline_with_mocks(monkeypatch, tmp_path):
    # Mock discs loading
    from caverdock_lowerbound import caver as caver_mod

    def fake_load_discs(_dir, tunnel_index=1):
        return [caver_mod.Disc(1, 0, 0, 0, 2.0, 0.0), caver_mod.Disc(2, 1, 0, 0, 2.0, 1.0)]

    monkeypatch.setattr(caver_mod, "load_discs", fake_load_discs)

    # Mock Vina to avoid external call
    import caverdock_lowerbound.vina as vina_mod

    def fake_run_vina(**kwargs):
        return (-6.5, "mocked")

    monkeypatch.setattr(vina_mod, "run_vina", fake_run_vina)

    # Run pipeline
    res = run_pipeline(
        receptor_input="rec.pdb",
        ligand_input="lig.pdb",
        caver_output_dir=str(tmp_path),
        out_dir=str(tmp_path / "out"),
        tunnel_index=1,
        charge_provider=DummyProvider(),
        config=PipelineConfig(vina_bin="vina", exhaustiveness=4, cpu=1),
    )

    assert "results_df" in res
    df = res["results_df"]
    assert len(df) == 2
    assert df["vina_score"].iloc[0] == -6.5
