from pathlib import Path
from caverdock_lowerbound.caver import parse_profile_csv, parse_centerline_pdb


def test_parse_profile_csv(tmp_path):
    p = tmp_path / "profile.csv"
    p.write_text("x,y,z,radius,distance\n0,0,0,2,0\n1,0,0,2,1\n")
    discs = parse_profile_csv(p)
    assert len(discs) == 2
    assert discs[1].x == 1.0
    assert discs[1].radius == 2.0


def test_parse_centerline_pdb(tmp_path):
    p = tmp_path / "centerline.pdb"
    # Two ATOM lines with B-factor as radius
    p.write_text(
        """
ATOM      1  C   XXX     1       0.000   0.000   0.000  1.00  2.00           C
ATOM      2  C   XXX     1       1.000   0.000   0.000  1.00  3.00           C
END
""".strip()
    )
    discs = parse_centerline_pdb(p)
    assert len(discs) == 2
    assert discs[0].radius == 2.0
    assert discs[1].radius == 3.0
    assert discs[1].distance > discs[0].distance
