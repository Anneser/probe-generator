from src.maker import get_amplifier, filter_homopolymers, generate_probes, downsample_probes, format_output
import pytest


def test_get_amplifier_valid():
    spacer_up, spacer_down, init_up, init_down = get_amplifier("B1")
    assert spacer_up == "aa"
    assert init_up.startswith("GAGGAG")


def test_get_amplifier_invalid():
    with pytest.raises(ValueError):
        get_amplifier("BAD")


def test_filter_homopolymers_removes_long_homopolymers():
    seq = "ATCG" * 13 + "AAAAAAA" + "CGTA" * 13
    regions = filter_homopolymers(seq, polyAT=5, polyCG=5)
    assert all("AAAAAAA" not in seq[s:e] for s, e in regions)


def test_generate_probes_spacing():
    pos = [(0, 52), (53, 105), (108, 160)]
    result = generate_probes(pos, pause=0, cdna_len=200, min_spacing=2)
    assert result == [(0, 52), (108, 160)]

def test_generate_probes_all_included_with_zero_spacing():
    pos = [(0, 52), (53, 105), (108, 160)]
    result = generate_probes(pos, pause=0, cdna_len=200, min_spacing=0)
    assert result == [(0, 52), (53, 105), (108, 160)]

def test_downsample_to_2_of_5():
    coords = [(i, i+52) for i in range(0, 260, 52)]
    downsampled = downsample_probes(coords, 'y', 2)
    assert len(downsampled) == 2
    assert downsampled[0] in coords and downsampled[1] in coords


def test_downsample_disabled():
    coords = [(i, i+52) for i in range(0, 260, 52)]
    result = downsample_probes(coords, 'n', 2)
    assert result == coords


def test_format_output_basic():
    probe_data = {
        0: (0, "A"*25 + "nn" + "T"*25, 52)
    }
    seq = "T"*1000
    lines, anti, sense, _ = format_output(probe_data, seq, 1000, "B1", 100, "TestGene")
    assert len(lines) == 2
    assert "Pool name" not in lines[0]  # Just the value
    assert isinstance(anti, str)