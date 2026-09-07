"""
population_coverage.py
HLA Population Coverage Calculator.
Estimates the percentage of individuals in target human populations (World, East Asia, Europe)
covered by selected HLA Class I alleles carrying binder epitopes.

Methodology: Bui HH, Sidney J, Dinh K, Southwood S, Newman MJ, Sette A.
"Predicting population coverage of T-cell epitope-based diagnostics and vaccines."
BMC Bioinformatics. 2006;7:153. PMID: 16545123.
(Same Hardy-Weinberg cumulative-coverage approach used by the IEDB Population
Coverage Analysis Resource, https://tools.iedb.org/population/.)

Allele frequency data: Gonzalez-Galarza FF, et al. "Allele frequency net database
(AFND) 2020 update: gold-standard data classification, open access genotype data
and new query tools." Nucleic Acids Research. 2020;48(D1):D783-D788.
(allelefrequencies.net, licensed under CC BY 4.0.)

Frequencies below are weighted averages (by reported sample size n) of 2-field
HLA-A/B/C allele records from a snapshot of the AFND dataset assembled by the
community project https://github.com/slowkow/allelefrequencies (MIT license;
full text bundled at THIRD_PARTY_LICENSES/LICENSE-slowkow-allelefrequencies.txt).
"East Asia" / "Europe" are our own region groupings (population entries whose
AFND label starts with a country in that group); "World" averages across every
population in the dataset. These are approximations across the countries/studies
present in AFND, not a demographically-representative census of true world
population proportions — treat cumulative coverage % as an illustrative estimate,
not a precise epidemiological figure.
"""

from typing import List, Dict, Any, Optional

ALLELE_FREQUENCIES: Dict[str, Dict[str, float]] = {
    "World": {
        "HLA-A24:02": 0.1611, "HLA-A*24:02": 0.1611,
        "HLA-A02:01": 0.1519, "HLA-A*02:01": 0.1519,
        "HLA-A11:01": 0.1263, "HLA-A*11:01": 0.1263,
        "HLA-C07:02": 0.1207, "HLA-C*07:02": 0.1207,
        "HLA-C04:01": 0.1176, "HLA-C*04:01": 0.1176,
        "HLA-C07:01": 0.0906, "HLA-C*07:01": 0.0906,
        "HLA-C06:02": 0.0804, "HLA-C*06:02": 0.0804,
        "HLA-C01:02": 0.0757, "HLA-C*01:02": 0.0757,
        "HLA-C03:04": 0.0671, "HLA-C*03:04": 0.0671,
        "HLA-A01:01": 0.0641, "HLA-A*01:01": 0.0641,
        "HLA-B35:01": 0.0579, "HLA-B*35:01": 0.0579,
        "HLA-B51:01": 0.056, "HLA-B*51:01": 0.056,
        "HLA-B46:01": 0.0553, "HLA-B*46:01": 0.0553,
        "HLA-A03:01": 0.0545, "HLA-A*03:01": 0.0545,
        "HLA-A33:03": 0.0506, "HLA-A*33:03": 0.0506,
        "HLA-B40:02": 0.0502, "HLA-B*40:02": 0.0502,
        "HLA-B07:02": 0.0481, "HLA-B*07:02": 0.0481,
        "HLA-B40:01": 0.0443, "HLA-B*40:01": 0.0443,
    },
    "East Asia": {
        "HLA-A11:01": 0.2429, "HLA-A*11:01": 0.2429,
        "HLA-A24:02": 0.2267, "HLA-A*24:02": 0.2267,
        "HLA-C07:02": 0.17, "HLA-C*07:02": 0.17,
        "HLA-C01:02": 0.1604, "HLA-C*01:02": 0.1604,
        "HLA-B40:01": 0.1411, "HLA-B*40:01": 0.1411,
        "HLA-C03:04": 0.1187, "HLA-C*03:04": 0.1187,
        "HLA-C08:01": 0.1085, "HLA-C*08:01": 0.1085,
        "HLA-B46:01": 0.1064, "HLA-B*46:01": 0.1064,
        "HLA-A02:01": 0.0949, "HLA-A*02:01": 0.0949,
        "HLA-A02:07": 0.0808, "HLA-A*02:07": 0.0808,
        "HLA-B15:01": 0.0786, "HLA-B*15:01": 0.0786,
        "HLA-A33:03": 0.078, "HLA-A*33:03": 0.078,
        "HLA-B15:02": 0.0706, "HLA-B*15:02": 0.0706,
        "HLA-C03:03": 0.0697, "HLA-C*03:03": 0.0697,
        "HLA-B13:01": 0.0683, "HLA-B*13:01": 0.0683,
        "HLA-B58:01": 0.0662, "HLA-B*58:01": 0.0662,
        "HLA-A02:03": 0.055, "HLA-A*02:03": 0.055,
        "HLA-C12:02": 0.0548, "HLA-C*12:02": 0.0548,
    },
    "Europe": {
        "HLA-A02:01": 0.2399, "HLA-A*02:01": 0.2399,
        "HLA-C07:01": 0.1405, "HLA-C*07:01": 0.1405,
        "HLA-C04:01": 0.1296, "HLA-C*04:01": 0.1296,
        "HLA-A03:01": 0.1224, "HLA-A*03:01": 0.1224,
        "HLA-A01:01": 0.1211, "HLA-A*01:01": 0.1211,
        "HLA-A24:02": 0.1156, "HLA-A*24:02": 0.1156,
        "HLA-C07:02": 0.1017, "HLA-C*07:02": 0.1017,
        "HLA-C06:02": 0.0952, "HLA-C*06:02": 0.0952,
        "HLA-B07:02": 0.0872, "HLA-B*07:02": 0.0872,
        "HLA-B51:01": 0.0821, "HLA-B*51:01": 0.0821,
        "HLA-B08:01": 0.0746, "HLA-B*08:01": 0.0746,
        "HLA-C05:01": 0.0691, "HLA-C*05:01": 0.0691,
        "HLA-C02:09": 0.069, "HLA-C*02:09": 0.069,
        "HLA-B44:02": 0.0614, "HLA-B*44:02": 0.0614,
        "HLA-B35:01": 0.0585, "HLA-B*35:01": 0.0585,
        "HLA-B18:01": 0.0569, "HLA-B*18:01": 0.0569,
        "HLA-A11:01": 0.0549, "HLA-A*11:01": 0.0549,
        "HLA-A32:01": 0.0438, "HLA-A*32:01": 0.0438,
    },
}


def calculate_population_coverage(alleles: List[str], region: str = "World") -> Dict[str, Any]:
    """
    Calculate cumulative population coverage % for a given set of HLA Class I alleles.

    Formula: Coverage = 1 - prod_{i}(1 - f(A_i))^2 (assuming Hardy-Weinberg equilibrium)
    """
    freq_db = ALLELE_FREQUENCIES.get(region, ALLELE_FREQUENCIES["World"])
    default_freq = 0.08  # Default average frequency for alleles not in the table above

    unique_alleles = list(dict.fromkeys(alleles))
    uncovered_prob = 1.0
    allele_details = []

    for a in unique_alleles:
        # Match allele with or without asterisk
        freq = freq_db.get(a, default_freq)
        # Probability an individual does NOT carry this allele on either chromosome
        uncovered_prob *= (1.0 - freq) ** 2
        allele_details.append({
            "allele": a,
            "frequency": round(freq, 4),
            "individual_coverage_pct": round((1.0 - (1.0 - freq) ** 2) * 100, 2)
        })

    cumulative_coverage = round((1.0 - uncovered_prob) * 100, 2)

    return {
        "region": region,
        "cumulative_coverage_pct": cumulative_coverage,
        "n_alleles": len(unique_alleles),
        "allele_details": allele_details,
    }


if __name__ == "__main__":
    test_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-A*24:02"]
    cov_world = calculate_population_coverage(test_alleles, "World")
    cov_asia = calculate_population_coverage(test_alleles, "East Asia")
    print("=== HLA Population Coverage Test ===")
    print("World Coverage:", cov_world["cumulative_coverage_pct"], "%")
    print("East Asia Coverage:", cov_asia["cumulative_coverage_pct"], "%")
