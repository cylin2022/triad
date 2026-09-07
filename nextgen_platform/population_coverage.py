"""
population_coverage.py
HLA Population Coverage Calculator.
Estimates the percentage of individuals in target human populations (World, East Asia, Europe, North America, etc.)
covered by selected HLA Class I alleles carrying binder epitopes.
"""

from typing import List, Dict, Any, Optional

# Representative Phenotype Frequencies for Common HLA Class I Alleles (source: AFND / IEDB)
ALLELE_FREQUENCIES: Dict[str, Dict[str, float]] = {
    "World": {
        "HLA-A02:01": 0.28, "HLA-A*02:01": 0.28,
        "HLA-A24:02": 0.18, "HLA-A*24:02": 0.18,
        "HLA-A01:01": 0.14, "HLA-A*01:01": 0.14,
        "HLA-A03:01": 0.13, "HLA-A*03:01": 0.13,
        "HLA-A11:01": 0.12, "HLA-A*11:01": 0.12,
        "HLA-B07:02": 0.12, "HLA-B*07:02": 0.12,
        "HLA-B08:01": 0.09, "HLA-B*08:01": 0.09,
        "HLA-B15:01": 0.08, "HLA-B*15:01": 0.08,
        "HLA-B40:01": 0.07, "HLA-B*40:01": 0.07,
        "HLA-C07:01": 0.26, "HLA-C*07:01": 0.26,
        "HLA-C04:01": 0.16, "HLA-C*04:01": 0.16,
    },
    "East Asia": {
        "HLA-A24:02": 0.36, "HLA-A*24:02": 0.36,
        "HLA-A11:01": 0.28, "HLA-A*11:01": 0.28,
        "HLA-A02:01": 0.20, "HLA-A*02:01": 0.20,
        "HLA-A33:03": 0.14, "HLA-A*33:03": 0.14,
        "HLA-B40:01": 0.15, "HLA-B*40:01": 0.15,
        "HLA-B51:01": 0.11, "HLA-B*51:01": 0.11,
        "HLA-B15:02": 0.09, "HLA-B*15:02": 0.09,
        "HLA-C07:02": 0.24, "HLA-C*07:02": 0.24,
        "HLA-C01:02": 0.18, "HLA-C*01:02": 0.18,
    },
    "Europe": {
        "HLA-A02:01": 0.45, "HLA-A*02:01": 0.45,
        "HLA-A01:01": 0.27, "HLA-A*01:01": 0.27,
        "HLA-A03:01": 0.22, "HLA-A*03:01": 0.22,
        "HLA-B07:02": 0.23, "HLA-B*07:02": 0.23,
        "HLA-B08:01": 0.18, "HLA-B*08:01": 0.18,
        "HLA-B44:02": 0.15, "HLA-B*44:02": 0.15,
        "HLA-C07:01": 0.38, "HLA-C*07:01": 0.38,
    },
}


def calculate_population_coverage(alleles: List[str], region: str = "World") -> Dict[str, Any]:
    """
    Calculate cumulative population coverage % for a given set of HLA Class I alleles.

    Formula: Coverage = 1 - prod_{i}(1 - f(A_i))^2 (assuming Hardy-Weinberg equilibrium)
    """
    freq_db = ALLELE_FREQUENCIES.get(region, ALLELE_FREQUENCIES["World"])
    default_freq = 0.08  # Default average frequency for unlisted alleles

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
