#!/usr/bin/env python3
"""
Example: Classify structure and write annotated mmCIF output.

This example demonstrates how to:
1. Load a structure with CIF data
2. Classify dinucleotide steps
3. Generate DNATCO-extended mmCIF output with annotations

This is the Python equivalent of the C++ classify_and_write_cif example.
"""

import sys
import pyllka


def format_angle(radians, normalize=True):
    """
    Convert radians to degrees and format to 1 decimal place.

    Args:
        radians: Angle in radians
        normalize: If True, normalize to 0-360 range. If False, keep sign (for differences).
    """
    if normalize:
        # Use LLKA functions: normalize to 0-2π then convert to degrees
        degrees = pyllka.rad2deg(pyllka.full_angle_from_rad(radians))
    else:
        # Just convert to degrees, keeping the sign
        degrees = pyllka.rad2deg(radians)
    return f"{degrees:.1f}"


def format_float(value, decimals=3):
    """Format float to specified decimal places."""
    return f"{value:.{decimals}f}"


def generate_details(cs):
    """
    Generate the details field for a classified step based on violations.
    This follows the logic from the C example.

    Args:
        cs: ClassifiedStep object

    Returns:
        String with semicolon-separated violation codes, or "." if no violations
    """
    TORSION_NAMES = ["d1", "e1", "z1", "a2", "b2", "g2", "d2", "ch1", "ch2"]
    details = []

    # Check average nearest neighbor torsions violations
    if cs.violations & pyllka.CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT:
        for idx in range(9):
            if cs.violating_torsions_average & (1 << idx):
                details.append(f"cAn{TORSION_NAMES[idx]}")

    # Check nearest neighbor torsions violations
    if cs.violations & pyllka.CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT:
        for idx in range(9):
            if cs.violating_torsions_nearest & (1 << idx):
                details.append(f"cNn{TORSION_NAMES[idx]}")

    # Check distance violations
    if cs.violations & (pyllka.CLASSIFICATION_E_NN_TOO_LOW | pyllka.CLASSIFICATION_E_NN_TOO_HIGH):
        details.append("cNN")

    if cs.violations & (pyllka.CLASSIFICATION_E_CC_TOO_LOW | pyllka.CLASSIFICATION_E_CC_TOO_HIGH):
        details.append("cCC")

    if cs.violations & (pyllka.CLASSIFICATION_E_MU_TOO_LOW | pyllka.CLASSIFICATION_E_MU_TOO_HIGH):
        details.append("cmu")

    if cs.violations & pyllka.CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH:
        details.append("cMB")

    # Check pseudorotation violations
    if cs.violations & pyllka.CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT:
        details.append("cP")

    if cs.violations & pyllka.CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT:
        details.append("cP1")

    # Check delta torsion rejection
    if cs.violations & pyllka.CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED:
        details.append("DELTA")

    return ";".join(details) if details else "."


def format_cif_table(rows):
    """
    Format a table of CIF data with proper column alignment.
    Each row is a list of string values.
    Returns a list of formatted strings.
    """
    if not rows:
        return []

    # Calculate column widths
    num_cols = len(rows[0])
    col_widths = [0] * num_cols

    for row in rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(str(val)))

    # Format rows with proper spacing
    formatted_rows = []
    for row in rows:
        formatted_cols = []
        for i, val in enumerate(row):
            formatted_cols.append(str(val).ljust(col_widths[i]))
        formatted_rows.append(" ".join(formatted_cols))

    return formatted_rows


def get_step_alt_ids(step):
    """
    Find alt_ids for both nucleotides in a step by iterating through all atoms.
    This follows the logic from the C example's stepAltIds() function.

    Returns:
        tuple: (alt_id_first, alt_id_second) - each is either a character or empty string
    """
    if len(step) == 0:
        return "", ""

    # Get sequence ID of first residue
    seq_id_first = step.get_atom(0).label_seq_id
    alt_id_first = ""

    # Iterate through first nucleotide's atoms to find alt_id
    idx = 0
    while idx < len(step):
        atom = step.get_atom(idx)

        # Stop if we moved to the second nucleotide
        if atom.label_seq_id != seq_id_first:
            break

        # Keep updating alt_id until we find a non-null one
        if atom.label_alt_id and atom.label_alt_id != '\x00':
            alt_id_first = atom.label_alt_id
            break

        idx += 1

    # Scroll to the second nucleotide
    while idx < len(step) and step.get_atom(idx).label_seq_id == seq_id_first:
        idx += 1

    if idx >= len(step):
        # Only one nucleotide found - this shouldn't happen for a proper step
        return "", ""

    # Iterate through second nucleotide's atoms to find alt_id
    alt_id_second = ""
    while idx < len(step):
        atom = step.get_atom(idx)

        # Keep updating alt_id until we find a non-null one
        if atom.label_alt_id and atom.label_alt_id != '\x00':
            alt_id_second = atom.label_alt_id
            break

        idx += 1

    return alt_id_first, alt_id_second


def generate_dnatco_categories(classified_steps, steps, entry_id, ctx):
    """
    Generate DNATCO annotation categories as mmCIF text according to mmcif_ndb_ntc.dic.

    Args:
        classified_steps: List of ClassifiedStep objects
        steps: List of Structure objects (dinucleotide steps)
        entry_id: Entry ID from CIF file
        ctx: ClassificationContext used for classification

    Returns:
        String containing DNATCO categories in mmCIF format
    """

    # Count statistics
    total_steps = len(classified_steps)
    num_classified = sum(1 for cs in classified_steps if cs.violations == pyllka.CLASSIFICATION_OK)
    num_unclassified = total_steps - num_classified
    num_unclassified_rmsd_close = sum(1 for cs in classified_steps
                                      if cs.violations & pyllka.CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH)

    # Calculate average confal score and percentile
    avg_confal_result = pyllka.average_confal_attempted(classified_steps, ctx)
    avg_confal = avg_confal_result.score
    confal_percentile = int(avg_confal_result.percentile)

    output = []

    # Overall statistics category (according to dictionary lines 104-242)
    # This is NOT a loop - it's a single data block
    output.append(f"_ndb_struct_ntc_overall.entry_id  {entry_id}")
    output.append(f"_ndb_struct_ntc_overall.confal_score {format_float(avg_confal, 1)}")
    output.append(f"_ndb_struct_ntc_overall.confal_percentile {confal_percentile}")
    output.append(f"_ndb_struct_ntc_overall.ntc_version {pyllka.NTC_VERSION}")
    output.append(f"_ndb_struct_ntc_overall.cana_version {pyllka.CANA_VERSION}")
    output.append(f"_ndb_struct_ntc_overall.num_steps {total_steps}")
    output.append(f"_ndb_struct_ntc_overall.num_classified {num_classified}")
    output.append(f"_ndb_struct_ntc_overall.num_unclassified {num_unclassified}")
    output.append(f"_ndb_struct_ntc_overall.num_unclassified_rmsd_close {num_unclassified_rmsd_close}")
    output.append("")

    # NtC step category - residue identifiers (according to dictionary lines 244-520)
    output.append("loop_")
    output.append("_ndb_struct_ntc_step.id")
    output.append("_ndb_struct_ntc_step.name")
    output.append("_ndb_struct_ntc_step.PDB_model_number")
    output.append("_ndb_struct_ntc_step.label_entity_id_1")
    output.append("_ndb_struct_ntc_step.label_asym_id_1")
    output.append("_ndb_struct_ntc_step.label_seq_id_1")
    output.append("_ndb_struct_ntc_step.label_comp_id_1")
    output.append("_ndb_struct_ntc_step.label_alt_id_1")
    output.append("_ndb_struct_ntc_step.label_entity_id_2")
    output.append("_ndb_struct_ntc_step.label_asym_id_2")
    output.append("_ndb_struct_ntc_step.label_seq_id_2")
    output.append("_ndb_struct_ntc_step.label_comp_id_2")
    output.append("_ndb_struct_ntc_step.label_alt_id_2")
    output.append("_ndb_struct_ntc_step.auth_asym_id_1")
    output.append("_ndb_struct_ntc_step.auth_seq_id_1")
    output.append("_ndb_struct_ntc_step.auth_asym_id_2")
    output.append("_ndb_struct_ntc_step.auth_seq_id_2")
    output.append("_ndb_struct_ntc_step.PDB_ins_code_1")
    output.append("_ndb_struct_ntc_step.PDB_ins_code_2")

    # Collect step data rows for proper alignment
    step_rows = []
    for i, (cs, step) in enumerate(zip(classified_steps, steps)):
        step_id = i + 1

        # Get first and last atoms
        first_atom = step.get_atom(0)
        last_atom = step.get_atom(len(step) - 1)

        # Get alt_ids for both nucleotides by checking all atoms (not just first/last)
        alt_id_first, alt_id_second = get_step_alt_ids(step)

        # Format alt IDs for table output (use the detected alt_ids from all atoms)
        alt_1 = alt_id_first if alt_id_first else "."
        alt_2 = alt_id_second if alt_id_second else "."
        ins_1 = first_atom.pdbx_PDB_ins_code if first_atom.pdbx_PDB_ins_code else "."
        ins_2 = last_atom.pdbx_PDB_ins_code if last_atom.pdbx_PDB_ins_code else "."

        # Generate step name with alt_id as suffix with dot (e.g., DG.A) only when alt_id is present
        # Use auth_comp_id for residue names and auth_seq_id for residue numbers
        comp_1 = first_atom.auth_comp_id
        if alt_id_first:
            comp_1 = f"{comp_1}.{alt_id_first}"
        comp_2 = last_atom.auth_comp_id
        if alt_id_second:
            comp_2 = f"{comp_2}.{alt_id_second}"

        step_name = f"{entry_id.lower()}_{first_atom.label_asym_id}_{comp_1}_{first_atom.auth_seq_id}_{comp_2}_{last_atom.auth_seq_id}"

        step_rows.append([
            str(step_id),
            step_name,
            str(first_atom.pdbx_PDB_model_num),
            first_atom.label_entity_id,
            first_atom.label_asym_id,
            str(first_atom.label_seq_id),
            first_atom.auth_comp_id,
            alt_1,
            last_atom.label_entity_id,
            last_atom.label_asym_id,
            str(last_atom.label_seq_id),
            last_atom.auth_comp_id,
            alt_2,
            first_atom.auth_asym_id,
            str(first_atom.auth_seq_id),
            last_atom.auth_asym_id,
            str(last_atom.auth_seq_id),
            ins_1,
            ins_2
        ])

    # Format and add rows
    output.extend(format_cif_table(step_rows))

    # Step summary category (according to dictionary lines 522-775)
    output.append("loop_")
    output.append("_ndb_struct_ntc_step_summary.step_id")
    output.append("_ndb_struct_ntc_step_summary.assigned_CANA")
    output.append("_ndb_struct_ntc_step_summary.assigned_NtC")
    output.append("_ndb_struct_ntc_step_summary.confal_score")
    output.append("_ndb_struct_ntc_step_summary.euclidean_distance_NtC_ideal")
    output.append("_ndb_struct_ntc_step_summary.cartesian_rmsd_closest_NtC_representative")
    output.append("_ndb_struct_ntc_step_summary.closest_CANA")
    output.append("_ndb_struct_ntc_step_summary.closest_NtC")
    output.append("_ndb_struct_ntc_step_summary.closest_step_golden")

    summary_rows = []
    for i, cs in enumerate(classified_steps):
        step_id = i + 1
        golden = cs.closest_golden_step if cs.closest_golden_step else "?"
        confal = format_float(cs.confal_score.total, 0) if cs.violations == pyllka.CLASSIFICATION_OK else "?"

        summary_rows.append([
            str(step_id),
            pyllka.cana_to_name(cs.assigned_cana),
            pyllka.ntc_to_name(cs.assigned_ntc),
            confal,
            format_float(cs.euclidean_distance_ntc_ideal, 1),
            format_float(cs.rmsd_to_closest_ntc, 3),
            pyllka.cana_to_name(cs.closest_cana),
            pyllka.ntc_to_name(cs.closest_ntc),
            golden
        ])

    output.extend(format_cif_table(summary_rows))
    output.append("#")

    # Step parameters category (according to dictionary lines 777-1206)
    output.append("loop_")
    output.append("_ndb_struct_ntc_step_parameters.step_id")
    output.append("_ndb_struct_ntc_step_parameters.tor_delta_1")
    output.append("_ndb_struct_ntc_step_parameters.tor_epsilon_1")
    output.append("_ndb_struct_ntc_step_parameters.tor_zeta_1")
    output.append("_ndb_struct_ntc_step_parameters.tor_alpha_2")
    output.append("_ndb_struct_ntc_step_parameters.tor_beta_2")
    output.append("_ndb_struct_ntc_step_parameters.tor_gamma_2")
    output.append("_ndb_struct_ntc_step_parameters.tor_delta_2")
    output.append("_ndb_struct_ntc_step_parameters.tor_chi_1")
    output.append("_ndb_struct_ntc_step_parameters.tor_chi_2")
    output.append("_ndb_struct_ntc_step_parameters.dist_NN")
    output.append("_ndb_struct_ntc_step_parameters.dist_CC")
    output.append("_ndb_struct_ntc_step_parameters.tor_NCCN")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_delta_1")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_epsilon_1")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_zeta_1")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_alpha_2")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_beta_2")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_gamma_2")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_delta_2")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_chi_1")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_chi_2")
    output.append("_ndb_struct_ntc_step_parameters.diff_dist_NN")
    output.append("_ndb_struct_ntc_step_parameters.diff_dist_CC")
    output.append("_ndb_struct_ntc_step_parameters.diff_tor_NCCN")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_delta_1")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_epsilon_1")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_zeta_1")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_alpha_2")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_beta_2")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_gamma_2")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_delta_2")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_chi_1")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_chi_2")
    output.append("_ndb_struct_ntc_step_parameters.confal_dist_NN")
    output.append("_ndb_struct_ntc_step_parameters.confal_dist_CC")
    output.append("_ndb_struct_ntc_step_parameters.confal_tor_NCCN")
    output.append("_ndb_struct_ntc_step_parameters.details")

    params_rows = []
    for i, cs in enumerate(classified_steps):
        step_id = i + 1
        m = cs.metrics
        diff = cs.differences_from_ntc_averages
        conf = cs.confal_score

        # Generate details string for violations
        details = generate_details(cs)

        params_rows.append([
            str(step_id),
            format_angle(m.delta_1), format_angle(m.epsilon_1), format_angle(m.zeta_1),
            format_angle(m.alpha_2), format_angle(m.beta_2), format_angle(m.gamma_2),
            format_angle(m.delta_2), format_angle(m.chi_1), format_angle(m.chi_2),
            format_float(m.NN, 2), format_float(m.CC, 2), format_angle(m.mu),
            format_angle(diff.delta_1, normalize=False), format_angle(diff.epsilon_1, normalize=False), format_angle(diff.zeta_1, normalize=False),
            format_angle(diff.alpha_2, normalize=False), format_angle(diff.beta_2, normalize=False), format_angle(diff.gamma_2, normalize=False),
            format_angle(diff.delta_2, normalize=False), format_angle(diff.chi_1, normalize=False), format_angle(diff.chi_2, normalize=False),
            format_float(diff.NN, 2), format_float(diff.CC, 2), format_angle(diff.mu, normalize=False),
            format_float(conf.delta_1, 1), format_float(conf.epsilon_1, 1), format_float(conf.zeta_1, 1),
            format_float(conf.alpha_2, 1), format_float(conf.beta_2, 1), format_float(conf.gamma_2, 1),
            format_float(conf.delta_2, 1), format_float(conf.chi_1, 1), format_float(conf.chi_2, 1),
            format_float(conf.NN, 1), format_float(conf.CC, 1), format_float(conf.mu, 1),
            details
        ])

    output.extend(format_cif_table(params_rows))

    # Sugar parameters category (according to dictionary lines 1208-1531)
    output.append("loop_")
    output.append("_ndb_struct_sugar_step_parameters.step_id")
    output.append("_ndb_struct_sugar_step_parameters.P_1")
    output.append("_ndb_struct_sugar_step_parameters.tau_1")
    output.append("_ndb_struct_sugar_step_parameters.Pn_1")
    output.append("_ndb_struct_sugar_step_parameters.P_2")
    output.append("_ndb_struct_sugar_step_parameters.tau_2")
    output.append("_ndb_struct_sugar_step_parameters.Pn_2")
    output.append("_ndb_struct_sugar_step_parameters.nu_1_1")
    output.append("_ndb_struct_sugar_step_parameters.nu_1_2")
    output.append("_ndb_struct_sugar_step_parameters.nu_1_3")
    output.append("_ndb_struct_sugar_step_parameters.nu_1_4")
    output.append("_ndb_struct_sugar_step_parameters.nu_1_5")
    output.append("_ndb_struct_sugar_step_parameters.nu_2_1")
    output.append("_ndb_struct_sugar_step_parameters.nu_2_2")
    output.append("_ndb_struct_sugar_step_parameters.nu_2_3")
    output.append("_ndb_struct_sugar_step_parameters.nu_2_4")
    output.append("_ndb_struct_sugar_step_parameters.nu_2_5")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_1_1")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_1_2")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_1_3")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_1_4")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_1_5")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_2_1")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_2_2")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_2_3")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_2_4")
    output.append("_ndb_struct_sugar_step_parameters.diff_nu_2_5")

    sugar_rows = []
    for i, cs in enumerate(classified_steps):
        step_id = i + 1
        nu1 = cs.nu_angles_1
        nu2 = cs.nu_angles_2
        diff_nu1 = cs.nu_angle_differences_1
        diff_nu2 = cs.nu_angle_differences_2

        # Use VERY_TERSE format for sugar pucker names (C2end, C4exo instead of C2endo, C4exo)
        pucker_1 = pyllka.sugar_pucker_to_name(cs.sugar_pucker_1, pyllka.SugarPuckerNameBrevity.VERY_TERSE)
        pucker_2 = pyllka.sugar_pucker_to_name(cs.sugar_pucker_2, pyllka.SugarPuckerNameBrevity.VERY_TERSE)

        sugar_rows.append([
            str(step_id),
            format_angle(cs.ribose_pseudorotation_1), format_angle(cs.tau_1), pucker_1,
            format_angle(cs.ribose_pseudorotation_2), format_angle(cs.tau_2), pucker_2,
            format_angle(nu1.nu_0), format_angle(nu1.nu_1), format_angle(nu1.nu_2),
            format_angle(nu1.nu_3), format_angle(nu1.nu_4),
            format_angle(nu2.nu_0), format_angle(nu2.nu_1), format_angle(nu2.nu_2),
            format_angle(nu2.nu_3), format_angle(nu2.nu_4),
            format_angle(diff_nu1.nu_0, normalize=False), format_angle(diff_nu1.nu_1, normalize=False), format_angle(diff_nu1.nu_2, normalize=False),
            format_angle(diff_nu1.nu_3, normalize=False), format_angle(diff_nu1.nu_4, normalize=False),
            format_angle(diff_nu2.nu_0, normalize=False), format_angle(diff_nu2.nu_1, normalize=False), format_angle(diff_nu2.nu_2, normalize=False),
            format_angle(diff_nu2.nu_3, normalize=False), format_angle(diff_nu2.nu_4, normalize=False)
        ])

    output.extend(format_cif_table(sugar_rows))
    output.append("#")

    return "\n".join(output)


def main():
    if len(sys.argv) < 2:
        print("Usage: python classify_and_write_cif.py <structure.cif>", file=sys.stderr)
        print("\nThis tool classifies nucleic acid conformations and outputs", file=sys.stderr)
        print("mmCIF with DNATCO extended annotations.", file=sys.stderr)
        sys.exit(1)

    cif_path = sys.argv[1]

    # Load structure with CIF data
    print(f"Loading structure from {cif_path}...", file=sys.stderr)
    imported = pyllka.load_structure_with_cif(cif_path)
    structure = imported.get_structure()
    entry_id = imported.get_entry_id()

    if not entry_id:
        entry_id = "XXXX"

    print(f"Entry ID: {entry_id}", file=sys.stderr)
    print(f"Number of atoms: {len(structure)}", file=sys.stderr)

    # Split to dinucleotide steps
    print("Splitting structure to dinucleotide steps...", file=sys.stderr)
    steps = pyllka.split_to_dinucleotides(structure)
    print(f"Number of steps: {len(steps)}", file=sys.stderr)

    # Create classification context
    print("Loading classification data...", file=sys.stderr)
    ctx = pyllka.create_classification_context()

    # Classify all steps
    print("Classifying steps...", file=sys.stderr)
    classified_steps = pyllka.classify_structure(structure, ctx)

    # Count statistics
    assigned = sum(1 for cs in classified_steps if cs.violations == pyllka.CLASSIFICATION_OK)
    print(f"Assigned steps: {assigned}/{len(classified_steps)}", file=sys.stderr)

    # Get original CIF as string
    print("Generating annotated CIF...", file=sys.stderr)
    original_cif = pyllka.cif_to_string(imported, pretty=True)

    # Generate DNATCO categories
    dnatco_annotations = generate_dnatco_categories(classified_steps, steps, entry_id, ctx)

    # Combine original CIF with DNATCO annotations
    # Remove the trailing "#" from original CIF if present, then append annotations
    cif_lines = original_cif.rstrip().rstrip('#').rstrip()
    annotated_cif = cif_lines + "\n\n" + dnatco_annotations

    # Write to stdout
    print(annotated_cif)

    print(f"\nDone! Classified {len(classified_steps)} steps, {assigned} fully assigned.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
