import importlib.util
import json
import os

import pandas as pd
import streamlit as st

DB_PATH = os.path.join(os.path.dirname(__file__), "onquadro-aligner.json")

spec = importlib.util.spec_from_file_location("find", os.path.join(os.path.dirname(__file__), "02-find.py"))
find = importlib.util.module_from_spec(spec)
spec.loader.exec_module(find)


@st.cache_data
def load_database(path):
    with open(path) as f:
        return json.load(f)


def detect_molecule_filter(raw_sequence):
    has_lower = any(c.islower() for c in raw_sequence if c.isalpha())
    has_upper = any(c.isupper() for c in raw_sequence if c.isalpha())
    if has_lower and not has_upper:
        return "DNA"
    if has_upper and not has_lower:
        return "RNA"
    return None


def build_g4composer_content(row):
    g4c = row.get("g4c")
    if not g4c:
        return None
    base = row["template"]
    sequence = row["Input sequence"]
    if row["Molecule"] == "RNA":
        seq = sequence.upper()
    else:
        seq = sequence.lower()
    structure = ["."] * len(sequence)
    for i in row["_combination"]:
        structure[i] = "^"
    structure = "".join(structure)
    score_line = (
        f"# tract_distance={row['Tract distance']} "
        f"linker_distance={row['Linker distance']} "
        f"viability={row['Viability']}"
    )
    if row["Loop lengths"]:
        score_line += f" loop_lengths={row['Loop lengths']}"
    if row["Topology"]:
        score_line += f" topology={row['Topology']}"
    score_line += f" template={base}"
    return find.serialize_g4composer_entry(
        base,
        seq,
        structure,
        "." * len(sequence),
        "." * len(sequence),
        g4c["orient"],
        g4c["rise"],
        g4c["twist"],
        g4c["path"],
        score_line,
    )


VIABILITY_COLORS = {
    "viable": "🟢",
    "marginal": "🟡",
    "unknown": "⚪",
    "not_viable": "🔴",
    "n/a": "⚪",
}


def main():
    st.set_page_config(page_title="onquadro-aligner", layout="wide")
    st.title("onquadro-aligner")
    st.caption("Search a sequence against known single-chain G-quadruplex structures")

    with st.form("search_form"):
        raw_sequence = st.text_input(
            "Sequence",
            placeholder="GGGTACCCGGGTGAGGTGCGGGGT",
            help="Uppercase = RNA, lowercase = DNA. Mixed case searches both.",
        )
        col1, col2 = st.columns(2)
        with col1:
            tetrad_count = st.number_input(
                "Tetrad count",
                min_value=0,
                value=0,
                help="0 = maximum possible (default)",
            )
        with col2:
            list_all = st.checkbox(
                "List all matching files",
                help="Show all matching file names instead of one per hit",
            )
        submitted = st.form_submit_button("Search")

    if submitted:
        if not raw_sequence.strip():
            st.warning("Please enter a sequence.")
            st.session_state.pop("search_result", None)
            return

        sequence = raw_sequence.upper()
        g_indices = [i for i, c in enumerate(sequence) if c == "G"]
        max_tetrad_count = len(g_indices) // 4

        if max_tetrad_count == 0:
            st.warning("Sequence contains fewer than 4 G's — cannot form a G-quadruplex.")
            st.session_state.pop("search_result", None)
            return

        tc_list = [int(tetrad_count)] if tetrad_count > 0 else [max_tetrad_count]

        data = load_database(DB_PATH)
        molecule_filter = detect_molecule_filter(raw_sequence)
        if molecule_filter:
            data = [d for d in data if d["molecule"] == molecule_filter]

        with st.spinner("Searching..."):
            result = find.match_quadruplexes(
                data,
                sequence,
                g_indices,
                tc_list,
                list_all,
                parallel=False,
            )

        if not result:
            st.warning("No matches found.")
            st.session_state.pop("search_result", None)
            return

        result.sort(
            key=lambda row: (
                0 if row["Tract distance"] == 0 and row["Linker distance"] == 0 else 1,
                find.VIABILITY_RANK.get(row["Viability"], 2),
                row["Tract distance"],
                row["Linker distance"],
            )
        )
        st.session_state.search_result = result

    result = st.session_state.get("search_result")
    if not result:
        st.info("Enter a sequence above and click Search.")
        return

    display_rows = []
    for i, row in enumerate(result):
        display_rows.append(
            {
                "": i,
                "Perfect": "✅" if row["Tract distance"] == 0 and row["Linker distance"] == 0 else "",
                "Viability": f"{VIABILITY_COLORS.get(row['Viability'], '⚪')} {row['Viability']}",
                "Tract": row["Tract distance"],
                "Linker": row["Linker distance"],
                "Tetrads": row["Tetrad count"],
                "Template": row["template"],
                "Loops": row["Loop lengths"],
                "Topology": row["Topology"],
                "QRS": f"`{row['QRS']}`",
            }
        )
    display_df = pd.DataFrame(display_rows)

    st.subheader(f"Results ({len(result)} matches)")
    event = st.dataframe(
        display_df,
        on_select="rerun",
        selection_mode="single-row",
        width="stretch",
        hide_index=True,
        column_config={
            "": st.column_config.NumberColumn(width="small"),
            "Perfect": st.column_config.TextColumn(width="small"),
            "Viability": st.column_config.TextColumn(width="medium"),
            "Tract": st.column_config.NumberColumn(width="small"),
            "Linker": st.column_config.NumberColumn(width="small"),
            "Tetrads": st.column_config.NumberColumn(width="small"),
            "Template": st.column_config.TextColumn(width="medium"),
            "Loops": st.column_config.TextColumn(width="small"),
            "Topology": st.column_config.TextColumn(width="small"),
            "QRS": st.column_config.TextColumn(width="large"),
        },
    )

    selected_rows = event.selection.rows
    if selected_rows:
        selected_idx = selected_rows[0]
        row = result[selected_idx]
        st.divider()
        st.subheader(f"g4composer output — {row['template']}")
        content = build_g4composer_content(row)
        if content:
            st.code(content, language="text")
            st.download_button(
                "Download .inp",
                data=content,
                file_name=f"{row['template']}.inp",
                mime="text/plain",
            )
        else:
            st.info("No g4composer data available for this match.")


if __name__ == "__main__":
    main()
