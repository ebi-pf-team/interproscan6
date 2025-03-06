import argparse
import hashlib
import re
import sqlite3
import sys


ESL_TRANSLATE = re.compile(r"^>orf\d+\s+source=(.+?)\s+coords=\d+\.\.\d+\s+length=\d+\s+frame=\d\s+desc=.+$")


def create_connection(db_path):
    """Connect to a local IPS6 sequence database."""
    try:
        conn = sqlite3.connect(db_path)
        return conn
    except sqlite3.Error as e:
        raise sqlite3.DatabaseError(f"Error connecting to the sequence database:\n{e}")


def get_md5(sequence):
    return hashlib.md5(sequence.encode()).hexdigest()


def populate_sequences(db_path, fasta, nucleic):
    """Parse input FASTA file and insert into the sequence database."""
    conn = create_connection(db_path)

    build_tables(conn)

    current_sequence = None
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if current_sequence:
                    add_sequence(conn, current_sequence, nucleic)
                seq_id, seq_desc = line.lstrip(">").rstrip("\n").split(" ", maxsplit=1)
                current_sequence = {"id": seq_id, "description": seq_desc, "sequence": ""}
            else:
                current_sequence["sequence"] += line.strip()

    if current_sequence:
        add_sequence(conn, current_sequence, nucleic)
    conn.close()


def insert_orfs(db_path, esl_output):
    conn = create_connection(db_path)
    current_sequence = None
    with open(esl_output, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if current_sequence:
                    add_sequence(conn, current_sequence, nucleic=False)
                    relate_orf_to_nucleotide(conn, current_sequence)
                header = re.match(ESL_TRANSLATE, line.strip())
                if not header: # mv to raise
                    raise ValueError(f"Incorrectly formatted ESL_TRANSLATE output:\nLine: {line}")
                orf_id, orf_desc = line.lstrip(">").rstrip("\n").split(" ", maxsplit=1)
                current_sequence = {
                    "id": orf_id,
                    "description": orf_desc, # the ID of the source
                    "sequence": "",
                    "source": header.group(1).lstrip("source=")  # id of the parent nt sequence
                }
            else:
                current_sequence["sequence"] += line.strip()

    if current_sequence:
        add_sequence(conn, current_sequence, nucleic=False)
        relate_orf_to_nucleotide(conn, current_sequence)

    conn.close()


def add_sequence(conn, current_sequence, nucleic):
    md5 = get_md5(current_sequence["sequence"])
    cur = conn.cursor()
    try:
        if nucleic:
            cur.execute("INSERT INTO NUCLEOTIDE_SEQUENCE (nt_md5, sequence) VALUES (?, ?)",
                        (md5, current_sequence['sequence']))
        else:
            cur.execute("INSERT INTO PROTEIN_SEQUENCE (protein_md5, sequence) VALUES (?, ?)",
                        (md5, current_sequence['sequence']))
    except sqlite3.Error as e:
        if str(e).strip().startswith("UNIQUE constraint failed"):
            pass
        else:
            raise sqlite3.DatabaseError(f"Error inserting sequence {current_sequence['id']} into the database:\n{e}")
    try:
        if nucleic:
            cur.execute("INSERT INTO NUCLEOTIDE (nt_md5, id, description) VALUES (?, ?, ?)",
                        (md5, current_sequence['id'], current_sequence['description']))
        else:
            cur.execute("INSERT INTO PROTEIN (protein_md5, id, description) VALUES (?, ?, ?)",
                        (md5, current_sequence['id'], current_sequence['description']))
    except sqlite3.Error as e:
        if str(e).strip().startswith("UNIQUE constraint failed"):
            pass
        else:
            raise sqlite3.DatabaseError(f"Error inserting sequence meta data for {current_sequence['id']} into the database:\n{e}")
    conn.commit()
    cur.close()


def relate_orf_to_nucleotide(conn, current_sequence):
    """Add the md5 of the translated ORF protein seq to the nt table.
    source = the parent nt seq id.
    We store the nt seq md5 and protein md5 in a relationship table because multiple
    protein seqs can originate from one nt seq, and multiple nt seqs can produce
    the same protein seq, and these relationships may overlap"""
    protein_md5 = get_md5(current_sequence["sequence"])
    cur = conn.cursor()

    # get the nt seq md5
    try:
        nt_md5 = ""
        cur.execute("SELECT nt_md5 FROM NUCLEOTIDE WHERE id = ?", (current_sequence["source"],))
        nt_md5 = cur.fetchone()[0]
    except TypeError as e:
        raise ValueError(f"Could not get md5 for nt id: {current_sequence['source']}")
    # add the nt seq md5 to the Protein table
    try:
        cur.execute("INSERT INTO PROTEIN_TO_NUCLEOTIDE (protein_md5, nt_md5) VALUES (?, ?)", (protein_md5, nt_md5))
    except sqlite3.Error as e:
        if str(e).strip().startswith("UNIQUE constraint failed"):
            pass
        else:
            raise Exception(f"Error updating PROTEIN_TO_NUCLEOTIDE table with protein {protein_md5} and nucleotide {nt_md5}:\n{e}")
    conn.commit()
    cur.close()


def build_tables(conn):
    cur = conn.cursor()
    try:
        cur.executescript("""
            CREATE TABLE IF NOT EXISTS PROTEIN_SEQUENCE (
                protein_md5  VARCHAR PRIMARY KEY,
                sequence     TEXT UNIQUE
            );
            CREATE TABLE IF NOT EXISTS PROTEIN (
                protein_md5  VARCHAR,
                id           VARCHAR,
                description  VARCHAR,
                FOREIGN KEY (protein_md5) REFERENCES PROTEIN_SEQUENCE(protein_md5),
                UNIQUE (protein_md5, id, description)
            );
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE_SEQUENCE (
                nt_md5       VARCHAR PRIMARY KEY,
                sequence     TEXT
            );
            CREATE TABLE IF NOT EXISTS NUCLEOTIDE (
                nt_md5       VARCHAR,
                id           VARCHAR,
                description  VARCHAR,
                FOREIGN KEY (nt_md5) REFERENCES NUCLEOTIDE_SEQUENCE(nt_md5),
                UNIQUE (nt_md5, id, description)
            );
            CREATE TABLE IF NOT EXISTS PROTEIN_TO_NUCLEOTIDE (
                protein_md5 VARCHAR,
                nt_md5 VARCHAR,
                FOREIGN KEY (protein_md5) REFERENCES PROTEIN(protein_md5),
                FOREIGN KEY (nt_md5) REFERENCES NUCLEOTIDE(nt_md5),
                UNIQUE (protein_md5, nt_md5)
            );
        """)
    except sqlite3.Error as e:
        raise sqlite3.DatabaseError(f"Error creating tables in the sequence database:\n{e}")


def build_batches(db_path, batch_size, nucleic):
    """Build batches of FASTA files for IPS6 sequence analysis, ensuring that all children protein seqs
    associated with a parent nt seq are stored within the same batch."""

    def write_fasta(batch, batch_index, nucleic, line_length=60):
        """Write batch to a FASTA file."""
        fasta_file = f"sequences.{batch_index}.fasta"
        with open(fasta_file, "w") as fh:
            for row in batch:
                md5_index = 0 if nucleic else 1
                md5 = row[md5_index]
                seq = '\n'.join([row[-1][i:i+line_length] for i in range(0, len(row[-1]), line_length)])
                fh.write(f">{md5}\n{seq}\n")

    conn = create_connection(db_path)
    cur = conn.cursor()

    batch = []
    current_md5 = None
    batch_index = 1

    query = """
        SELECT P2N.nt_md5, S.protein_md5, S.sequence
        FROM PROTEIN_SEQUENCE AS S
        LEFT JOIN PROTEIN AS P ON S.protein_md5 = P.protein_md5
        LEFT JOIN PROTEIN_TO_NUCLEOTIDE AS P2N ON P.protein_md5 = P2N.protein_md5
        WHERE P2N.nt_md5 IS NOT NULL
        ORDER BY P2N.nt_md5
    """ if nucleic else "SELECT NULL AS nt_md5, protein_md5, sequence FROM PROTEIN_SEQUENCE ORDER BY protein_md5"

    cur.execute(query)

    for nt_md5, protein_md5, sequence in cur:
        # if we encounter a new nt_md5 and the batch is full --> write to a new batch
        if current_md5 and current_md5 != nt_md5 and len(batch) >= batch_size:
            write_fasta(batch, batch_index, nucleic)
            batch.clear()
            batch_index += 1

        batch.append((nt_md5, protein_md5, sequence))
        current_md5 = nt_md5 if nucleic else protein_md5

    # Write remaining batch
    if batch:
        write_fasta(batch, batch_index, nucleic)

    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='ips6.database.py',
        description='Handle interactions with the ISP6 internal seq database',
    )
    parser.add_argument(
        "db_path",
        type=str,
        help="Path to local sqlite db"
    )
    parser.add_argument(
        "method",
        choices=["populate_sequences", "update_orfs", "build_batches"],
        help=(
            "Select one method.\n"
            "populate_sequences: Insert sequences from a --fasta file to the db\n"
            "update_orfs: Insert translated ORFs from the --fasta esl.translate output to the db\n"
            "build_batches: Build sequence analysis FASTA batch files of --batch_size sequences per file"
        )
    )
    parser.add_argument(
        "--fasta",
        type=str,
        default=None,
        help="Str repr of the path to the FASTA files of seqs to be inserted into the db"
    )
    parser.add_argument(
        "--nucleic",
        action="store_true",
        help="Inserting nucleotide sequences into the db. Used with --populate_sequences"
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=None,
        help="Number of seqs per FASTA batch file"
    )
    args = parser.parse_args()

    if args.method == "populate_sequences":
        if args.fasta:
            populate_sequences(args.db_path, args.fasta, args.nucleic)
        else:
            raise ValueError("The --fasta flag must be used with the populate_sequences method to provide the FASTA file of sequences")
    elif args.method == "update_orfs":
        if args.fasta:
            for fasta in args.fasta.split(","):
                insert_orfs(args.db_path, fasta.lstrip("[").rstrip("]").strip())
        else:
            raise ValueError("The --fasta flag must be used with the update_orfs method to provide the FASTA file of sequences")
    elif args.method == "build_batches":
        if args.batch_size:
            build_batches(args.db_path, args.batch_size, args.nucleic)
        else:
            raise ValueError("The --batch_size flag must be used with the build_batches method")
