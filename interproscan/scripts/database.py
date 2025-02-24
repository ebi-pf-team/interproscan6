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
        print(f"Error connecting to the sequence database:\n{e}")
        sys.exit(1)


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
                seq_id, seq_desc = line.lstrip(">").split(" ", maxsplit=1)
                current_sequence = {"id": seq_id, "description": seq_desc, "sequence": ""}
            else:
                current_sequence["sequence"] += line.strip()

    if current_sequence:
        add_sequence(conn, current_sequence, nucleic)
    conn.close()


def insert_translated_seqs(db_path, esl_output):
    conn = create_connection(db_path)
    current_sequence = None
    with open(esl_output, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if current_sequence:
                    add_sequence(conn, current_sequence, nucleic=False)
                    relate_orf_to_nucleotide(conn, current_sequence)
                header = re.match(ESL_TRANSLATE, line.strip())
                if not header:
                    print(f"Incorrectly formatted ESL_TRANSLATE output:\nLine: {line}")
                    sys.exit(1)
                orf_id, orf_desc = line.lstrip(">").split(" ", maxsplit=1)
                current_sequence = {
                    "id": orf_id,
                    "description": orf_desc, # the ID of the source
                    "sequence": "",
                    "source": header.group(1)  # id of the parent nt sequence
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
            print(f"Error inserting sequence {current_sequence['id']} into the database:\n{e}")
            sys.exit(1)
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
            print(f"Error inserting sequence meta data for {current_sequence['id']} into the database:\n{e}")
            sys.exit(1)
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
    nt_md5 = ""
    cur.execute("SELECT nt_md5 FROM NUCLEOTIDE WHERE id = ?", (current_sequence["source"],))
    cur.fetchone()

    # add the nt seq md5 to the Protein table
    try:
        cur.execute("UPDATE PROTEIN_TO_NUCLEOTIDE SET nt_md5 = ? WHERE protein_md5 = ?", (protein_md5, nt_md5))
    except sqlite3.Error as e:
        if str(e).strip().startswith("UNIQUE constraint failed"):
            pass
        else:
            print(f"Error updating PROTEIN_TO_NUCLEOTIDE table with protein {protein_md5} and nucleotide {nt_md5}:\n{e}")
            sys.exit(1)
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
                nt_md5       VARCHAR PRIMARY KEY,
                id           VARCHAR,
                description  VARCHAR,
                FOREIGN KEY (nt_md5) REFERENCES NUCLEOTIDE_SEQUENCE(nt_md5)
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
        print(f"Error creating tables in the sequence database:\n{e}")
        sys.exit(1)


def build_batches(db_path, batch_size):
    """Build batches of FASTA file for the IPS6 sequence analysis."""
    def write_fasta(batch, batch_index, line_length=60):
        fasta = f"sequences.{batch_index}.fasta"
        with open(fasta, "w") as fh:
            for row in batch:
                md5 = row[0]
                sequence = '\n'.join([row[1][i:i+line_length] for i in range(0, len(row[1]), line_length)])
                fh.write(f">{md5}\n{sequence}\n")

    conn = create_connection(db_path)
    cur = conn.cursor()
    offset, batch_index = 0, 1
    while True:
        cur.execute("SELECT protein_md5, sequence FROM PROTEIN_SEQUENCE LIMIT ? OFFSET ?", (batch_size, offset))
        batch = cur.fetchall()

        if not batch:
            break

        write_fasta(batch, batch_index)
        offset += batch_size
        batch_index += 1

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
            print("The --fasta flag must be used with the populate_sequences method to provide the FASTA file of sequences")
            sys.exit(1)
    elif args.method == "update_orfs":
        if args.fasta:
            insert_translated_seqs(args.db_path, args.fasta)
        else:
            print("The --fasta flag must be used with the update_orfs method to provide the FASTA file of sequences")
            sys.exit(1)
    elif args.method == "build_batches":
        if args.batch_size:
            build_batches(args.db_path, args.batch_size)
        else:
            print("The --batch_size flag must be used with the build_batches method")
            sys.exit(1)
