import hashlib
import re
import sqlite3
import sys


ESL_TRANSLATE = re.compile(r"^>(orf\d+)\s+source=(.+?)\s+coords=\d+\.\.\d+\s+length=\d+\s+frame=\d\s+desc=.+$")


def create_connection(db_path):
    """Connect to a local IPS6 sequence database."""
    try:
        conn = sqlite3.connect(db_path)
        return conn
    except sqlite3.Error as e:
        print(f"Error connecting to the sequence database:\n{e}")
        sys.exit(1)


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
                current_sequence = {
                    "id": header.group(1),
                    "description": header.group(2).rstrip("\n"),
                    "sequence": "",
                }
            else:
                current_sequence["sequence"] += line.strip()

    if current_sequence:
        add_sequence(conn, current_sequence, nucleic=False)
        relate_orf_to_nucleotide(conn, current_sequence)

    conn.close()


def add_sequence(conn, current_sequence, nucleic):
    def get_md5(sequence):
        return hashlib.md5(sequence.encode()).hexdigest()

    md5 = get_md5(current_sequence["sequence"])
    cur = conn.cursor()
    try:
        if nucleic:
            cur.execute("INSERT INTO NUCLEOTIDE_SEQUENCE (nt_md5, sequence) VALUES (?, ?)",
                        (md5, current_sequence['sequence']))
            cur.execute("INSERT INTO NUCLEOTIDE (protein_md5, nt_md5, id, description) VALUES (?, ?, ?, ?)",
                        (None, md5, current_sequence['id'], current_sequence['description']))
        else:
            cur.execute("INSERT INTO PROTEIN_SEQUENCE (protein_md5, sequence) VALUES (?, ?)",
                        (md5, current_sequence['sequence']))
            cur.execute("INSERT INTO PROTEIN (protein_md5, id, description) VALUES (?, ?, ?)",
                        (md5, current_sequence['id'], current_sequence['description']))
    except sqlite3.Error as e:
        if str(e).strip().startswith("UNIQUE constraint failed"):
            pass
        else:
            print(f"Error inserting sequence {current_sequence['id']} into the database:\n{e}")
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
                protein_md5  VARCHAR,
                nt_md5       VARCHAR PRIMARY KEY,
                id           VARCHAR,
                description  VARCHAR,
                FOREIGN KEY (protein_md5) REFERENCES PROTEIN(protein_md5),
                FOREIGN KEY (nt_md5) REFERENCES NUCLEOTIDE_SEQUENCE(nt_md5)
            );
        """)
    except sqlite3.Error as e:
        print(f"Error creating tables in the sequence database:\n{e}")
        sys.exit(1)


if __name__ == "__main__":
    """
    Sys.argv:
    [1] Path to local sequence db
    [2] Path to fasta file
    [3] Method to run
    [4] Pipeline input is nucleic sequences
    """
    if len(sys.argv) != 5:
        print(f"Expected 4 arguments but received {len(sys.argv) - 1}.")
        sys.exit(1)
    db_path, fasta, method, nucleic = sys.argv[1:5]
    nucleic = nucleic.lower() == "true"

    if method == "populate_sequences":
        populate_sequences(db_path, fasta, nucleic)
    elif method == "insert_translated_seqs":
        insert_translated_seqs(db_path, fasta)
    else:
        print(f"Method not recognised: {method}")
        sys.exit(1)
