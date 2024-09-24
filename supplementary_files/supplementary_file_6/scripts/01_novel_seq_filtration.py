from Bio import SeqIO
import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
from pathlib import Path
import hashlib

#Pass user specifiedarguments to variables
novel_sequences = sys.argv[1]
current_marker_alleles = sys.argv[2]
save_new_marker_alleles = sys.argv[3]
save_subsequences_to = sys.argv[4]

#Extracting sequences of interest 
def seq_extraction(fastafile, key_search):
    """Extract sequences from a fasta file if
    key_search found in header. 
    Sequences with ambiguity symbols are discarded. 

    :param fastafile: multi FASTA file
    :param key_search: str
    :return: list of SeqRecords
    """
    seqlist = list(SeqIO.parse(fastafile, "fasta"))
    found_seq = []

    for seq_record in seqlist:
        if key_search in seq_record.description:
            if re.search(r"[^ATGC]", str(seq_record.seq.upper())):
                print(f"Ambiguity symbol found in {seq_record.description.split('/')} for {key_search}")
            else:
                found_seq.append(seq_record)
    

    return list(found_seq)




rRNA = seq_extraction(f"{novel_sequences}", "16S")
atpD = seq_extraction(f"{novel_sequences}", "atpD")
gyrB = seq_extraction(f"{novel_sequences}", "gyrB")
recA = seq_extraction(f"{novel_sequences}", "recA")
rpoB = seq_extraction(f"{novel_sequences}", "rpoB")
trpB = seq_extraction(f"{novel_sequences}", "trpB")




#Delete novel sequences if they already exist in the database
def remove_duplicates(novel_seq, gene):
    """Remove novel sequences when already present in the database.

    :param novel_seq: a SeqIO objects with novel sequences
    :param database: a SeqIO objectes with sequences already in db
    """

    pubMLST_alleles = list(SeqIO.parse(f'{current_marker_alleles}/{gene}.tfa', "fasta"))

    remove = []
    novel = []
    for records in pubMLST_alleles:
        for record in novel_seq:
            if record.seq.upper() == records.seq:
                print(f'{record.description} and {records.description} are identical')
                remove.append(record.description)
                
    for records in novel_seq:
        if records.description not in remove:
            novel.append(records)
            
    return novel


rRNA = remove_duplicates(rRNA, '16S')
atpD = remove_duplicates(atpD, 'atpD')
gyrB = remove_duplicates(gyrB, 'gyrB')
recA = remove_duplicates(recA, 'recA')
rpoB = remove_duplicates(rpoB, 'rpoB')
trpB = remove_duplicates(trpB, 'trpB')


#delete novel sequences if longer varainats are found!

def remove_substrings_novel(novel):
    """Remove shorter novel subsequences, when longer novel
    variant has been found by mlst tool. 

    :param novel: a SeqIO objects with novel sequences
    """

    try:
        remove = []
        novel_no_substr = []
        for index, seq in enumerate(novel):
            for seq2 in novel[index+1:]:
                if seq.seq in seq2.seq.upper():
                    print(f"Longer variant of {seq.description} has been found!")
                    remove.append(seq.description)
                if seq2.seq in seq.seq:
                    print(f"Longer variant of {seq2.description} has been found!")
                    remove.append(seq2.description)

        for records in novel:
            if records.description not in remove:
                novel_no_substr.append(records)
        

    except IndexError:
        print()

    return novel_no_substr

rRNA = remove_substrings_novel(rRNA)
atpD = remove_substrings_novel(atpD)
gyrB = remove_substrings_novel(gyrB)
recA = remove_substrings_novel(recA)
rpoB = remove_substrings_novel(rpoB)
trpB = remove_substrings_novel(trpB)



def remove_substrings_db(novel_seq, gene):
    """Remove novel sequences that are substings of sequences in db
    or when records in db are substrings of a novel allele.

    :param novel: a SeqIO objects with novel sequences
    :param database: a SeqIO objectes with sequences already in db
    """

    pubMLST_alleles = list(SeqIO.parse(f'{current_marker_alleles}/{gene}.tfa', "fasta"))

    remove = []
    novel_db = []
    subsequences = []

    for pubmlst_record in pubMLST_alleles:
        for novel_record in list(novel_seq):
            if pubmlst_record.seq in novel_record.seq.upper():
                print(f'{pubmlst_record.description} is a subsequence of {novel_record.description}!')
                subsequences.append(novel_record)
                subsequences.append(pubmlst_record)
                remove.append(novel_record.description)
            if novel_record.seq.upper() in pubmlst_record.seq:
                print(f'{novel_record.description} is a subsequene of {pubmlst_record.description}!')
                subsequences.append(novel_record)
                subsequences.append(pubmlst_record)
                remove.append(novel_record.description)
    for records in list(novel_seq):
        if records.description not in remove:
            novel_db.append(records)

    if len(subsequences) != 0:
        SeqIO.write(subsequences, f'{save_subsequences_to}/{gene}_subsequences.fasta', 'fasta')

    return novel_db


rRNA = remove_substrings_db(rRNA, '16S')
atpD = remove_substrings_db(atpD, 'atpD')
gyrB = remove_substrings_db(gyrB, 'gyrB')
recA = remove_substrings_db(recA, 'recA')
rpoB = remove_substrings_db(rpoB, 'rpoB')
trpB = remove_substrings_db(trpB, 'trpB')



# Assign allele numbers to novel genes

def assign_allelic_num(fasta, gene):
    """Return fasta file with existing and novel alleles, 
    with MLST acceptable format 

    :param fasta: fasta files with novel alleles 
    """
    
    seqlist = fasta

    pubMLST_alleles = list(SeqIO.parse(Path(f'{current_marker_alleles}/{gene}.tfa').expanduser(), "fasta"))

    current_allele = len(pubMLST_alleles)


    seen = []

    for _ in seqlist:
        if hashlib.md5(str(_.seq).encode('utf-8')).hexdigest() not in seen:
            seen.append(hashlib.md5(str(_.seq).encode('utf-8')).hexdigest())
            current_allele +=1
            _.description = f'{gene}_{current_allele}'
            _.name = f'{gene}_{current_allele}'
            _.id = f'{gene}_{current_allele}'
            _.seq = _.seq.upper()
            pubMLST_alleles.append(_)

        else:
            print(f"{hashlib.md5(str(_.seq).encode('utf-8')).hexdigest()} already present!")


    updated = SeqIO.write(pubMLST_alleles, f'{save_new_marker_alleles}/{gene}.tfa', 'fasta')

    return pubMLST_alleles

assign_allelic_num(rRNA, "16S")
assign_allelic_num(atpD, "atpD")
assign_allelic_num(gyrB, "gyrB")
assign_allelic_num(recA, "recA")
assign_allelic_num(rpoB, "rpoB")
assign_allelic_num(trpB, "trpB")

