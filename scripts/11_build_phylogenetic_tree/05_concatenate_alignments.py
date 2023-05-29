"""This script was used to concatenate CDS sequences, and generate a single 
aligmnents file that can be used as an output for RAXML-ng. Additionally, this
scipt was used to generate partitional file for Modeltest-ng. 

The partition file format is as follow:

DATA_TYPE, PARTITION_NAME = PARTITION_SITES

To run the script you need to pass the following arguments

`arg1` = path to alignments that need to be concatenated
`arg2` = FASTA output file path
`arg3` = partition output file path


We run:
`python 05_concatenate_alignments.py ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/build_tree/alignments/backthreaded_nt_alignments_and_16S ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/build_tree/alignments/concatenated_alignment/concantenated_alignment.fasta ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/build_tree/alignments/concatenated_alignment/partition_file.txt`

"""

#SetUp
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
from Bio import AlignIO

#pass given arguments to variables
input_alignments = sys.argv[1]
output_concat_alignment = sys.argv[2]
output_part_file = sys.argv[3]

#Path to SCO cds alignments
datadir = Path(f"{input_alignments}").expanduser()
filenames = sorted(datadir.glob("*"))


#Concatenate sequence and generate partition file
concatenated = {}
partitions = []
last_position = 0
for filename in filenames:
    if '.DS' not in str(filename):
        OG = str(filename).split('/')[-1].split('_')[0]
        records = list(SeqIO.parse(filename, "fasta"))
        records_lengths = [len(_.seq) for _ in records]
        next_position = last_position + records_lengths[0]
        partitions.append(f"DNA, {OG} = {last_position +1}-{next_position}")
        last_position = next_position
        for record in records:

            accession = record.description
            sequence = record.seq
            if accession not in concatenated:
                concatenated[accession] = sequence
            else:
                concatenated[accession] += sequence



#Save concatenate sequences
sequence_list = []
for key, value in concatenated.items():
    record = SeqRecord(
                seq = value, 
                description = key,
                id = key

            )
    sequence_list.append(record)
SeqIO.write(sequence_list, Path(f"{output_concat_alignment}").expanduser(), "fasta")

#Save partition file

with open(Path(f"{output_part_file}").expanduser(), "w") as text_file:
    text_file.write('\n'.join(partitions))