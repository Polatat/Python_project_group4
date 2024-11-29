from Bio import SeqIO
import gzip
import os

class FastqParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.grouped_sequences = {}
        

    def parse_fastq(self):
        if not os.path.isfile(self.file_path):
            raise FileNotFoundError("FASTQ file not found:{self.file_path}")
        
        if self.file_path.endswith('.gz'):
            handle = gzip.open(self.file_path, "rt")
        else:
            handle = open(self.file_path, "r")

        for record in SeqIO.parse(handle, "fastq"):
            
            group, barcode_seq = self.extract_barcode(record)
            
            if barcode_seq and group:
                combined_key = f"{barcode_seq}_{group}"
            else:
                combined_key = "No Barcode_No Group"

            if combined_key not in self.grouped_sequences:
                self.grouped_sequences[combined_key] = {
                    'barcode_seq': barcode_seq,
                    'group': group if group else "No Group",
                    'sequences': []
                }

            self.grouped_sequences[combined_key]['sequences'].append(str(record.seq))

        handle.close()

    def extract_barcode(self, record):
        
        barcode_seq = str(record.seq[:6])
        
        if 'barcode=' in record.description:
            group = record.description.split('barcode=')[1].split()[0]
        else:
            group = None 
        
        return group, barcode_seq

    def print_grouped_sequences(self):
        output = []
        for key, data in self.grouped_sequences.items():
            barcode_seq = data['barcode_seq']
            group = data['group']
            output.append(f"Barcode: {barcode_seq}, Number of Sequences: {len(data['sequences'])}, group={group}")
            for seq in data['sequences']:
                output.append(f"  Sequence: {seq}") 
        output.append("\n")
        return "\n".join(output)
    
    def get_grouped_sequences(self):
        return self.grouped_sequences



