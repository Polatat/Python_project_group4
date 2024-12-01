from collections import Counter
from Bio import SeqIO
import re
import gzip

class FastqStat:
    def __init__(self, filename):
        self.filename = filename

    def read_fastq(self):
        """Reads the FASTQ file and returns a list of SeqRecord objects."""
        if self.filename.endswith('.gz'):
            with gzip.open(self.filename,"rt") as fastq_file:
                return list(SeqIO.parse(fastq_file,"fastq"))
        else:
            with open(self.filename, 'r') as fastq_file:
                return list(SeqIO.parse(fastq_file, "fastq"))

    def sequence(self):
        """Return a list of sequence IDs from the FASTQ file."""
        records = self.read_fastq()
        return [record.id for record in records]

    def calculate_sequence_lengths(self):
        """Calculate lengths of all sequences."""
        records = self.read_fastq()
        return [len(record.seq) for record in records]

    def calculate_gc_content_per_sequence(self):
        """Calculate GC content for each sequence."""
        gc_contents = []
        records = self.read_fastq()
        for record in records:
            gc_count = record.seq.count('G') + record.seq.count('C')
            gc_contents.append((gc_count / len(record.seq)) * 100 if len(record.seq) > 0 else 0)
        return gc_contents

    def calculate_mean_quality_scores(self):
        """Calculate mean quality scores for each sequence."""
        quality_scores = []
        records = self.read_fastq()
        for record in records:
            scores = record.letter_annotations["phred_quality"]
            mean_quality = sum(scores) / len(scores) if scores else 0
            quality_scores.append(mean_quality)
        return quality_scores

    def calculate_barcode_distribution(self, barcode_length, header_barcode=False):
        """Calculate the distribution of barcodes."""
        barcode_counter = Counter()
        records = self.read_fastq()
        for record in records:
            if header_barcode:
                match = re.search(r'barcode=(\S+)', record.description)
                if match:
                    barcode = match.group(1)
                else:
                    continue
            else:
                barcode = str(record.seq)[:barcode_length]  # Assuming barcode is at the start
            barcode_counter[barcode] += 1

        total_barcodes = sum(barcode_counter.values())
        distribution = {barcode: (count / total_barcodes * 100) for barcode, count in barcode_counter.items()} if total_barcodes > 0 else {}
        return distribution

    def calculate_metrics(self, barcode_length, header_barcode= False):
        """Calculate and return all metrics in a structured format."""
        sequence_ids = self.sequence()
        lengths = self.calculate_sequence_lengths()
        gc_contents = self.calculate_gc_content_per_sequence()
        mean_quality_scores = self.calculate_mean_quality_scores()
        barcode_distribution = self.calculate_barcode_distribution(barcode_length, header_barcode)

        metrics = []
        records = self.read_fastq()

        for idx, record in enumerate(records):
            # Extract the barcode from the description (e.g., 'barcode=barcode01')
            match = re.search(r'barcode=(\S+)', record.description)
            barcode_group = match.group(1) if match else None

            # # Extract the barcode from the sequence
            barcode = str(record.seq)[:barcode_length]  # Assuming the barcode is at the start of the sequence

            metrics.append({
                "sequence_id": sequence_ids[idx],
                "length": lengths[idx],
                "gc_content​ (%)": round(gc_contents[idx],2) if idx < len(gc_contents) else None,
                "mean_quality_score": round(mean_quality_scores[idx],2) if idx < len(mean_quality_scores) else None,
                "barcode": barcode,
                "barcode_group": barcode_group,
            })
        return metrics
    # Using in filtering
    def get_mean_quality_scores(self):
        return sum(self.calculate_mean_quality_scores()) /len(self.records)
    def get_gc_content_stats(self):
        gc_contents = self.calculate_gc_content_per_sequence()
        return {
              'min_gc': min(gc_contents) if gc_contents else 0,
              'max_gc': max(gc_contents) if gc_contents else 0,
              'avg_gc': sum(gc_contents)/ len(gc_contents) if gc_contents else 0
        }

# Usage Example
if __name__ == "__main__":
    mock_gene = "mock_gene_ex002.fastq"  # Replace with your FASTQ file path
    stat_cal = FastqStat(mock_gene)
    barcode_length = 6  # Specify the length of the barcode

    metrics = stat_cal.calculate_metrics(barcode_length)
    # Print the metrics in a structured format
    print(f"{'Seq ID':<45} {'Length':<15} {'GC Content (%)':<15} {'Mean Quality Score':<20} {'Barcode Group':<10} {'Dist(%)':<25}")
    for metric in metrics:
        print(f"{metric['sequence_id']:<45} {metric['length']:<15} {metric['gc_content']:<15.2f} {metric['mean_quality_score']:<20.2f} {metric['barcode_group']:<15} {metric['barcode_distribution']:<13.2f}")

# class BarcodeGroupSum:
#      def __init__(self, metrics):
#         self.metrics = group.
